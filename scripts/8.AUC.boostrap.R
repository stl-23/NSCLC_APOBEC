library(openxlsx)
library(tidyverse)
library(pROC)

setwd("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/")

# ============================================================================
# 1. 数据加载与准备
# ============================================================================

# Combine 26 samples
load("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/DNA_and_ctDNA_26samples/combine_group.rdata")

# Tissue group (34 samples)
tissue_group <- read.xlsx("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/DNA/01.baseline_DNA_oncoplot/clinical_dna.xlsx",
                          sheet = 1)
tissue_group$Signature_APOBEC_ml[is.na(tissue_group$Signature_APOBEC_ml)] <- 0
# 统一 response_new: response->MPR, non.response->NonMPR
tissue_group$response_new <- ifelse(tissue_group$response_new == "response", "MPR", "NonMPR")
dna_apobec <- tissue_group

# ctDNA group (35 samples)
load("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/ctDNA/01.PL0_oncoplot/1.ctDNA_oncoplot_rm_NSCLC/add_2_apobec_neg_samples/group_2025.5.19.rdata")
ctdna_apobec <- group
ctdna_apobec$patient_id <- as.character(ctdna_apobec$patient_id)
# 统一 response_new: Non.MPR->NonMPR
ctdna_apobec$response_new <- as.character(ctdna_apobec$response_new)
ctdna_apobec$response_new[ctdna_apobec$response_new == "Non.MPR"] <- "NonMPR"

# Combine 数据准备: merge TMB/PDL1
combine_group.order <- combine_group.order %>%
  dplyr::left_join(dna_apobec[, c("patient_id", "Signature_APOBEC_ml", "TMB", "PDL1")],
                   by = c("ID" = "patient_id"))
combine_group.order <- combine_group.order %>%
  dplyr::left_join(ctdna_apobec[, c("patient_id", "Signature_APOBEC_ml", "TMB", "PDL1.fraction")],
                   by = c("ID" = "patient_id"))

# 用列名重命名 (dplyr 自动加 .x/.y 后缀)
colnames(combine_group.order)[colnames(combine_group.order) == "Signature_APOBEC_ml.x"] <- "DNA.APOBEC"
colnames(combine_group.order)[colnames(combine_group.order) == "Signature_APOBEC_ml.y"] <- "ctDNA.APOBEC"
combine_group.order$DNA.APOBEC[is.na(combine_group.order$DNA.APOBEC)] <- 0
combine_group.order$combine_APOBEC <- rowMeans(combine_group.order[, c("DNA.APOBEC", "ctDNA.APOBEC")])

# TMB: join 已带来 TMB.x (tissue) 和 TMB.y (ctDNA), 统一用 tissue TMB
if ("TMB.x" %in% colnames(combine_group.order)) {
  combine_group.order$TMB <- combine_group.order$`TMB.x`
} else {
  combine_group.order$TMB <- tissue_group$TMB[match(combine_group.order$ID, tissue_group$patient_id)]
}

# PDL1 binary: join 已带来 PDL1 (tissue, 0/1)
combine_group.order$PDL1_binary <- combine_group.order$PDL1

# PDL1 fraction (from ctDNA), dichotomize at 1%
combine_group.order$PDL1_fraction <- combine_group.order$PDL1.fraction
combine_group.order$PDL1_ge1 <- ifelse(!is.na(combine_group.order$PDL1_fraction) &
                                        combine_group.order$PDL1_fraction >= 0.01, 1, 0)

# ============================================================================
# 2. Z-score 辅助函数 & 组合 biomarker
# ============================================================================
z_score <- function(x) {
  x_num <- as.numeric(x)
  (x_num - mean(x_num, na.rm = TRUE)) / sd(x_num, na.rm = TRUE)
}

z_mean <- function(...) {
  cols <- list(...)
  z_cols <- lapply(cols, z_score)
  rowMeans(do.call(cbind, z_cols), na.rm = TRUE)
}

# --- Tissue (34) 组合 ---
tissue_group$APOBEC_TMB   <- z_mean(tissue_group$Signature_APOBEC_ml, tissue_group$TMB)
tissue_group$APOBEC_PDL1  <- z_mean(tissue_group$Signature_APOBEC_ml, tissue_group$PDL1)
tissue_group$Triplet      <- z_mean(tissue_group$Signature_APOBEC_ml, tissue_group$TMB, tissue_group$PDL1)

# --- ctDNA (35) 组合 ---
ctdna_apobec$APOBEC_TMB   <- z_mean(ctdna_apobec$Signature_APOBEC_ml, ctdna_apobec$TMB)
ctdna_apobec$APOBEC_PDL1  <- z_mean(ctdna_apobec$Signature_APOBEC_ml, ctdna_apobec$PDL1.fraction)
ctdna_apobec$Triplet      <- z_mean(ctdna_apobec$Signature_APOBEC_ml, ctdna_apobec$TMB, ctdna_apobec$PDL1.fraction)

# --- Combine (26) 组合 ---
combine_group.order$APOBEC_TMB   <- z_mean(combine_group.order$combine_APOBEC, combine_group.order$TMB)
combine_group.order$APOBEC_PDL1  <- z_mean(combine_group.order$combine_APOBEC, combine_group.order$PDL1_ge1)
combine_group.order$Triplet      <- z_mean(combine_group.order$combine_APOBEC, combine_group.order$TMB, combine_group.order$PDL1_ge1)

# ============================================================================
# 3. 分析参数
# ============================================================================
pos_level <- "MPR"
neg_level <- "NonMPR"
n_boot    <- 1000
set.seed(42)

# ============================================================================
# 4. 辅助函数
# ============================================================================
calc_metrics <- function(truth, predicted, pos_label, neg_label) {
  tp <- sum(truth == pos_label & predicted == pos_label, na.rm = TRUE)
  fp <- sum(truth == neg_label & predicted == pos_label, na.rm = TRUE)
  fn <- sum(truth == pos_label & predicted == neg_label, na.rm = TRUE)
  tn <- sum(truth == neg_label & predicted == neg_label, na.rm = TRUE)
  c(Sensitivity = ifelse(tp + fn > 0, tp / (tp + fn), NA),
    Specificity = ifelse(tn + fp > 0, tn / (tn + fp), NA),
    PPV         = ifelse(tp + fp > 0, tp / (tp + fp), NA),
    NPV         = ifelse(tn + fn > 0, tn / (tn + fn), NA))
}

safe_quantile <- function(x, probs) {
  x_valid <- x[!is.na(x)]
  if (length(x_valid) < 2) return(rep(NA, length(probs)))
  quantile(x_valid, probs = probs, na.rm = TRUE)
}

# 从 bootstrap 矩阵 (n_boot x 4) 计算 4x2 CI 矩阵
compute_metrics_ci <- function(boot_matrix) {
  metric_names <- colnames(boot_matrix)
  ci <- t(vapply(seq_along(metric_names), function(j) {
    safe_quantile(boot_matrix[, j], c(0.025, 0.975))
  }, numeric(2)))
  rownames(ci) <- metric_names
  colnames(ci) <- c("2.5%", "97.5%")
  ci[ci < 0] <- 0; ci[ci > 1] <- 1
  ci
}

# ============================================================================
# 5. Bootstrap 分析函数
#    - 连续 predictor: AUC + optimism-corrected metrics (Youden cutoff)
#    - 二分类 predictor: 无 AUC, 直接以 predictor 值作为分类结果
# ============================================================================
bootstrap_analysis <- function(data, predictor, response, pos_level, neg_level,
                               n_boot = 1000, is_binary = FALSE) {
  # 去除 predictor 缺失
  resp_vec <- data[[response]]
  pred_vec <- data[[predictor]]
  complete <- !is.na(resp_vec) & !is.na(pred_vec)
  data <- data[complete, ]
  n <- nrow(data)
  truth <- data[[response]]
  pred  <- data[[predictor]]

  if (is_binary) {
    # ===== 二分类 predictor =====
    pred_class <- ifelse(pred >= 0.5, pos_level, neg_level)
    metrics_apparent <- calc_metrics(truth, pred_class, pos_level, neg_level)

    metrics_boot <- matrix(NA, nrow = n_boot, ncol = 4,
                           dimnames = list(NULL, c("Sensitivity", "Specificity", "PPV", "NPV")))
    for (i in seq_len(n_boot)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      metrics_boot[i, ] <- calc_metrics(truth[idx], pred_class[idx], pos_level, neg_level)
    }
    metrics_ci <- compute_metrics_ci(metrics_boot)

    return(list(
      auc_apparent = NA, auc_ci = c(NA, NA),
      metrics_apparent = metrics_apparent, metrics_corrected = metrics_apparent,
      metrics_ci = metrics_ci, optimism_metrics = rep(0, 4),
      best_threshold = 0.5, threshold_ci = c(NA, NA),
      n_samples = n, is_binary = TRUE
    ))

  } else {
    # ===== 连续 predictor =====
    roc_apparent <- roc(response = truth, predictor = pred,
                        levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
    auc_apparent <- as.numeric(auc(roc_apparent))
    best_thresh <- coords(roc_apparent, "best", best.method = "youden",
                          ret = "threshold", drop = TRUE)
    pred_class_apparent <- ifelse(pred >= best_thresh, pos_level, neg_level)
    metrics_apparent <- calc_metrics(truth, pred_class_apparent, pos_level, neg_level)

    auc_boot_each <- numeric(n_boot)
    boot_thresholds <- numeric(n_boot)
    optimism_metrics <- matrix(NA, nrow = n_boot, ncol = 4,
                               dimnames = list(NULL, c("Sensitivity", "Specificity", "PPV", "NPV")))
    metrics_corrected_each <- matrix(NA, nrow = n_boot, ncol = 4,
                                     dimnames = list(NULL, c("Sensitivity", "Specificity", "PPV", "NPV")))

    for (i in seq_len(n_boot)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      boot_data <- data[idx, ]
      boot_truth <- boot_data[[response]]; boot_pred <- boot_data[[predictor]]

      roc_boot <- roc(response = boot_truth, predictor = boot_pred,
                      levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
      auc_boot_each[i] <- as.numeric(auc(roc_boot))

      boot_thresh <- coords(roc_boot, "best", best.method = "youden",
                            ret = "threshold", drop = TRUE)
      boot_thresholds[i] <- boot_thresh
      pred_class_boot <- ifelse(boot_pred >= boot_thresh, pos_level, neg_level)
      metrics_boot <- calc_metrics(boot_truth, pred_class_boot, pos_level, neg_level)

      pred_class_orig <- ifelse(pred >= boot_thresh, pos_level, neg_level)
      metrics_orig <- calc_metrics(truth, pred_class_orig, pos_level, neg_level)

      optimism_metrics[i, ] <- metrics_boot - metrics_orig
      metrics_corrected_each[i, ] <- metrics_apparent - optimism_metrics[i, ]
    }

    auc_ci <- safe_quantile(auc_boot_each, c(0.025, 0.975))
    auc_ci <- pmin(pmax(auc_ci, 0), 1)

    threshold_ci <- safe_quantile(boot_thresholds, c(0.025, 0.975))

    mean_optimism <- colMeans(optimism_metrics, na.rm = TRUE)
    metrics_corrected <- metrics_apparent - mean_optimism

    metrics_ci <- compute_metrics_ci(metrics_corrected_each)

    return(list(
      auc_apparent = auc_apparent, auc_ci = auc_ci,
      metrics_apparent = metrics_apparent, metrics_corrected = metrics_corrected,
      metrics_ci = metrics_ci, optimism_metrics = mean_optimism,
      best_threshold = best_thresh, threshold_ci = threshold_ci,
      n_samples = n, is_binary = FALSE
    ))
  }
}

# ============================================================================
# 6. 定义全部分析任务
# ============================================================================
analysis_tasks <- list(
  # --- Tissue (34 samples) ---
  list(label="Tissue_APOBEC",     dataset="Tissue",  marker="APOBEC",        data=tissue_group,   pred="Signature_APOBEC_ml", binary=FALSE),
  list(label="Tissue_TMB",        dataset="Tissue",  marker="TMB",           data=tissue_group,   pred="TMB",        binary=FALSE),
  list(label="Tissue_PDL1",       dataset="Tissue",  marker="PDL1",          data=tissue_group,   pred="PDL1",       binary=TRUE),
  list(label="Tissue_APOBEC_TMB", dataset="Tissue",  marker="APOBEC+TMB",    data=tissue_group,   pred="APOBEC_TMB", binary=FALSE),
  list(label="Tissue_APOBEC_PDL1",dataset="Tissue",  marker="APOBEC+PDL1",   data=tissue_group,   pred="APOBEC_PDL1",binary=FALSE),
  list(label="Tissue_Triplet",    dataset="Tissue",  marker="APOBEC+TMB+PDL1",data=tissue_group,  pred="Triplet",    binary=FALSE),
  # --- ctDNA (35 samples) ---
  list(label="ctDNA_APOBEC",      dataset="ctDNA",   marker="APOBEC",        data=ctdna_apobec,   pred="Signature_APOBEC_ml", binary=FALSE),
  list(label="ctDNA_TMB",         dataset="ctDNA",   marker="TMB",           data=ctdna_apobec,   pred="TMB",            binary=FALSE),
  list(label="ctDNA_PDL1",        dataset="ctDNA",   marker="PDL1",          data=ctdna_apobec,   pred="PDL1.fraction",   binary=FALSE),
  list(label="ctDNA_APOBEC_TMB",  dataset="ctDNA",   marker="APOBEC+TMB",    data=ctdna_apobec,   pred="APOBEC_TMB",      binary=FALSE),
  list(label="ctDNA_APOBEC_PDL1", dataset="ctDNA",   marker="APOBEC+PDL1",   data=ctdna_apobec,   pred="APOBEC_PDL1",     binary=FALSE),
  list(label="ctDNA_Triplet",     dataset="ctDNA",   marker="APOBEC+TMB+PDL1",data=ctdna_apobec,  pred="Triplet",         binary=FALSE),
  # --- Combine (26 samples) ---
  list(label="Combine_APOBEC",     dataset="Combine", marker="APOBEC",        data=combine_group.order, pred="combine_APOBEC", binary=FALSE),
  list(label="Combine_TMB",         dataset="Combine", marker="TMB",           data=combine_group.order, pred="TMB",            binary=FALSE),
  list(label="Combine_PDL1",        dataset="Combine", marker="PDL1",          data=combine_group.order, pred="PDL1_binary",    binary=TRUE),
  list(label="Combine_APOBEC_TMB",  dataset="Combine", marker="APOBEC+TMB",    data=combine_group.order, pred="APOBEC_TMB",     binary=FALSE),
  list(label="Combine_APOBEC_PDL1", dataset="Combine", marker="APOBEC+PDL1",   data=combine_group.order, pred="APOBEC_PDL1",    binary=FALSE),
  list(label="Combine_Triplet",     dataset="Combine", marker="APOBEC+TMB+PDL1",data=combine_group.order, pred="Triplet",        binary=FALSE)
)

# ============================================================================
# 7. 执行全部分析
# ============================================================================
results <- list()

for (i in seq_along(analysis_tasks)) {
  task <- analysis_tasks[[i]]
  cat(sprintf("\n========== %s / %s (n=%d) ==========\n",
              task$dataset, task$marker, nrow(task$data)))

  res <- bootstrap_analysis(
    data       = task$data,
    predictor  = task$pred,
    response   = "response_new",
    pos_level  = pos_level,
    neg_level  = neg_level,
    n_boot     = n_boot,
    is_binary  = task$binary
  )
  results[[task$label]] <- res

  if (res$is_binary) {
    cat(sprintf("  [Binary predictor, no AUC]\n"))
  } else {
    cat(sprintf("  AUC:               %.4f\n", res$auc_apparent))
    cat(sprintf("  95%% CI (bootstrap): [%.4f, %.4f]\n", res$auc_ci[1], res$auc_ci[2]))
  }
  cat(sprintf("  Best threshold:    %.4f", res$best_threshold))
  if (!is.na(res$threshold_ci[1])) {
    cat(sprintf("  (95%% CI: [%.4f, %.4f])", res$threshold_ci[1], res$threshold_ci[2]))
  }
  cat("\n")
  cat("  --- Classification metrics ---\n")
  for (m in c("Sensitivity", "Specificity", "PPV", "NPV")) {
    if (res$is_binary) {
      cat(sprintf("    %s: %.4f  (95%% CI: [%.4f, %.4f])\n",
                  m, res$metrics_corrected[m], res$metrics_ci[m, 1], res$metrics_ci[m, 2]))
    } else {
      cat(sprintf("    %s: %.4f  (optimism=%.4f, 95%% CI: [%.4f, %.4f])\n",
                  m, res$metrics_corrected[m], res$optimism_metrics[m],
                  res$metrics_ci[m, 1], res$metrics_ci[m, 2]))
    }
  }
}

# ============================================================================
# 8. 汇总结果表
# ============================================================================
summary_rows <- lapply(analysis_tasks, function(task) {
  res <- results[[task$label]]
  data.frame(
    Dataset     = task$dataset,
    N_samples   = res$n_samples,
    Marker      = task$marker,
    Predictor   = task$pred,
    Binary      = res$is_binary,
    AUC         = ifelse(res$is_binary, NA, round(res$auc_apparent, 4)),
    AUC_CI      = ifelse(res$is_binary, NA,
                         sprintf("[%.4f, %.4f]", res$auc_ci[1], res$auc_ci[2])),
    Sensitivity = round(res$metrics_corrected["Sensitivity"], 4),
    Sens_CI     = sprintf("[%.4f, %.4f]", res$metrics_ci["Sensitivity", 1], res$metrics_ci["Sensitivity", 2]),
    Specificity = round(res$metrics_corrected["Specificity"], 4),
    Spec_CI     = sprintf("[%.4f, %.4f]", res$metrics_ci["Specificity", 1], res$metrics_ci["Specificity", 2]),
    PPV         = round(res$metrics_corrected["PPV"], 4),
    PPV_CI      = sprintf("[%.4f, %.4f]", res$metrics_ci["PPV", 1], res$metrics_ci["PPV", 2]),
    NPV         = round(res$metrics_corrected["NPV"], 4),
    NPV_CI      = sprintf("[%.4f, %.4f]", res$metrics_ci["NPV", 1], res$metrics_ci["NPV", 2]),
    Threshold   = round(res$best_threshold, 4),
    Threshold_CI = ifelse(is.na(res$threshold_ci[1]), NA,
                          sprintf("[%.4f, %.4f]", res$threshold_ci[1], res$threshold_ci[2])),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)

cat("\n\n========== SUMMARY TABLE ==========\n")
print(summary_table, right = FALSE, row.names = FALSE)

write.xlsx(summary_table, file = "Bootstrap_ROC_results_summary.xlsx", overwrite = TRUE)

# ============================================================================
# 9. 绘制 ROC 曲线 (按 dataset 分组, 仅连续 predictor)
# ============================================================================
colors5 <- c("APOBEC" = "#333333",
             "TMB" = "#9C27B0", "PDL1" = "#FF9800",
             "APOBEC+TMB" = "#E53935", "APOBEC+PDL1" = "#1E88E5",
             "APOBEC+TMB+PDL1" = "#43A047")

for (ds in c("Tissue", "ctDNA", "Combine")) {
  ds_tasks <- Filter(function(t) t$dataset == ds && !t$binary, analysis_tasks)
  if (length(ds_tasks) == 0) next

  pdf(sprintf("ROC_curves_%s_all_markers.pdf", ds), width = 6, height = 5, family = "ArialMT")
  plot(NULL, xlim = c(1, 0), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = sprintf("%s: ROC Curves (All Markers)", ds))
  abline(a = 0, b = 1, lty = 2, col = "gray60")

  y_offset <- 0
  for (tk in ds_tasks) {
    df <- tk$data
    resp <- df[["response_new"]]
    pred_vals <- df[[tk$pred]]
    roc_obj <- roc(response = resp, predictor = pred_vals,
                   levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
    plot.roc(roc_obj, add = TRUE,
             col = colors5[tk$marker], lwd = 2,
             legacy.axes = TRUE,
             print.auc = TRUE, print.auc.y = 0.1 + y_offset * 0.08)
    y_offset <- y_offset + 1
  }

  legend_labels <- sapply(ds_tasks, function(t) {
    res <- results[[t$label]]
    sprintf("%s (AUC=%.3f)", t$marker, res$auc_apparent)
  })
  legend("bottomright", legend = legend_labels,
         col = unname(colors5[names(ds_tasks)]),
         lwd = 2, bty = "n", cex = 0.7)
  dev.off()
}

# 单独 ROC 图
for (tk in analysis_tasks) {
  if (tk$binary) next
  res <- results[[tk$label]]
  df <- tk$data
  resp <- df[["response_new"]]
  pred_vals <- df[[tk$pred]]
  roc_obj <- roc(response = resp, predictor = pred_vals,
                 levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
  pdf(sprintf("ROC_%s_%s.pdf", tk$dataset, gsub("\\+", "plus", tk$marker)),
      width = 4, height = 4, family = "ArialMT")
  plot.roc(roc_obj,
           legacy.axes = TRUE, print.auc = TRUE,
           col = colors5[tk$marker],
           print.thres = "best", print.thres.best.method = "youden",
           print.thres.col = colors5[tk$marker],
           main = sprintf("%s %s (AUC=%.3f, n=%d)",
                          tk$dataset, tk$marker, res$auc_apparent, res$n_samples))
  dev.off()
}

# ============================================================================
# 10. APOBEC 连续变量 Logistic 回归: 预测概率以 0.5 阈值分类 & AUC 95% CI
#     评价 APOBEC 作为连续变量通过 logistic 回归预测 MPR 的能力
#     以预测概率 0.5 为阈值, 比较预测分类与实际 MPR/NonMPR 的一致性
# ============================================================================

logistic_auc_analysis <- function(data, apobec_col, response_col, dataset_name,
                                  n_boot = 1000) {
  cat(sprintf("\n========== Logistic-Predicted AUC: %s ==========\n", dataset_name))

  # 去除缺失
  df <- data[, c(apobec_col, response_col)]
  colnames(df) <- c("APOBEC", "response")
  df <- df[complete.cases(df), ]
  n <- nrow(df)

  # 二值化 response: MPR=1, NonMPR=0
  df$MPR <- ifelse(df$response == pos_level, 1, 0)

  # --- 拟合 logistic 回归 ---
  fit <- glm(MPR ~ APOBEC, data = df, family = binomial())
  cat(sprintf("  Logistic model: MPR ~ APOBEC  (n=%d)\n", n))
  cat(sprintf("  Coefficient: %.4f, OR=%.4f, P=%.4f\n",
              coef(fit)["APOBEC"],
              exp(coef(fit)["APOBEC"]),
              summary(fit)$coefficients["APOBEC", "Pr(>|z|)"]))

  # 预测概率
  df$pred_prob <- predict(fit, type = "response")

  # 以 0.5 为阈值分类
  df$pred_class <- ifelse(df$pred_prob >= 0.5, "MPR", "NonMPR")
  df$pred_class <- factor(df$pred_class, levels = c(neg_level, pos_level))

  # --- 表观分类指标 ---
  truth <- factor(df$response, levels = c(neg_level, pos_level))
  metrics_apparent <- calc_metrics(truth, df$pred_class, pos_level, neg_level)

  # 混淆矩阵
  TP <- sum(truth == pos_level & df$pred_class == pos_level)
  FP <- sum(truth == neg_level & df$pred_class == pos_level)
  FN <- sum(truth == pos_level & df$pred_class == neg_level)
  TN <- sum(truth == neg_level & df$pred_class == neg_level)
  cat(sprintf("  Confusion (threshold=0.5): TP=%d, FP=%d, FN=%d, TN=%d\n", TP, FP, FN, TN))
  cat(sprintf("  Sensitivity=%.4f, Specificity=%.4f, PPV=%.4f, NPV=%.4f\n",
              metrics_apparent["Sensitivity"], metrics_apparent["Specificity"],
              metrics_apparent["PPV"], metrics_apparent["NPV"]))

  # --- AUC (以预测概率为连续 predictor) ---
  roc_apparent <- roc(response = truth, predictor = df$pred_prob,
                      levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
  auc_apparent <- as.numeric(auc(roc_apparent))
  cat(sprintf("  AUC (apparent): %.4f\n", auc_apparent))

  # --- Bootstrap 95% CI for AUC ---
  auc_boot_each <- numeric(n_boot)
  metrics_boot <- matrix(NA, nrow = n_boot, ncol = 4,
                         dimnames = list(NULL, c("Sensitivity", "Specificity", "PPV", "NPV")))

  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_df <- df[idx, ]
    boot_truth <- factor(boot_df$response, levels = c(neg_level, pos_level))

    # bootstrap AUC (预测概率直接重采样)
    roc_boot <- roc(response = boot_truth, predictor = boot_df$pred_prob,
                    levels = c(neg_level, pos_level), direction = "<", quiet = TRUE)
    auc_boot_each[i] <- as.numeric(auc(roc_boot))

    # bootstrap 分类指标 (固定 0.5 阈值)
    boot_pred_class <- ifelse(boot_df$pred_prob >= 0.5, pos_level, neg_level)
    metrics_boot[i, ] <- calc_metrics(boot_truth, boot_pred_class, pos_level, neg_level)
  }

  auc_ci <- safe_quantile(auc_boot_each, c(0.025, 0.975))
  auc_ci <- pmin(pmax(auc_ci, 0), 1)
  cat(sprintf("  95%% CI (bootstrap): [%.4f, %.4f]\n", auc_ci[1], auc_ci[2]))

  # 分类指标 95% CI
  metrics_ci <- compute_metrics_ci(metrics_boot)
  for (m in c("Sensitivity", "Specificity", "PPV", "NPV")) {
    cat(sprintf("    %s: %.4f  (95%% CI: [%.4f, %.4f])\n",
                m, metrics_apparent[m], metrics_ci[m, 1], metrics_ci[m, 2]))
  }

  # --- 汇总 ---
  summary_df <- data.frame(
    Dataset          = dataset_name,
    N_samples        = n,
    APOBEC_OR        = exp(coef(fit)["APOBEC"]),
    APOBEC_P         = summary(fit)$coefficients["APOBEC", "Pr(>|z|)"],
    AUC              = round(auc_apparent, 4),
    AUC_CI           = sprintf("[%.4f, %.4f]", auc_ci[1], auc_ci[2]),
    TP               = TP,
    FP               = FP,
    FN               = FN,
    TN               = TN,
    Sensitivity      = round(metrics_apparent["Sensitivity"], 4),
    Sens_CI          = sprintf("[%.4f, %.4f]", metrics_ci["Sensitivity", 1], metrics_ci["Sensitivity", 2]),
    Specificity      = round(metrics_apparent["Specificity"], 4),
    Spec_CI          = sprintf("[%.4f, %.4f]", metrics_ci["Specificity", 1], metrics_ci["Specificity", 2]),
    PPV              = round(metrics_apparent["PPV"], 4),
    PPV_CI           = sprintf("[%.4f, %.4f]", metrics_ci["PPV", 1], metrics_ci["PPV", 2]),
    NPV              = round(metrics_apparent["NPV"], 4),
    NPV_CI           = sprintf("[%.4f, %.4f]", metrics_ci["NPV", 1], metrics_ci["NPV", 2]),
    Threshold        = 0.5,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  # 样本级别明细
  detail_df <- data.frame(
    Dataset       = dataset_name,
    SampleID      = if ("patient_id" %in% colnames(data)) data[complete.cases(data[, c(apobec_col, response_col)]), "patient_id"]
                    else if ("ID" %in% colnames(data)) data[complete.cases(data[, c(apobec_col, response_col)]), "ID"]
                    else seq_len(n),
    Actual        = df$response,
    APOBEC        = df$APOBEC,
    PredProb_MPR  = round(df$pred_prob, 4),
    PredClass_0.5 = df$pred_class,
    Concordant    = df$pred_class == df$response,
    stringsAsFactors = FALSE
  )

  list(summary = summary_df, detail = detail_df,
       auc_apparent = auc_apparent, auc_ci = auc_ci,
       metrics_apparent = metrics_apparent, metrics_ci = metrics_ci)
}

# --- 执行三组分析 ---
logistic_results <- list()

logistic_results[["Tissue"]] <- logistic_auc_analysis(
  tissue_group, "Signature_APOBEC_ml", "response_new", "Tissue", n_boot)

logistic_results[["ctDNA"]] <- logistic_auc_analysis(
  ctdna_apobec, "Signature_APOBEC_ml", "response_new", "ctDNA", n_boot)

logistic_results[["Combine"]] <- logistic_auc_analysis(
  combine_group.order, "combine_APOBEC", "response_new", "Combine", n_boot)

# --- 汇总表格 & 输出 Excel ---
logistic_summary_table <- do.call(rbind, lapply(logistic_results, function(x) x$summary))
logistic_detail_table  <- do.call(rbind, lapply(logistic_results, function(x) x$detail))

cat("\n\n========== LOGISTIC-PREDICTED AUC SUMMARY ==========\n")
print(logistic_summary_table, right = FALSE, row.names = FALSE)

write.xlsx(list(
  Summary        = logistic_summary_table,
  Sample_Detail  = logistic_detail_table
), file = "APOBEC_logistic_predicted_AUC.xlsx", overwrite = TRUE)

cat("\nLogistic predicted AUC results saved to APOBEC_logistic_predicted_AUC.xlsx\n")

cat("\nAll analyses completed. Results saved.\n")
