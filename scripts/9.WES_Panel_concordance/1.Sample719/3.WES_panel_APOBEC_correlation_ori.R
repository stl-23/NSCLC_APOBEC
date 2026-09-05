#!/usr/bin/env Rscript
# ============================================================================
# WES vs Panel APOBEC 突变特征一致性分析
# 比较 3 种 WES 方法 (deconstructSigs / MutationalPatterns / SigProfiler)
# 与 Panel (SigMA) 的 APOBEC 信号一致性
# ============================================================================

library(openxlsx)
library(corrplot)
library(ggplot2)

base_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation"
out_dir <- file.path(base_dir, "corr_gt5")
dir.create(out_dir, showWarnings = FALSE)

# ============================================================================
# 辅助函数
# ============================================================================

# 从各种样本名中提取 D 编号 (统一为 "D123" 格式)
# 支持: MP-D123, D123, D123_ND123.xxx, MPx-...-MP-D123-... 等
extract_dnum <- function(x) {
  out <- rep(NA_character_, length(x))
  # 优先匹配 MP-D[0-9]+
  m <- regexpr("MP-D[0-9]+", x, perl = TRUE)
  hit <- m > 0
  if (any(hit)) {
    out[hit] <- substring(x[hit], m[hit], m[hit] + attr(m, "match.length")[hit] - 1)
  }
  # 对未匹配的，尝试 D[0-9]+
  miss <- which(!hit)
  if (length(miss) > 0) {
    m2 <- regexpr("D[0-9]+", x[miss], perl = TRUE)
    hit2 <- m2 > 0
    if (any(hit2)) {
      out[miss[hit2]] <- substring(x[miss][hit2], m2[hit2], m2[hit2] + attr(m2, "match.length")[hit2] - 1)
    }
  }
  # 归一化: MP-D116 -> D116
  out <- sub("^MP-", "", out)
  return(out)
}

# ============================================================================
# 0. WES 样本去重 (D 格式 vs MP 格式, 保留 Total.SNV 最大的)
# ============================================================================
cat("===== 0. WES 样本去重 =====\n")
wes_stats <- read.csv("/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES_filter_statistics.csv")
wes_stats$dnum <- extract_dnum(wes_stats$Sample)
cat(sprintf("  WES 总样本: %d, 有效 D 编号: %d\n", nrow(wes_stats), sum(!is.na(wes_stats$dnum))))

# 按 D 编号去重: 保留 Total.SNV 最大的
wes_stats <- wes_stats[!is.na(wes_stats$dnum), ]
wes_dedup <- do.call(rbind, lapply(split(wes_stats, wes_stats$dnum), function(g) {
  g[which.max(g$Total.SNV), ]
}))
cat(sprintf("  去重后: %d 样本 (去除 %d 重复)\n", nrow(wes_dedup), nrow(wes_stats) - nrow(wes_dedup)))

# 记录被保留的样本名 (用于后续 WES 数据过滤)
kept_samples <- setNames(wes_dedup$Sample, wes_dedup$dnum)

# 展示几个去重例子
cat("  去重示例:\n")
for (dn in head(unique(wes_stats$dnum[duplicated(wes_stats$dnum)]), 3)) {
  g <- wes_stats[wes_stats$dnum == dn, c("Sample", "Total.SNV")]
  kept <- kept_samples[dn]
  cat(sprintf("    %s: 保留 [%s] (SNV=%d), 丢弃其他\n", dn, kept, max(g$Total.SNV)))
}

# ============================================================================
# 1. 读取 Panel SigMA APOBEC 值
# ============================================================================
cat("\n===== 1. 读取 Panel SigMA 结果 =====\n")
sig <- read.xlsx(file.path(base_dir, "MasterPanel/sigMA_signature_lite.xlsx"))
sig$dnum <- extract_dnum(sig$tumor)
cat(sprintf("  SigMA 样本数: %d, 有 D 编号: %d\n", nrow(sig), sum(!is.na(sig$dnum))))

# ============================================================================
# 2. 读取 QC 并过滤 (Coverage >= 99%, Duplication <= 0.7)
# ============================================================================
cat("\n===== 2. QC 过滤 =====\n")
qc <- read.xlsx(file.path(base_dir, "MasterPanel/QC_dedup.xlsx"))
qc$Sample <- gsub('xml:space="preserve">', '', qc$Sample, fixed = TRUE)
qc$dnum <- extract_dnum(qc$Sample)
qc$Coverage    <- as.numeric(qc$Coverage)
qc$Duplication <- as.numeric(qc$Duplication)

cat(sprintf("  QC 总样本: %d\n", nrow(qc)))
cat(sprintf("  Coverage < 99%%: %d\n", sum(qc$Coverage < 0.99, na.rm = TRUE)))
cat(sprintf("  Duplication > 0.7: %d\n", sum(qc$Duplication > 0.7, na.rm = TRUE)))

# Panel QC 过滤后
qc_pass <- qc$dnum[qc$Coverage >= 0.99 & qc$Duplication <= 0.7]
qc_pass <- qc_pass[!is.na(qc_pass)]
cat(sprintf("  QC 通过 (Coverage>=99%% & Dup<=0.7): %d 样本\n", length(qc_pass)))

# # Panel APOBEC (QC 过滤后 + total_snvs >= 5)
# panel_apobec <- sig[!is.na(sig$dnum) & sig$dnum %in% qc_pass, ]
# cat(sprintf("  Panel APOBEC (QC后): %d 样本\n", nrow(panel_apobec)))
# # 过滤 total_snvs < 5
# n_before <- nrow(panel_apobec)
# panel_apobec <- panel_apobec[panel_apobec$total_snvs >= 5, ]
# cat(sprintf("  过滤 total_snvs < 5: %d -> %d 样本 (去除 %d)\n", n_before, nrow(panel_apobec), n_before - nrow(panel_apobec)))
# panel_apobec <- panel_apobec[, c("dnum", "Signature_APOBEC_ml")]
# colnames(panel_apobec) <- c("dnum", "Panel_SigMA")

# ============================================================================
# 3. 读取 WES 各方法 APOBEC (SBS2 + SBS13), 仅保留去重后保留的样本
# ============================================================================
cat("\n===== 3. 读取 WES 各方法 APOBEC 信号 =====\n")

# --- deconstructSigs ---
# 格式: 行=样本 (D116, D167...), 列=SBS签名, 值为分数
read_decon <- function(f) {
  d <- read.csv(f, row.names = 1)
  data.frame(
    dnum = extract_dnum(rownames(d)),
    APOBEC = d$SBS2 + d$SBS13,
    stringsAsFactors = FALSE
  )
}
des_primary <- read_decon(file.path(base_dir, "WES/deconstructSigs/deconstructSigs_primary_weights.csv"))
des_strict  <- read_decon(file.path(base_dir, "WES/deconstructSigs/deconstructSigs_strict_weights.csv"))

# --- MutationalPatterns ---
# 格式: 行=签名(SBS1-SBS60), 列=样本
read_mutpat <- function(f) {
  d <- read.csv(f, row.names = 1)
  sbs2  <- as.numeric(d["SBS2", ])
  sbs13 <- as.numeric(d["SBS13", ])
  data.frame(
    dnum = extract_dnum(colnames(d)),
    APOBEC = sbs2 + sbs13,
    stringsAsFactors = FALSE
  )
}
mp_primary_reg <- read_mutpat(file.path(base_dir, "WES/MutationalPatterns/MutPatterns_primary_regular_weights.csv"))
mp_primary_str <- read_mutpat(file.path(base_dir, "WES/MutationalPatterns/MutPatterns_primary_strict_weights.csv"))

# --- SigProfiler ---
# 格式: 行=样本, 列=SBS签名, 值为突变数 (需归一化)
read_sigprof <- function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  colnames(d)[1] <- "sample"
  sig_cols <- grep("^SBS", colnames(d), value = TRUE)
  d$total  <- rowSums(d[, sig_cols])
  data.frame(
    dnum = extract_dnum(d$sample),
    APOBEC = (d$SBS2 + d$SBS13) / d$total,  # 归一化为分数
    stringsAsFactors = FALSE
  )
}
sp_primary <- read_sigprof(file.path(base_dir, "WES/SigProfiler/primary/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"))
sp_strict  <- read_sigprof(file.path(base_dir, "WES/SigProfiler/strict/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"))

# --- 对每个 WES 方法: 去重 + 过滤 ---
# 策略: 每个方法内可能有 D 和 MP 两种格式的重复
# 对重复的 D 编号, 保留 APOBEC 值最大的 (与 WES stats 中按 Total.SNV 选择一致)
filter_wes_dedup <- function(wes_df, method_label) {
  wes_df <- wes_df[!is.na(wes_df$dnum), ]
  # 仅保留在 WES stats 去重中存在的 D 编号
  wes_df <- wes_df[wes_df$dnum %in% wes_dedup$dnum, ]
  # 对重复的 D 编号 (同一方法内 D 和 MP 都在), 保留 APOBEC 值最大的
  if (any(duplicated(wes_df$dnum))) {
    wes_df <- do.call(rbind, lapply(split(wes_df, wes_df$dnum), function(g) {
      g[which.max(g$APOBEC), ]
    }))
    rownames(wes_df) <- NULL
  }
  cat(sprintf("  %s: %d 样本 (去重后)\n", method_label, nrow(wes_df)))
  return(wes_df)
}

cat("\n--- 去重过滤后 ---\n")
des_primary <- filter_wes_dedup(des_primary, "deconstructSigs primary")
des_strict  <- filter_wes_dedup(des_strict,  "deconstructSigs strict")
mp_primary_reg <- filter_wes_dedup(mp_primary_reg, "MutPatterns primary reg")
mp_primary_str <- filter_wes_dedup(mp_primary_str, "MutPatterns primary strict")
sp_primary  <- filter_wes_dedup(sp_primary, "SigProfiler primary")
sp_strict   <- filter_wes_dedup(sp_strict,  "SigProfiler strict")

# ============================================================================
# 4. 合并并计算一致性
# ============================================================================
cat("\n===== 4. 一致性分析 =====\n")

wes_methods <- list(
  "deconstructSigs\n(primary)"    = des_primary,
  "deconstructSigs\n(strict)"     = des_strict,
  "MutPatterns\n(primary,reg)"    = mp_primary_reg,
  "MutPatterns\n(primary,strict)" = mp_primary_str,
  "SigProfiler\n(primary)"        = sp_primary,
  "SigProfiler\n(strict)"         = sp_strict
)

results <- list()
for (method_name in names(wes_methods)) {
  wes <- wes_methods[[method_name]]
  colnames(wes)[2] <- "WES_APOBEC"
  
  merged <- merge(panel_apobec, wes, by = "dnum")
  cat(sprintf("  %s: %d 共有样本\n", gsub("\n", " ", method_name), nrow(merged)))
  
  if (nrow(merged) >= 5) {
    pearson  <- cor.test(merged$WES_APOBEC, merged$Panel_SigMA, method = "pearson")
    spearman <- cor.test(merged$WES_APOBEC, merged$Panel_SigMA, method = "spearman")
    cat(sprintf("    Pearson:  r=%.4f, p=%.2e\n", pearson$estimate, pearson$p.value))
    cat(sprintf("    Spearman: rho=%.4f, p=%.2e\n", spearman$estimate, spearman$p.value))
    results[[method_name]] <- list(
      data = merged,
      pearson = pearson,
      spearman = spearman
    )
  } else {
    cat("    样本太少，跳过\n")
  }
}

# ============================================================================
# 5. corrplot 可视化 (参考 corr_plot.R 风格)
# ============================================================================
cat("\n===== 5. 绘制 corrplot =====\n")

# 构建所有方法的全矩阵 (7×7)
method_labels <- c("Panel(SigMA)", "deS(primary)", "deS(strict)",
                   "MP(reg)", "MP(strict)", "SP(primary)", "SP(strict)")

# 收集所有方法的 APOBEC 值, 按 dnum 合并
all_apobec <- panel_apobec
colnames(all_apobec)[2] <- "Panel(SigMA)"

wes_label_map <- list(
  "deconstructSigs\n(primary)"    = "deS(primary)",
  "deconstructSigs\n(strict)"     = "deS(strict)",
  "MutPatterns\n(primary,reg)"    = "MP(reg)",
  "MutPatterns\n(primary,strict)" = "MP(strict)",
  "SigProfiler\n(primary)"        = "SP(primary)",
  "SigProfiler\n(strict)"         = "SP(strict)"
)

for (method_name in names(wes_methods)) {
  wes <- wes_methods[[method_name]]
  short <- wes_label_map[[method_name]]
  d <- wes[!is.na(wes$dnum), c("dnum", "APOBEC")]
  colnames(d)[2] <- short
  all_apobec <- merge(all_apobec, d, by = "dnum")
}

# 排序列确保一致
all_apobec <- all_apobec[, c("dnum", method_labels)]
cat(sprintf("  全矩阵样本数: %d\n", nrow(all_apobec)))

mat_data <- all_apobec[, method_labels]

# 计算相关矩阵和 p 值矩阵
cor_pearson <- cor(mat_data, method = "pearson", use = "pairwise.complete.obs")
cor_spearman <- cor(mat_data, method = "spearman", use = "pairwise.complete.obs")

pearson_pmat  <- cor.mtest(mat_data, method = "pearson")$p
spearman_pmat <- cor.mtest(mat_data, method = "spearman")$p

# --- Pearson corrplot ---
pdf(file.path(base_dir, "APOBEC_corrplot_Pearson.pdf"), width = 5, height = 5, family = "ArialMT")
corrplot(corr = cor_pearson,
         p.mat = pearson_pmat,
         sig.level = c(0.01, 0.05, 0.1),
         pch.cex = 0.5,
         pch.col = "white",
         insig = "pch",
         type = "lower",
         col = rev(COL2('RdBu', 200)),
         tl.col = "black",
         tl.cex = 0.7,
         tl.srt = 45,
         cl.ratio = 0.2,
         mar = c(0, 0, 0, 0),
         method = "number",
         number.cex = 0.5,
         number.digits = 2)
dev.off()

# --- Spearman corrplot ---
pdf(file.path(base_dir, "APOBEC_corrplot_Spearman.pdf"), width = 5, height = 5, family = "ArialMT")
corrplot(corr = cor_spearman,
         p.mat = spearman_pmat,
         sig.level = c(0.01, 0.05, 0.1),
         pch.cex = 0.5,
         pch.col = "white",
         insig = "pch",
         type = "lower",
         col = rev(COL2('RdBu', 200)),
         tl.col = "black",
         tl.cex = 0.7,
         tl.srt = 45,
         cl.ratio = 0.2,
         mar = c(0, 0, 0, 0),
         method = "number",
         number.cex = 0.5,
         number.digits = 2)
dev.off()

cat("  corrplot 已保存\n")

# ============================================================================
# 5b. 输出 719 共有样本 QC 结果及统计
# ============================================================================
cat("\n===== 5b. 输出共有样本 QC 结果 =====\n")

# 获取 719 共有样本的 dnum
common_dnums <- all_apobec$dnum
cat(sprintf("  共有样本数: %d\n", length(common_dnums)))

# --- WES QC ---
wes_qc <- wes_stats[wes_stats$dnum %in% common_dnums, c("Sample", "dnum", "Total.SNV", "PASS.SNV", "ALT.ge3.SNV", "DP.ge20.SNV", "ExAC_EAS.lt001.SNV", "Filtered.SNV", "ALT.ge5.SNV", "DP.ge30.SNV", "Filtered.strict.SNV")]
# 去重: 每个 dnum 保留 Total.SNV 最大的
wes_qc <- do.call(rbind, lapply(split(wes_qc, wes_qc$dnum), function(g) g[which.max(g$Total.SNV), ]))
rownames(wes_qc) <- NULL

# --- Panel QC ---
panel_qc <- qc[qc$dnum %in% common_dnums, c("dnum", "Sample", "Mapping", "cleanQ30", "TotalData", "cleanData", "Coverage", "Depth", "Duplication")]
# 去重: 每个 dnum 保留 Coverage 最大的
panel_qc <- do.call(rbind, lapply(split(panel_qc, panel_qc$dnum), function(g) g[which.max(as.numeric(g$Coverage)), ]))
rownames(panel_qc) <- NULL

# 合并 WES + Panel QC (先重命名 Sample 列避免冲突)
wes_qc_sub <- wes_qc[, c("dnum", "Sample", "Total.SNV", "PASS.SNV", "ALT.ge3.SNV", "DP.ge20.SNV", "ExAC_EAS.lt001.SNV", "Filtered.SNV", "ALT.ge5.SNV", "DP.ge30.SNV", "Filtered.strict.SNV")]
colnames(wes_qc_sub)[2] <- "Sample_WES"
panel_qc_sub <- panel_qc[, c("dnum", "Sample", "Mapping", "cleanQ30", "TotalData", "cleanData", "Coverage", "Depth", "Duplication")]
colnames(panel_qc_sub)[2] <- "Sample_Panel"

merged_qc <- merge(wes_qc_sub, panel_qc_sub, by = "dnum")

write.csv(merged_qc, file.path(base_dir, "matched_719_QC_detail.csv"), row.names = FALSE, quote = FALSE)
cat(sprintf("  QC 明细: matched_719_QC_detail.csv (%d 样本)\n", nrow(merged_qc)))

# --- QC 统计 (类似 QC_statistics.csv) ---
qc_metrics_wes <- c("Total.SNV", "PASS.SNV", "ALT.ge3.SNV", "DP.ge20.SNV", "ExAC_EAS.lt001.SNV", "Filtered.SNV", "ALT.ge5.SNV", "DP.ge30.SNV", "Filtered.strict.SNV")
qc_metrics_panel <- c("Mapping", "cleanQ30", "TotalData", "cleanData", "Coverage", "Depth", "Duplication")

stats_list <- list()
for (m in qc_metrics_wes) {
  vals <- as.numeric(merged_qc[[m]])
  qs <- as.numeric(quantile(vals, c(0.25, 0.75), na.rm = TRUE))
  stats_list[[m]] <- data.frame(Metric = m, Min = min(vals, na.rm=TRUE), Q1 = qs[1],
                                Median = median(vals, na.rm=TRUE), Mean = mean(vals, na.rm=TRUE),
                                Q3 = qs[2], Max = max(vals, na.rm=TRUE),
                                stringsAsFactors = FALSE)
}
for (m in qc_metrics_panel) {
  vals <- as.numeric(merged_qc[[m]])
  qs <- as.numeric(quantile(vals, c(0.25, 0.75), na.rm = TRUE))
  stats_list[[m]] <- data.frame(Metric = m, Min = min(vals, na.rm=TRUE), Q1 = qs[1],
                                Median = median(vals, na.rm=TRUE), Mean = mean(vals, na.rm=TRUE),
                                Q3 = qs[2], Max = max(vals, na.rm=TRUE),
                                stringsAsFactors = FALSE)
}

qc_stats_df <- do.call(rbind, stats_list)
rownames(qc_stats_df) <- NULL

write.csv(qc_stats_df, file.path(base_dir, "matched_719_QC_statistics.csv"), row.names = FALSE, quote = FALSE)
cat("  QC 统计: matched_719_QC_statistics.csv\n")
print(qc_stats_df)

# ============================================================================
# 6. 散点图
# ============================================================================
cat("\n===== 6. 绘制散点图 =====\n")

for (method_name in names(results)) {
  r <- results[[method_name]]
  d <- r$data
  method_clean <- gsub("[\n(), ]", "_", method_name)
  
  p <- ggplot(d, aes(x = WES_APOBEC, y = Panel_SigMA)) +
    geom_point(alpha = 0.4, size = 1.2, color = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8) +
    labs(
      title = method_name,
      subtitle = sprintf("Pearson r=%.3f (p=%.2e)  |  Spearman rho=%.3f (p=%.2e)  |  n=%d",
        r$pearson$estimate, r$pearson$p.value,
        r$spearman$estimate, r$spearman$p.value, nrow(d)),
      x = "WES APOBEC (SBS2 + SBS13)",
      y = "Panel APOBEC (SigMA)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  
  ggsave(file.path(base_dir, sprintf("APOBEC_scatter_%s.pdf", method_clean)),
         p, width = 7, height = 5)
}

# ============================================================================
# 7. 保存汇总表格
# ============================================================================
cat("\n===== 7. 保存汇总 =====\n")

summary_df <- data.frame(
  Method = gsub("\n", " ", names(wes_methods)),
  N_samples = sapply(names(wes_methods), function(m) {
    if (m %in% names(results)) nrow(results[[m]]$data) else NA
  }),
  Pearson_r = sapply(names(wes_methods), function(m) {
    if (m %in% names(results)) results[[m]]$pearson$estimate else NA
  }),
  Pearson_p = sapply(names(wes_methods), function(m) {
    if (m %in% names(results)) results[[m]]$pearson$p.value else NA
  }),
  Spearman_rho = sapply(names(wes_methods), function(m) {
    if (m %in% names(results)) results[[m]]$spearman$estimate else NA
  }),
  Spearman_p = sapply(names(wes_methods), function(m) {
    if (m %in% names(results)) results[[m]]$spearman$p.value else NA
  }),
  row.names = NULL
)

write.csv(summary_df, file.path(base_dir, "APOBEC_correlation_summary.csv"),
          row.names = FALSE, quote = FALSE)
print(summary_df)

cat("\n===== 完成 =====\n")
cat(sprintf("  汇总: %s/APOBEC_correlation_summary.csv\n", base_dir))
cat(sprintf("  Pearson:  %s/APOBEC_corrplot_Pearson.pdf\n", base_dir))
cat(sprintf("  Spearman: %s/APOBEC_corrplot_Spearman.pdf\n", base_dir))
cat(sprintf("  散点图: %s/APOBEC_scatter_*.pdf\n", base_dir))
