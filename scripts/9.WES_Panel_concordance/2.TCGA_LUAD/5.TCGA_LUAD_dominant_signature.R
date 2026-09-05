#!/usr/bin/env Rscript
# ============================================================================
# TCGA-LUAD WES vs Panel Dominant Signature 一致性分析
# 比较 WES 3 种方法 (deconstructSigs / MutationalPatterns / SigProfiler)
# 与 Panel (SigMA) 的 dominant signature 分类一致性
#
# Dominant signature 定义:
#   APOBEC3 = SBS2 + SBS13
#   HRD     = SBS3 + SBS8
#   Clock   = SBS1 + SBS5
#   Other   = 其余所有签名
#   dominant = which.max(c(APOBEC3, HRD, Clock, Other))
#
# APOBEC dominant 分组: dominant == APOBEC → 1, 否则 → 0
# 以 WES 为金标准, Panel 为预测, 评估 Sensitivity/Specificity/Accuracy
# ============================================================================

library(openxlsx)
library(ggplot2)

base_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
out_dir  <- base_dir

# ============================================================================
# 辅助函数
# ============================================================================

# 从 SigMA tumor 路径中提取样本名
extract_sample_from_path <- function(x) {
  gsub(".*/(.+)\\.vcf", "\\1", x)
}

# 计算 dominant signature 分类
compute_dominant_wes <- function(d, sbs_prefix = "SBS") {
  get_sbs <- function(name) {
    col <- paste0(sbs_prefix, name)
    if (col %in% colnames(d)) as.numeric(d[[col]]) else rep(0, nrow(d))
  }
  
  apobec <- get_sbs("2") + get_sbs("13")
  hrd    <- get_sbs("3") + get_sbs("8")
  clock  <- get_sbs("1") + get_sbs("5")
  total  <- rowSums(d[, grep(paste0("^", sbs_prefix), colnames(d)), drop = FALSE])
  other  <- pmax(total - apobec - hrd - clock, 0)
  
  cats <- cbind(APOBEC = apobec, HRD = hrd, Clock = clock, Other = other)
  dominant_idx <- apply(cats, 1, which.max)
  dominant <- colnames(cats)[dominant_idx]
  
  data.frame(
    sample = d$sample,
    APOBEC = apobec, HRD = hrd, Clock = clock, Other = other,
    dominant = dominant,
    is_apobec = as.integer(dominant == "APOBEC"),
    stringsAsFactors = FALSE
  )
}

# 性能评估
eval_performance <- function(wes_binary, panel_binary, label) {
  tp <- sum(wes_binary == 1 & panel_binary == 1)
  fp <- sum(wes_binary == 0 & panel_binary == 1)
  fn <- sum(wes_binary == 1 & panel_binary == 0)
  tn <- sum(wes_binary == 0 & panel_binary == 0)
  
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA)
  accuracy    <- ifelse(tp + tn + fp + fn > 0, (tp + tn) / (tp + tn + fp + fn), NA)
  precision   <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
  
  data.frame(
    Method = label,
    N = length(wes_binary),
    TP = tp, FP = fp, FN = fn, TN = tn,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    Precision = precision,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# 1. 读取 Panel SigMA 结果
# ============================================================================
cat("\n===== 1. 读取 Panel SigMA 结果 =====\n")
panel_sig <- read.xlsx(file.path(base_dir, "TCGA_LUAD_sigMA_signature_lite.xlsx"))
panel_sig$sample <- extract_sample_from_path(panel_sig$tumor)
cat(sprintf("  Panel SigMA: %d 样本\n", nrow(panel_sig)))

# Panel dominant
panel_cats <- data.frame(
  APOBEC = panel_sig$Signature_APOBEC_ml,
  HRD    = panel_sig$Signature_3_ml,
  Clock  = panel_sig$Signature_clock_ml,
  Other  = panel_sig$Signature_28_ml + panel_sig$Signature_18_ml +
           panel_sig$Signature_4_ml + panel_sig$Signature_16_ml
)
panel_dominant_idx <- apply(panel_cats, 1, which.max)
panel_sig$dominant  <- c("APOBEC", "HRD", "Clock", "Other")[panel_dominant_idx]
panel_sig$is_apobec <- as.integer(panel_sig$dominant == "APOBEC")

panel_apobec_df <- panel_sig[, c("sample", "dominant", "is_apobec")]
cat(sprintf("  Panel APOBEC dominant: %d / %d (%.1f%%)\n",
    sum(panel_sig$is_apobec), nrow(panel_sig),
    100 * mean(panel_sig$is_apobec)))
cat("  Panel dominant 分布:\n")
print(table(panel_sig$dominant))

# ============================================================================
# 2. 读取 WES 各方法结果并计算 dominant
# ============================================================================
cat("\n===== 2. 读取 WES 各方法并计算 dominant =====\n")

wes_methods <- list()

# --- deconstructSigs ---
# 格式: 行=样本(TCGA-05-4249-01), 列=SBS1-SBS60, 值=权重分数
{
  f <- file.path(base_dir, "WES_signatures/deconstructSigs/deconstructSigs_primary_weights.csv")
  d <- read.csv(f, row.names = 1)
  d$sample <- rownames(d)
  wes_methods[["deconstructSigs"]] <- compute_dominant_wes(d)
  cat(sprintf("  deconstructSigs: %d 样本\n", nrow(d)))
}

# --- MutationalPatterns ---
# 格式: 行=签名(SBS1-SBS60), 列=样本, 值=突变数
{
  f <- file.path(base_dir, "WES_signatures/MutationalPatterns/MutPatterns_primary_weights.csv")
  d_raw <- read.csv(f, row.names = 1, check.names = FALSE)
  # 转置: 行=样本, 列=签名
  d <- as.data.frame(t(d_raw))
  d[] <- lapply(d, as.numeric)
  d$sample <- rownames(d)
  wes_methods[["MutPatterns"]] <- compute_dominant_wes(d)
  cat(sprintf("  MutPatterns: %d 样本\n", nrow(d)))
}

# --- SigProfiler ---
# 格式: 行=样本, 列=SBS签名, 值=突变数
{
  f <- file.path(base_dir, "WES_signatures/SigProfiler/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$sample <- d$Samples
  wes_methods[["SigProfiler"]] <- compute_dominant_wes(d)
  cat(sprintf("  SigProfiler: %d 样本\n", nrow(d)))
}

# ============================================================================
# 3. 合并并评估性能
# ============================================================================
cat("\n===== 3. 一致性评估 =====\n")

perf_list <- list()
detail_list <- list()

for (method_name in names(wes_methods)) {
  wes_dom <- wes_methods[[method_name]]
  
  # 按 sample 合并
  merged <- merge(panel_apobec_df, wes_dom[, c("sample", "dominant", "is_apobec")],
                  by = "sample", suffixes = c("_panel", "_wes"))
  
  cat(sprintf("\n  %s: %d 共有样本\n", method_name, nrow(merged)))
  cat(sprintf("    WES APOBEC dominant: %d (%.1f%%)\n",
    sum(merged$is_apobec_wes), 100 * mean(merged$is_apobec_wes)))
  cat(sprintf("    Panel APOBEC dominant: %d (%.1f%%)\n",
    sum(merged$is_apobec_panel), 100 * mean(merged$is_apobec_panel)))
  
  # 性能评估
  perf <- eval_performance(merged$is_apobec_wes, merged$is_apobec_panel, method_name)
  perf_list[[method_name]] <- perf
  cat(sprintf("    Sensitivity=%.3f, Specificity=%.3f, Accuracy=%.3f\n",
    perf$Sensitivity, perf$Specificity, perf$Accuracy))
  
  detail_list[[method_name]] <- merged
}

# ============================================================================
# 4. 保存汇总
# ============================================================================
cat("\n===== 4. 保存结果 =====\n")

# 性能汇总表
perf_df <- do.call(rbind, perf_list)
rownames(perf_df) <- NULL
write.csv(perf_df, file.path(out_dir, "TCGA_LUAD_dominant_APOBEC_performance.csv"),
          row.names = FALSE, quote = FALSE)
cat("\n  性能汇总:\n")
print(perf_df)

# 各方法 dominant 分布对比
cat("\n  Dominant 分布对比:\n")
for (method_name in names(detail_list)) {
  d <- detail_list[[method_name]]
  cat(sprintf("\n  %s:\n", method_name))
  cat(sprintf("    WES:   APOBEC=%d  HRD=%d  Clock=%d  Other=%d\n",
    sum(d$dominant_wes=="APOBEC"), sum(d$dominant_wes=="HRD"),
    sum(d$dominant_wes=="Clock"), sum(d$dominant_wes=="Other")))
  cat(sprintf("    Panel: APOBEC=%d  HRD=%d  Clock=%d  Other=%d\n",
    sum(d$dominant_panel=="APOBEC"), sum(d$dominant_panel=="HRD"),
    sum(d$dominant_panel=="Clock"), sum(d$dominant_panel=="Other")))
}

# 保存详细数据
for (method_name in names(detail_list)) {
  fname <- sprintf("TCGA_LUAD_dominant_detail_%s.csv", method_name)
  write.csv(detail_list[[method_name]], file.path(out_dir, fname),
            row.names = FALSE, quote = FALSE)
}

# ============================================================================
# 5. 可视化
# ============================================================================
cat("\n===== 5. 可视化 =====\n")

# 5a. Dominant signature 分布堆叠图
all_details <- do.call(rbind, lapply(names(detail_list), function(m) {
  d <- detail_list[[m]]
  data.frame(Method = m, Source = "WES", dominant = d$dominant_wes, stringsAsFactors = FALSE)
}))
all_details_panel <- do.call(rbind, lapply(names(detail_list), function(m) {
  d <- detail_list[[m]]
  data.frame(Method = m, Source = "Panel", dominant = d$dominant_panel, stringsAsFactors = FALSE)
}))
all_dom <- rbind(all_details, all_details_panel)
all_dom$dominant <- factor(all_dom$dominant, levels = c("APOBEC", "HRD", "Clock", "Other"))

p1 <- ggplot(all_dom, aes(x = Method, fill = dominant)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Source, ncol = 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("APOBEC" = "#E41A1C", "HRD" = "#377EB8",
                                "Clock" = "#4DAF4A", "Other" = "#999999")) +
  labs(title = "TCGA-LUAD: Dominant Signature Distribution",
       x = "Method", y = "Proportion", fill = "Dominant") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, "TCGA_LUAD_dominant_signature_distribution.pdf"),
       p1, width = 8, height = 6)

# 5b. 性能指标柱状图
perf_long <- data.frame(
  Method = rep(perf_df$Method, 3),
  Metric = rep(c("Sensitivity", "Specificity", "Accuracy"), each = nrow(perf_df)),
  Value = c(perf_df$Sensitivity, perf_df$Specificity, perf_df$Accuracy)
)
perf_long$Metric <- factor(perf_long$Metric, levels = c("Sensitivity", "Specificity", "Accuracy"))

p2 <- ggplot(perf_long, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Sensitivity" = "#E41A1C", "Specificity" = "#377EB8",
                                "Accuracy" = "#4DAF4A")) +
  labs(title = "TCGA-LUAD: APOBEC Dominant Classification Performance",
       x = "WES Method", y = "Score", fill = "Metric") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylim(0, 1.15)

ggsave(file.path(out_dir, "TCGA_LUAD_dominant_APOBEC_performance.pdf"),
       p2, width = 8, height = 6)

cat("  图表已保存\n")

cat("\n===== 完成 =====\n")
cat(sprintf("  性能汇总: %s/TCGA_LUAD_dominant_APOBEC_performance.csv\n", out_dir))
cat(sprintf("  分布图: %s/TCGA_LUAD_dominant_signature_distribution.pdf\n", out_dir))
cat(sprintf("  性能图: %s/TCGA_LUAD_dominant_APOBEC_performance.pdf\n", out_dir))
