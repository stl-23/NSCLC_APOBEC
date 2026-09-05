#!/usr/bin/env Rscript
# ============================================================================
# TCGA-LUAD WES vs Panel 突变特征一致性分析
# 比较 3 种 WES 方法 (deconstructSigs / MutationalPatterns / SigProfiler)
# 与 Panel (SigMA) 的 APOBEC 信号一致性
#
# 参考: 3.WES_panel_APOBEC_correlation_ori.R
# ============================================================================

library(openxlsx)
library(corrplot)
library(ggplot2)

base_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
out_dir  <- base_dir

# ============================================================================
# 辅助函数
# ============================================================================

# 从 SigMA tumor 路径中提取样本名 (TCGA-XX-XXXX-XX)
extract_sample_from_path <- function(x) {
  # 从路径中提取 .vcf 前的样本名
  gsub(".*/(.+)\\.vcf", "\\1", x)
}

# ============================================================================
# 1. 读取 Panel SigMA APOBEC 值
# ============================================================================
cat("\n===== 1. 读取 Panel SigMA 结果 =====\n")
sig <- read.xlsx(file.path(base_dir, "TCGA_LUAD_sigMA_signature_lite.xlsx"))
sig$sample <- extract_sample_from_path(sig$tumor)
cat(sprintf("  SigMA 样本数: %d\n", nrow(sig)))
cat(sprintf("  total_snvs 范围: %d - %d\n", min(sig$total_snvs), max(sig$total_snvs)))

# Panel APOBEC
panel_apobec <- sig[, c("sample", "Signature_APOBEC_ml", "total_snvs")]
colnames(panel_apobec) <- c("sample", "Panel_SigMA", "total_snvs")

# ============================================================================
# 2. 读取 WES 各方法 APOBEC (SBS2 + SBS13)
# ============================================================================
cat("\n===== 2. 读取 WES 各方法 APOBEC 信号 =====\n")

# --- deconstructSigs ---
# 格式: 行=样本 (TCGA-05-4249-01), 列=SBS签名, 值为分数
read_decon <- function(f) {
  d <- read.csv(f, row.names = 1)
  data.frame(
    sample = rownames(d),
    APOBEC = d$SBS2 + d$SBS13,
    stringsAsFactors = FALSE
  )
}
des_primary <- read_decon(file.path(base_dir, "WES_signatures/deconstructSigs/deconstructSigs_primary_weights.csv"))
cat(sprintf("  deconstructSigs: %d 样本\n", nrow(des_primary)))

# --- MutationalPatterns ---
# 格式: 行=签名(SBS1-SBS60), 列=样本
read_mutpat <- function(f) {
  d <- read.csv(f, row.names = 1, check.names = FALSE)
  sbs2  <- as.numeric(d["SBS2", ])
  sbs13 <- as.numeric(d["SBS13", ])
  data.frame(
    sample = colnames(d),
    APOBEC = sbs2 + sbs13,
    stringsAsFactors = FALSE
  )
}
mp_primary <- read_mutpat(file.path(base_dir, "WES_signatures/MutationalPatterns/MutPatterns_primary_weights.csv"))
cat(sprintf("  MutationalPatterns: %d 样本\n", nrow(mp_primary)))

# --- SigProfiler ---
# 格式: 行=样本, 列=SBS签名, 值为突变数 (需归一化)
read_sigprof <- function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  colnames(d)[1] <- "sample"
  sig_cols <- grep("^SBS", colnames(d), value = TRUE)
  d$total  <- rowSums(d[, sig_cols])
  data.frame(
    sample = d$sample,
    APOBEC = (d$SBS2 + d$SBS13) / d$total,  # 归一化为分数
    stringsAsFactors = FALSE
  )
}
sp_primary <- read_sigprof(file.path(base_dir, "WES_signatures/SigProfiler/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"))
cat(sprintf("  SigProfiler: %d 样本\n", nrow(sp_primary)))

# ============================================================================
# 3. 合并并计算一致性
# ============================================================================
cat("\n===== 3. 一致性分析 =====\n")

wes_methods <- list(
  "deconstructSigs"    = des_primary,
  "MutPatterns"        = mp_primary,
  "SigProfiler"        = sp_primary
)

results <- list()
for (method_name in names(wes_methods)) {
  wes <- wes_methods[[method_name]]
  colnames(wes)[2] <- "WES_APOBEC"
  
  merged <- merge(panel_apobec, wes[, c("sample", "WES_APOBEC")], by = "sample")
  cat(sprintf("  %s: %d 共有样本\n", method_name, nrow(merged)))
  
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
# 4. corrplot 可视化
# ============================================================================
cat("\n===== 4. 绘制 corrplot =====\n")

# 构建所有方法的全矩阵
method_labels <- c("Panel(SigMA)", "deS", "MP", "SP")

# 收集所有方法的 APOBEC 值, 按 sample 合并
all_apobec <- panel_apobec[, c("sample", "Panel_SigMA")]
colnames(all_apobec)[2] <- "Panel(SigMA)"

wes_label_map <- list(
  "deconstructSigs" = "deS",
  "MutPatterns"     = "MP",
  "SigProfiler"     = "SP"
)

for (method_name in names(wes_methods)) {
  wes <- wes_methods[[method_name]]
  short <- wes_label_map[[method_name]]
  d <- wes[!is.na(wes$sample), c("sample", "APOBEC")]
  colnames(d)[2] <- short
  all_apobec <- merge(all_apobec, d, by = "sample")
}

# 排序列确保一致
all_apobec <- all_apobec[, c("sample", method_labels)]
cat(sprintf("  全矩阵样本数: %d\n", nrow(all_apobec)))

mat_data <- as.matrix(all_apobec[, method_labels])

# 计算相关矩阵和 p 值矩阵
cor_pearson  <- cor(mat_data, method = "pearson",  use = "pairwise.complete.obs")
cor_spearman <- cor(mat_data, method = "spearman", use = "pairwise.complete.obs")

pearson_pmat  <- cor.mtest(mat_data, method = "pearson")$p
spearman_pmat <- cor.mtest(mat_data, method = "spearman")$p

# --- Pearson corrplot ---
pdf(file.path(out_dir, "TCGA_LUAD_APOBEC_corrplot_Pearson.pdf"), width = 5, height = 5, family = "ArialMT")
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
         number.cex = 0.7,
         number.digits = 2)
dev.off()

# --- Spearman corrplot ---
pdf(file.path(out_dir, "TCGA_LUAD_APOBEC_corrplot_Spearman.pdf"), width = 5, height = 5, family = "ArialMT")
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
         number.cex = 0.7,
         number.digits = 2)
dev.off()

cat("  corrplot 已保存\n")

# ============================================================================
# 5. 散点图
# ============================================================================
cat("\n===== 5. 绘制散点图 =====\n")

for (method_name in names(results)) {
  r <- results[[method_name]]
  d <- r$data
  
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
  
  ggsave(file.path(out_dir, sprintf("TCGA_LUAD_APOBEC_scatter_%s.pdf", method_name)),
         p, width = 7, height = 5)
}

cat("  散点图已保存\n")

# ============================================================================
# 6. 保存汇总表格
# ============================================================================
cat("\n===== 6. 保存汇总 =====\n")

summary_df <- data.frame(
  Method = names(wes_methods),
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

write.csv(summary_df, file.path(out_dir, "TCGA_LUAD_APOBEC_correlation_summary.csv"),
          row.names = FALSE, quote = FALSE)
cat("  汇总表:\n")
print(summary_df)

cat("\n===== 完成 =====\n")
cat(sprintf("  汇总: %s/TCGA_LUAD_APOBEC_correlation_summary.csv\n", out_dir))
cat(sprintf("  Pearson:  %s/TCGA_LUAD_APOBEC_corrplot_Pearson.pdf\n", out_dir))
cat(sprintf("  Spearman: %s/TCGA_LUAD_APOBEC_corrplot_Spearman.pdf\n", out_dir))
cat(sprintf("  散点图: %s/TCGA_LUAD_APOBEC_scatter_*.pdf\n", out_dir))
