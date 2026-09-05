#!/usr/bin/env Rscript
# ============================================================================
# TCGA-LUAD Panel 突变特征分析 - SigMA
# 基于 panel_vcfs/ 中的 VCF 文件
# ============================================================================

library(openxlsx)
devtools::load_all("/Users/stl/Documents/pipelines/sigMA/SigMA/")

# --- 路径设置 ---
base_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
vcf_dir  <- file.path(base_dir, "panel_vcfs")

cat("===== TCGA-LUAD SigMA 突变特征分析 =====\n")
cat(sprintf("  VCF 目录: %s\n", vcf_dir))
cat(sprintf("  VCF 文件数: %d\n", length(list.files(vcf_dir, pattern = "[.]vcf$"))))

# ============================================================================
# 1. 构建 96 通道突变矩阵
# ============================================================================
cat("\n===== 1. 构建 96 通道突变矩阵 =====\n")
genomes_matrix <- make_matrix(vcf_dir, file_type = 'vcf', ref_genome_name = 'hg19')
genomes <- conv_snv_matrix_to_df(genomes_matrix)
cat(sprintf("  矩阵维度: %d x %d\n", nrow(genomes), ncol(genomes)))

# 保存 96 通道矩阵
genome_file <- file.path(base_dir, 'genome_matrix_96_panel.csv')
write.table(genomes, genome_file, sep = ',', row.names = FALSE,
            col.names = TRUE, quote = FALSE)
cat(sprintf("  96 通道矩阵已保存: genome_matrix_96_panel.csv\n"))

# ============================================================================
# 2. 运行 SigMA (lung 类型, MSK platform)
# ============================================================================
cat("\n===== 2. 运行 SigMA =====\n")
cat("  tumor_type=lung, data=msk, snv_cutoff=5\n")
output_file <- run(genome_file,
                   data = "msk",
                   do_assign = TRUE,
                   do_mva = FALSE,
                   tumor_type = 'lung',
                   #snv_cutoff = 5,
                   lite_format = FALSE)

cat(sprintf("  SigMA 输出: %s\n", output_file))

# ============================================================================
# 3. 保存结果
# ============================================================================
cat("\n===== 3. 保存结果 =====\n")
m <- read.csv(output_file)
lite <- lite_df(m)

write.xlsx(m,    file = file.path(base_dir, "TCGA_LUAD_sigMA_signature_all.xlsx"))
write.xlsx(lite, file = file.path(base_dir, "TCGA_LUAD_sigMA_signature_lite.xlsx"))

cat(sprintf("  通过样本: %d / %d\n", sum(m$pass_ml, na.rm = TRUE), nrow(m)))
cat(sprintf("  完整结果: TCGA_LUAD_sigMA_signature_all.xlsx\n"))
cat(sprintf("  精简结果: TCGA_LUAD_sigMA_signature_lite.xlsx\n"))
cat("\n===== 完成 =====\n")
