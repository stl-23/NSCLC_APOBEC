#!/usr/bin/env Rscript
# ============================================================================
# Step 2: 根据 Panel BED 文件提取 WES VCF 中的 Panel 区域突变
# BED 文件: chr开头, 需去掉 chr 与 VCF 的 CHROM 对应
# 判断 VCF 的 POS 是否落在 BED 区间 [start, end] 内
# ============================================================================

library(data.table)

# --- 路径配置 ---
base_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
bed_file <- "/Users/stl/Documents/AmoyDX/1学习/MasterPanel/master.snvindel.20251214.renameanno.manualupdate.bed"
vcf_dir  <- file.path(base_dir, "vcfs")
panel_dir <- file.path(base_dir, "panel_vcfs")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# --- 读取 BED 文件 ---
cat("===== 1. 读取 Panel BED =====\n")
bed <- fread(bed_file, sep = "\t", header = FALSE)
colnames(bed) <- c("chrom_raw", "start", "end", "gene_info")
# 去掉 chr 前缀
bed$chrom <- sub("^chr", "", bed$chrom_raw)
cat(sprintf("  BED 区间数: %d\n", nrow(bed)))
cat(sprintf("  染色体: %s\n", paste(sort(unique(bed$chrom)), collapse = ", ")))

# 按染色体分组 BED 区间
bed_list <- split(bed, bed$chrom)

# --- 读取过滤后 SNV 数据 ---
cat("\n===== 2. 读取过滤后 SNV =====\n")
snv_file <- file.path(base_dir, "TCGA_LUAD_filtered_SNV.tsv")
snvs <- fread(snv_file, sep = "\t", header = TRUE)
cat(sprintf("  总 SNV 数: %d\n", nrow(snvs)))
cat(sprintf("  总样本数: %d\n", length(unique(snvs$Tumor_Sample_Barcode))))

# --- 提取 Panel 区域突变 ---
cat("\n===== 3. 提取 Panel 区域突变 =====\n")

# 方法: 对每个 SNV, 检查其 CHROM+POS 是否落在任意 BED 区间内
# 先按染色体分组处理提高效率
snvs$in_panel <- FALSE

for (chr_name in unique(snvs$Chromosome)) {
  if (!(chr_name %in% names(bed_list))) next

  chr_snvs <- snvs[snvs$Chromosome == chr_name, ]
  chr_bed  <- bed_list[[chr_name]]

  # 对每个 SNV 检查是否落在任一 BED 区间
  for (i in seq_len(nrow(chr_snvs))) {
    pos <- chr_snvs$Start_Position[i]
    # 检查 POS 是否在任意 [start, end] 区间内
    hits <- chr_bed$start <= pos & chr_bed$end >= pos
    if (any(hits)) {
      idx <- which(snvs$Chromosome == chr_name)[i]
      snvs$in_panel[idx] <- TRUE
    }
  }
}

n_panel <- sum(snvs$in_panel)
cat(sprintf("  Panel 区域内 SNV: %d / %d (%.1f%%)\n",
            n_panel, nrow(snvs), 100 * n_panel / nrow(snvs)))

panel_snvs <- snvs[snvs$in_panel == TRUE, ]

# --- 按样本输出 Panel VCF ---
cat("\n===== 4. 生成 Panel VCF =====\n")
samples <- unique(panel_snvs$Tumor_Sample_Barcode)
snv_counts <- data.table(Sample = character(), Panel_SNV_Count = integer())

for (i in seq_along(samples)) {
  samp <- samples[i]
  sub <- panel_snvs[Tumor_Sample_Barcode == samp, ]

  vcf_file <- file.path(panel_dir, sprintf("%s.vcf", samp))

  con <- file(vcf_file, "w")
  writeLines("##fileformat=VCFv4.2", con)
  writeLines(sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", samp), con)

  for (j in seq_len(nrow(sub))) {
    info_str <- sprintf("Gene.refGene=%s;GeneDetail.refGene=%s;AAChange.refGene=%s:%s:%s",
                        sub$Hugo_Symbol[j], sub$Variant_Classification[j],
                        sub$Hugo_Symbol[j], sub$HGVSc[j], sub$HGVSp_Short[j])
    fmt_str <- sprintf("0/1:%d:%.4f:%d",
                       sub$t_alt_count[j], sub$t_vaf[j], sub$t_depth[j])
    line <- sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t%s\tGT:AD:AF:DP\t%s",
                    sub$Chromosome[j], sub$Start_Position[j],
                    sub$Reference_Allele[j], sub$Tumor_Seq_Allele2[j],
                    info_str, fmt_str)
    writeLines(line, con)
  }
  close(con)

  snv_counts <- rbind(snv_counts, data.table(Sample = samp, Panel_SNV_Count = nrow(sub)))

  if (i %% 50 == 0) cat(sprintf("  已处理 %d/%d 样本\n", i, length(samples)))
}

cat(sprintf("  共生成 %d 个 Panel VCF 文件\n", length(samples)))

# --- 统计 ---
cat("\n===== 5. Panel SNV 统计 =====\n")
snv_counts <- snv_counts[order(-Panel_SNV_Count), ]
cat(sprintf("  样本数: %d\n", nrow(snv_counts)))
cat(sprintf("  Panel SNV 总数: %d\n", sum(snv_counts$Panel_SNV_Count)))
cat(sprintf("  每样本: min=%d, Q1=%d, median=%d, mean=%.1f, Q3=%d, max=%d\n",
            min(snv_counts$Panel_SNV_Count),
            quantile(snv_counts$Panel_SNV_Count, 0.25),
            median(snv_counts$Panel_SNV_Count),
            mean(snv_counts$Panel_SNV_Count),
            quantile(snv_counts$Panel_SNV_Count, 0.75),
            max(snv_counts$Panel_SNV_Count)))

fwrite(snv_counts, file.path(base_dir, "TCGA_LUAD_Panel_SNV_statistics.csv"))
cat("  统计文件: TCGA_LUAD_Panel_SNV_statistics.csv\n")

# 保存 Panel 区域合并数据
fwrite(panel_snvs[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
                       Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
                       HGVSc, HGVSp_Short, t_alt_count, t_ref_count, t_depth, t_vaf)],
       file.path(base_dir, "TCGA_LUAD_panel_SNV.tsv"), sep = "\t")
cat(sprintf("  Panel 数据: TCGA_LUAD_panel_SNV.tsv (%d 位点)\n", nrow(panel_snvs)))

cat("\n===== 完成 =====\n")
