#!/usr/bin/env Rscript
# ============================================================================
# Step 1: TCGA-LUAD MAF -> 按样本生成过滤后的 VCF 文件
# 过滤条件: t_alt_count >= 3 & (t_alt_count + t_ref_count) >= 20 &
#           t_alt_count / (t_alt_count + t_ref_count) >= 0.05
# 仅保留 SNV (Variant_Type == "SNP", Start == End, nchar(REF)==1 & nchar(ALT)==1)
# ============================================================================

library(data.table)

# --- 路径配置 ---
maf_file  <- "/Volumes/Elements/Amoydx/public_datasets/luad_tcga/data_mutations_extended.txt"
out_dir   <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
vcf_dir   <- file.path(out_dir, "vcfs")
dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

# --- 读取 MAF ---
cat("===== 1. 读取 MAF =====\n")
maf <- fread(maf_file, sep = "\t", header = TRUE, quote = "")
cat(sprintf("  总突变数: %d\n", nrow(maf)))
cat(sprintf("  总样本数: %d\n", length(unique(maf$Tumor_Sample_Barcode))))

# --- 过滤 ---
cat("\n===== 2. 过滤 =====\n")
maf$t_depth   <- maf$t_alt_count + maf$t_ref_count
maf$t_vaf     <- maf$t_alt_count / maf$t_depth

n_total <- nrow(maf)

# 仅保留 SNP
maf_snp <- maf[maf$Variant_Type == "SNP", ]
cat(sprintf("  SNP: %d (去除 %d non-SNP)\n", nrow(maf_snp), n_total - nrow(maf_snp)))

# t_alt_count >= 3
maf_filt <- maf_snp[maf_snp$t_alt_count >= 3, ]
cat(sprintf("  t_alt_count >= 3: %d\n", nrow(maf_filt)))

# t_depth >= 20
maf_filt <- maf_filt[maf_filt$t_depth >= 20, ]
cat(sprintf("  t_depth >= 20: %d\n", nrow(maf_filt)))

# VAF >= 5%
maf_filt <- maf_filt[maf_filt$t_vaf >= 0.05, ]
cat(sprintf("  VAF >= 5%%: %d\n", nrow(maf_filt)))

# 仅保留单碱基替换 (Start == End, nchar(REF)==1, nchar(ALT)==1)
# 注意: Tumor_Seq_Allele2 是真正的变异等位基因 (Tumor_Seq_Allele1 = Reference)
maf_filt <- maf_filt[maf_filt$Start_Position == maf_filt$End_Position, ]
maf_filt <- maf_filt[nchar(maf_filt$Reference_Allele) == 1 & nchar(maf_filt$Tumor_Seq_Allele2) == 1, ]
cat(sprintf("  单碱基 SNV: %d\n", nrow(maf_filt)))

# --- 构建 VCF 字段 ---
cat("\n===== 3. 生成 VCF =====\n")

# INFO 字段
maf_filt[, INFO := sprintf("Gene.refGene=%s;GeneDetail.refGene=%s;AAChange.refGene=%s:%s:%s",
                           Hugo_Symbol,
                           Variant_Classification,
                           Hugo_Symbol, HGVSc, HGVSp_Short)]

# FORMAT 字段: GT:AD:AF:DP
# 样本列: 0/1:{t_alt_count}:{VAF}:{t_depth}
maf_filt[, sample_fmt := sprintf("0/1:%d:%.4f:%d",
                                  t_alt_count, t_vaf, t_depth)]

# --- 按样本输出 VCF ---
samples <- unique(maf_filt$Tumor_Sample_Barcode)
snv_counts <- data.table(Sample = character(), SNV_Count = integer())

for (i in seq_along(samples)) {
  samp <- samples[i]
  sub <- maf_filt[Tumor_Sample_Barcode == samp, ]

  vcf_file <- file.path(vcf_dir, sprintf("%s.vcf", samp))

  # 写入 VCF
  con <- file(vcf_file, "w")
  writeLines("##fileformat=VCFv4.2", con)
  writeLines(sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", samp), con)

  for (j in seq_len(nrow(sub))) {
    line <- sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t%s\tGT:AD:AF:DP\t%s",
                    sub$Chromosome[j],
                    sub$Start_Position[j],
                    sub$Reference_Allele[j],
                    sub$Tumor_Seq_Allele2[j],
                    sub$INFO[j],
                    sub$sample_fmt[j])
    writeLines(line, con)
  }
  close(con)

  snv_counts <- rbind(snv_counts, data.table(Sample = samp, SNV_Count = nrow(sub)))

  if (i %% 50 == 0) cat(sprintf("  已处理 %d/%d 样本\n", i, length(samples)))
}

cat(sprintf("  共生成 %d 个 VCF 文件\n", length(samples)))

# --- SNV 统计 ---
cat("\n===== 4. SNV 统计 =====\n")
snv_counts <- snv_counts[order(-SNV_Count), ]
cat(sprintf("  样本数: %d\n", nrow(snv_counts)))
cat(sprintf("  SNV 总数: %d\n", sum(snv_counts$SNV_Count)))
cat(sprintf("  每样本 SNV: min=%d, Q1=%d, median=%d, mean=%.1f, Q3=%d, max=%d\n",
            min(snv_counts$SNV_Count),
            quantile(snv_counts$SNV_Count, 0.25),
            median(snv_counts$SNV_Count),
            mean(snv_counts$SNV_Count),
            quantile(snv_counts$SNV_Count, 0.75),
            max(snv_counts$SNV_Count)))

# 保存统计
fwrite(snv_counts, file.path(out_dir, "TCGA_LUAD_SNV_statistics.csv"))
cat(sprintf("  统计文件: TCGA_LUAD_SNV_statistics.csv\n"))

# 保存过滤后合并数据 (供 Step 2 使用)
fwrite(maf_filt[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
                     Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
                     HGVSc, HGVSp_Short, t_alt_count, t_ref_count, t_depth, t_vaf, INFO)],
       file.path(out_dir, "TCGA_LUAD_filtered_SNV.tsv"), sep = "\t")
cat(sprintf("  过滤后数据: TCGA_LUAD_filtered_SNV.tsv (%d 位点)\n", nrow(maf_filt)))

cat("\n===== 完成 =====\n")
