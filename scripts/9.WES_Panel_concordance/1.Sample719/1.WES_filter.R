#!/usr/bin/env Rscript
# ============================================================================
# WES 过滤脚本 - 两种严格程度
# 
# 宽松过滤 (primary):  PASS & ALT>=3 & DP>=20 & ExAC_EAS<0.01
# 严格过滤 (strict):   PASS & ALT>=5 & DP>=30 & ExAC_EAS<0.01
#
# VCF 格式说明:
#   - 倒数第二列为肿瘤样本: GT:AD:AF:DP:F1R2:F2R1:FAD:SB
#   - ALT reads = AF * DP
#   - population AF: INFO 列中 ExAC_EAS=xxx (分号前或行尾的数字)
#
# 输出:
#   - WES.filter.final.primary/  宽松过滤后的 VCF
#   - WES.filter.final.strict/   严格过滤后的 VCF
#   - WES_filter_statistics.csv  各步骤保留 SNV 数统计表
# ============================================================================

library(openxlsx)

# ============================================================================
# 1. 路径设置
# ============================================================================
wes_dir       <- "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter"
primary_dir   <- "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.primary"
strict_dir    <- "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.strict"

dir.create(primary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(strict_dir,  showWarnings = FALSE, recursive = TRUE)

vcf_files <- list.files(wes_dir, pattern = "\\.vcf$", full.names = TRUE)
cat(sprintf("发现 %d 个 VCF 文件\n", length(vcf_files)))

# ============================================================================
# 2. 初始化统计表
# ============================================================================
stats_list <- data.frame(
  Sample          = character(),
  Total.SNV       = integer(),
  PASS.SNV        = integer(),
  ALT.ge3.SNV     = integer(),
  DP.ge20.SNV     = integer(),
  ExAC_EAS.lt001.SNV = integer(),
  Filtered.SNV    = integer(),
  ALT.ge5.SNV     = integer(),
  DP.ge30.SNV     = integer(),
  Filtered.strict.SNV = integer(),
  stringsAsFactors = FALSE
)

# ============================================================================
# 3. 处理单个 VCF 文件
# ============================================================================
process_vcf <- function(vcf_path) {
  # ---------- 读取文件 ----------
  all_lines <- readLines(vcf_path, warn = FALSE)
  header_idx <- grep("^#CHROM", all_lines)
  if (length(header_idx) == 0) {
    warning(sprintf("跳过 %s: 未找到 #CHROM 行", basename(vcf_path)))
    return(NULL)
  }
  header_idx <- header_idx[1]
  header_line <- all_lines[header_idx]
  data_lines  <- all_lines[(header_idx + 1):length(all_lines)]
  data_lines  <- data_lines[nchar(trimws(data_lines)) > 0]

  if (length(data_lines) == 0) {
    # 无变异位点
    cols <- strsplit(header_line, "\t")[[1]]
    tumor_idx <- length(cols) - 1
    sample_name <- cols[tumor_idx]
    header_lines <- all_lines[1:header_idx]
    return(list(
      sample = sample_name,
      total_snv = 0L, pass_snv = 0L,
      alt3_snv = 0L, dp20_snv = 0L, exac_snv = 0L, filtered_snv = 0L,
      alt5_snv = 0L, dp30_snv = 0L, filtered_strict_snv = 0L,
      primary_lines = character(0), strict_lines = character(0),
      header_lines = header_lines
    ))
  }

  # ---------- 解析表头 ----------
  cols      <- strsplit(header_line, "\t")[[1]]
  tumor_idx <- length(cols) - 1   # 倒数第二列 = 肿瘤样本
  sample_name <- cols[tumor_idx]

  # ---------- 逐行解析 ----------
  n <- length(data_lines)
  is_snv       <- logical(n)
  is_pass      <- logical(n)
  alt_reads    <- numeric(n)
  dp_vals      <- numeric(n)
  exac_eas_val <- numeric(n)   # NA = 未注释

  for (i in seq_len(n)) {
    fields <- strsplit(data_lines[i], "\t")[[1]]

    # --- SNV 判断: REF 和 ALT 均为单碱基 ---
    ref <- fields[4]
    alt <- fields[5]
    is_snv[i] <- (nchar(ref) == 1) && (nchar(alt) == 1)

    # --- FILTER ---
    is_pass[i] <- (fields[7] == "PASS")

    # --- 肿瘤样本 FORMAT 字段 ---
    tumor_field <- fields[tumor_idx]
    format_parts <- strsplit(fields[9], ":")[[1]]
    tumor_parts  <- strsplit(tumor_field, ":")[[1]]

    # AF: FORMAT 中第 3 个字段
    af_idx <- which(format_parts == "AF")
    dp_idx <- which(format_parts == "DP")

    af <- NA_real_
    dp <- NA_real_
    if (length(af_idx) > 0 && length(tumor_parts) >= af_idx[1]) {
      af <- suppressWarnings(as.numeric(tumor_parts[af_idx[1]]))
    }
    if (length(dp_idx) > 0 && length(tumor_parts) >= dp_idx[1]) {
      dp <- suppressWarnings(as.numeric(tumor_parts[dp_idx[1]]))
    }

    dp_vals[i]   <- dp
    alt_reads[i] <- if (!is.na(af) && !is.na(dp)) af * dp else NA_real_

    # --- ExAC_EAS: INFO 列中提取 ---
    info <- fields[8]
    m <- regmatches(info, regexpr("ExAC_EAS=([0-9.]+)", info))
    if (length(m) > 0 && nchar(m) > 0) {
      exac_eas_val[i] <- as.numeric(sub("ExAC_EAS=", "", m))
    } else {
      exac_eas_val[i] <- NA_real_
    }
  }

  # ---------- 统计 (仅 SNV) ----------
  snv_idx <- which(is_snv)
  total_snv <- length(snv_idx)

  if (total_snv == 0) {
    header_lines <- all_lines[1:header_idx]
    return(list(
      sample = sample_name,
      total_snv = 0L, pass_snv = 0L,
      alt3_snv = 0L, dp20_snv = 0L, exac_snv = 0L, filtered_snv = 0L,
      alt5_snv = 0L, dp30_snv = 0L, filtered_strict_snv = 0L,
      primary_lines = character(0), strict_lines = character(0),
      header_lines = header_lines
    ))
  }

  # 各步过滤 (累积)
  pass_idx   <- snv_idx[is_pass[snv_idx]]
  pass_snv   <- length(pass_idx)

  # 宽松: ALT >= 3
  alt3_idx   <- pass_idx[!is.na(alt_reads[pass_idx]) & alt_reads[pass_idx] >= 3]
  alt3_snv   <- length(alt3_idx)

  # 宽松: DP >= 20
  dp20_idx   <- alt3_idx[!is.na(dp_vals[alt3_idx]) & dp_vals[alt3_idx] >= 20]
  dp20_snv   <- length(dp20_idx)

  # 宽松: ExAC_EAS < 0.01 (NA 视为通过)
  exac_idx   <- dp20_idx[is.na(exac_eas_val[dp20_idx]) | exac_eas_val[dp20_idx] < 0.01]
  exac_snv   <- length(exac_idx)
  filtered_snv <- exac_snv

  # 严格: ALT >= 5
  alt5_idx   <- pass_idx[!is.na(alt_reads[pass_idx]) & alt_reads[pass_idx] >= 5]
  alt5_snv   <- length(alt5_idx)

  # 严格: DP >= 30
  dp30_idx   <- alt5_idx[!is.na(dp_vals[alt5_idx]) & dp_vals[alt5_idx] >= 30]
  dp30_snv   <- length(dp30_idx)

  # 严格: ExAC_EAS < 0.01
  exac_strict_idx <- dp30_idx[is.na(exac_eas_val[dp30_idx]) | exac_eas_val[dp30_idx] < 0.01]
  exac_strict_snv <- length(exac_strict_idx)
  filtered_strict_snv <- exac_strict_snv

  # ---------- 收集输出行 ----------
  primary_idx <- exac_idx
  strict_idx  <- exac_strict_idx
  primary_out <- data_lines[primary_idx]
  strict_out  <- data_lines[strict_idx]

  list(
    sample = sample_name,
    total_snv = total_snv, pass_snv = pass_snv,
    alt3_snv = alt3_snv, dp20_snv = dp20_snv, exac_snv = exac_snv,
    filtered_snv = filtered_snv,
    alt5_snv = alt5_snv, dp30_snv = dp30_snv,
    filtered_strict_snv = filtered_strict_snv,
    primary_lines = primary_out, strict_lines = strict_out,
    header_lines = all_lines[1:header_idx]
  )
}

# ============================================================================
# 4. 批量处理所有 VCF
# ============================================================================
all_stats <- stats_list
n_files   <- length(vcf_files)

for (i in seq_along(vcf_files)) {
  vcf_path <- vcf_files[i]
  cat(sprintf("[%d/%d] 处理: %s\n", i, n_files, basename(vcf_path)))

  result <- process_vcf(vcf_path)
  if (is.null(result)) next

  # 写入统计表
  row <- data.frame(
    Sample             = result$sample,
    Total.SNV          = result$total_snv,
    PASS.SNV           = result$pass_snv,
    ALT.ge3.SNV        = result$alt3_snv,
    DP.ge20.SNV        = result$dp20_snv,
    ExAC_EAS.lt001.SNV = result$exac_snv,
    Filtered.SNV       = result$filtered_snv,
    ALT.ge5.SNV        = result$alt5_snv,
    DP.ge30.SNV        = result$dp30_snv,
    Filtered.strict.SNV = result$filtered_strict_snv,
    stringsAsFactors   = FALSE
  )
  all_stats <- rbind(all_stats, row)

  # 写入过滤后 VCF (宽松)
  out_primary <- file.path(primary_dir, basename(vcf_path))
  writeLines(c(result$header_lines, result$primary_lines), out_primary)
  
  # 写入过滤后 VCF (严格)
  out_strict <- file.path(strict_dir, basename(vcf_path))
  writeLines(c(result$header_lines, result$strict_lines), out_strict)
}

# ============================================================================
# 5. 输出统计表
# ============================================================================
stats_path <- file.path(dirname(primary_dir), "WES_filter_statistics.csv")
write.csv(all_stats, stats_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\n统计表已保存: %s\n", stats_path))

# 打印汇总
cat("\n========== 汇总 ==========\n")
cat(sprintf("总样本数: %d\n", nrow(all_stats)))
cat(sprintf("宽松过滤 - 总 SNV: %d → 过滤后: %d (%.1f%%)\n",
            sum(all_stats$Total.SNV), sum(all_stats$Filtered.SNV),
            100 * sum(all_stats$Filtered.SNV) / max(sum(all_stats$Total.SNV), 1)))
cat(sprintf("严格过滤 - 总 SNV: %d → 过滤后: %d (%.1f%%)\n",
            sum(all_stats$Total.SNV), sum(all_stats$Filtered.strict.SNV),
            100 * sum(all_stats$Filtered.strict.SNV) / max(sum(all_stats$Total.SNV), 1)))

cat("\n前 10 个样本统计:\n")
print(head(all_stats, 10), row.names = FALSE)

cat("\nWES 过滤完成!\n")
