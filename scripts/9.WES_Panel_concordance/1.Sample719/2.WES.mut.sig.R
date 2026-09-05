#!/usr/bin/env Rscript
# ============================================================================
# WES 突变特征信号分析 - deconstructSigs & MutationalPatterns
# 对两种过滤标准 (primary/strict) 的 VCF 分别进行 3 种方法的特征搜索
#
# 方法 1: deconstructSigs (R)
# 方法 2: MutationalPatterns - regular & strict (R)
# 方法 3: SigProfilerAssignment (见 2.WES_SigProfiler.py)
#
# 参考基因组: hg19 (BSgenome.Hsapiens.UCSC.hg19)
# ============================================================================

suppressPackageStartupMessages({
  library(deconstructSigs)
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
})

# ----------------------------------------------------------------------------
# 修补 MutationalPatterns 的 fit_to_signatures_strict bug:
# 内部函数 .plot_sim_decay 中 S4Vectors::isEmpty() 不能处理 NULL 对象
# 用 vapply + is.null 替代
# ----------------------------------------------------------------------------
.fixed_plot_sim_decay <- function(sims, removed_sigs, max_delta,
                                   method = c("backwards", "best_subset")) {
  Removed_signatures <- Cosine_similarity <- NULL
  method <- match.arg(method)
  sims_empty <- vapply(sims, is.null, logical(1))
  sims <- unlist(sims[!sims_empty])
  rs_empty <- vapply(removed_sigs, is.null, logical(1))
  removed_sigs <- unlist(removed_sigs[!rs_empty])
  tb <- tibble::tibble(
    Cosine_similarity = sims,
    Removed_signatures = factor(removed_sigs, levels = removed_sigs)
  )
  sims_l <- length(sims)
  col <- rep("low_delta", sims_l)
  if (sims_l > 1) {
    final_delta <- sims[sims_l - 1] - sims[sims_l]
    if (final_delta > max_delta) col[sims_l] <- "high_delta"
  }
  if (method == "backwards") {
    xlab <- "Removed signatures"
    my_theme <- ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
      text = ggplot2::element_text(size = 12)
    )
  } else {
    xlab <- "Nr. signatures used"
    my_theme <- ggplot2::theme(text = ggplot2::element_text(size = 12))
  }
  ggplot2::ggplot(tb, ggplot2::aes(Removed_signatures, Cosine_similarity, fill = col)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(limits = c("low_delta", "high_delta"),
                               values = c("grey", "red"), guide = "none") +
    ggplot2::labs(x = xlab,
                  y = paste0("Cosine similarity (max delta: ", max_delta, ")")) +
    ggplot2::theme_classic() + my_theme
}
assignInNamespace(".plot_sim_decay", .fixed_plot_sim_decay, ns = "MutationalPatterns")

# ============================================================================
# 1. 路径设置
# ============================================================================
primary_vcf_dir <- "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.primary"
strict_vcf_dir  <- "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.strict"

out_base <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/WES"
dir.create(file.path(out_base, "deconstructSigs"),            recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "MutationalPatterns"),         recursive = TRUE, showWarnings = FALSE)

bsg <- BSgenome.Hsapiens.UCSC.hg19

# ============================================================================
# 2. 通用函数: 从 VCF 提取 SNV 突变
# ============================================================================
parse_vcf_snvs <- function(vcf_path) {
  all_lines   <- readLines(vcf_path, warn = FALSE)
  header_idx  <- grep("^#CHROM", all_lines)[1]
  data_lines  <- all_lines[(header_idx + 1):length(all_lines)]
  data_lines  <- data_lines[nchar(trimws(data_lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  cols        <- strsplit(all_lines[header_idx], "\t")[[1]]
  tumor_idx   <- length(cols) - 1
  sample_name <- cols[tumor_idx]

  fields_list <- strsplit(data_lines, "\t")
  chroms  <- vapply(fields_list, `[`, character(1), 1)
  pos_str <- vapply(fields_list, `[`, character(1), 2)
  refs    <- vapply(fields_list, `[`, character(1), 4)
  alts    <- vapply(fields_list, `[`, character(1), 5)

  # 仅保留 SNV (REF/ALT 均为单碱基)
  keep <- nchar(refs) == 1 & nchar(alts) == 1 &
          refs %in% c("A","C","G","T") & alts %in% c("A","C","G","T")
  if (!any(keep)) return(NULL)

  # 标准化染色体: 添加 chr 前缀
  ch <- chroms[keep]
  if (!any(grepl("^chr", ch))) ch <- paste0("chr", ch)

  data.frame(
    Sample = sample_name,
    chr    = ch,
    pos    = as.numeric(pos_str[keep]),
    ref    = refs[keep],
    alt    = alts[keep],
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# 3. deconstructSigs 分析
# ============================================================================
run_deconstructSigs <- function(vcf_dir, label) {
  cat(sprintf("\n===== deconstructSigs: %s =====\n", label))

  vcf_files <- list.files(vcf_dir, pattern = "[.]vcf$", full.names = TRUE)
  cat(sprintf("  读取 %d 个 VCF 文件...\n", length(vcf_files)))

  # 合并所有样本突变
  mut_list <- lapply(vcf_files, parse_vcf_snvs)
  mut_df   <- do.call(rbind, mut_list[!vapply(mut_list, is.null, logical(1))])
  cat(sprintf("  共 %d 个样本, %d 个 SNV\n",
              length(unique(mut_df$Sample)), nrow(mut_df)))

  # 转换为 deconstructSigs 输入格式
  # mut.to.sigs.input 返回: 样本(行) x 96context(列) 的 matrix
  # whichSignatures 需要同样方向的 data frame
  sigs.input <- mut.to.sigs.input(
    mut.ref    = mut_df,
    sample.id  = "Sample",
    chr        = "chr",
    pos        = "pos",
    ref        = "ref",
    alt        = "alt",
    bsg        = bsg
  )
  sigs.input <- as.data.frame(sigs.input)  # matrix → data frame

  # 逐样本计算特征权重
  samples <- rownames(sigs.input)
  cat(sprintf("  计算 %d 个样本的特征权重...\n", length(samples)))

  weights_list <- lapply(seq_along(samples), function(i) {
    s <- samples[i]
    if (i %% 100 == 0) cat(sprintf("    [%d/%d] %s\n", i, length(samples), s))
    result <- whichSignatures(
      tumor.ref        = sigs.input[s, , drop = FALSE],
      signatures.ref   = signatures.exome.cosmic.v3.may2019,
      sample.id        = s,
      contexts.needed  = TRUE,
      tri.counts.method = "exome2genome"
    )
    result$weights
  })

  w <- do.call(rbind, weights_list)
  cat(sprintf("  权重矩阵: %d 样本 x %d 签名\n", nrow(w), ncol(w)))

  # 保存
  out_dir <- file.path(out_base, "deconstructSigs")
  save(w, file = file.path(out_dir, sprintf("deconstructSigs_%s.Rdata", label)))
  write.csv(w, file = file.path(out_dir, sprintf("deconstructSigs_%s_weights.csv", label)))
  cat(sprintf("  结果已保存: %s/\n", out_dir))
}

# ============================================================================
# 4. MutationalPatterns 分析
# ============================================================================
run_MutationalPatterns <- function(vcf_dir, label) {
  cat(sprintf("\n===== MutationalPatterns: %s =====\n", label))

  vcf_files <- list.files(vcf_dir, pattern = "[.]vcf$", full.names = TRUE)
  cat(sprintf("  读取 %d 个 VCF 文件为 GRanges...\n", length(vcf_files)))

  # # 逐个 VCF 构建 GRanges
  # gr_list <- list()
  # for (i in seq_along(vcf_files)) {
  #   vcf_path <- vcf_files[i]
  #   if (i %% 100 == 0) cat(sprintf("    [%d/%d]\n", i, length(vcf_files)))
  # 
  #   mut <- parse_vcf_snvs(vcf_path)
  #   if (is.null(mut) || nrow(mut) == 0) next
  # 
  #   sample_name <- mut$Sample[1]
  #   gr <- GRanges(seqnames = mut$chr, ranges = IRanges(start = mut$pos, width = 1),
  #                 ALT = mut$alt, REF = mut$ref)
  #   # 仅保留标准染色体
  #   std_chr <- paste0("chr", c(as.character(1:22), "X", "Y"))
  #   gr <- gr[seqnames(gr) %in% std_chr]
  #   if (length(gr) == 0) next
  #   # 先设置 seqlevels 再赋值 seqinfo
  #   seqlevels(gr) <- std_chr
  #   seqinfo(gr) <- seqinfo(bsg)[std_chr]
  #   gr <- sortSeqlevels(gr)
  #   names(gr) <- NULL
  #   gr_list[[sample_name]] <- gr
  # }
  # 
  # if (length(gr_list) == 0) {
  #   cat("  未找到突变，跳过\n")
  #   return(invisible(NULL))
  # }
  # cat(sprintf("  成功构建 %d 个样本的 GRangesList\n", length(gr_list)))
  # 
  # # 转换为 GRangesList
  # grl <- GRangesList(gr_list)

  sample_names <- tools::file_path_sans_ext(
    basename(vcf_files)
  )
  
  grl <- read_vcfs_as_granges(
    vcf_files = vcf_files,
    sample_names = sample_names,
    genome = bsg,
    type = "snv",
    change_seqnames = TRUE,
    remove_duplicate_variants = TRUE
  )
  
  # 提取三核苷酸上下文并构建 96 通道矩阵
  cat("  提取三核苷酸上下文...\n")
  mm <- mut_matrix(vcf_list = grl, ref_genome = bsg)
  cat(sprintf("  突变矩阵维度: %d x %d\n", nrow(mm), ncol(mm)))

  # # 拟合 COSMIC 签名 (使用 MutationalPatterns 内置的已知签名)
  # sigs_for_mp <- get_known_signatures()
  # 
  # cat("  拟合 COSMIC 签名 (regular)...\n")
  # fit_reg <- fit_to_signatures(mm, sigs_for_mp)
  # 
  # # strict 拟合: 根据 regular 结果筛选签名后用 fit_to_signatures_strict
  # # 策略: 过滤掉在 >70% 样本中贡献为 0 的签名，保留剩余签名
  # # 如果 SBS2/SBS13 不在其中则强制加入 (APOBEC 相关签名)
  # cat("  筛选签名用于 strict 拟合...\n")
  # reg_contrib <- fit_reg$contribution  # n_sigs x n_samples
  # n_samples <- ncol(reg_contrib)
  # zero_frac <- rowMeans(reg_contrib == 0)  # 每个签名在多少比例样本中为 0
  # 
  # # 从 70% 开始，如果签名仍太多则逐步提高到 80%、90%
  # threshold <- 0.70
  # selected_sigs_names <- rownames(reg_contrib)[zero_frac < threshold]
  # if (length(selected_sigs_names) > 30) {
  #   threshold <- 0.80
  #   selected_sigs_names <- rownames(reg_contrib)[zero_frac < threshold]
  # }
  # if (length(selected_sigs_names) > 30) {
  #   threshold <- 0.90
  #   selected_sigs_names <- rownames(reg_contrib)[zero_frac < threshold]
  # }
  # 
  # # 强制加入 SBS2 和 SBS13 (APOBEC 相关签名)
  # for (s in c("SBS2", "SBS13")) {
  #   if (!s %in% selected_sigs_names) {
  #     selected_sigs_names <- c(selected_sigs_names, s)
  #     cat(sprintf("    强制加入 %s\n", s))
  #   }
  # }
  # cat(sprintf("    阈值: %.0f%%, 筛选出 %d 个签名\n", threshold * 100, length(selected_sigs_names)))
  # 
  # selected_sigs <- sigs_for_mp[, selected_sigs_names, drop = FALSE]

  # ============================================================
  # 1. Regular fitting
  # ============================================================
  
  sigs_for_mp <- get_known_signatures()
  
  cat("  拟合 COSMIC 签名 (regular)...\n")
  
  fit_reg <- fit_to_signatures(
    mm,
    sigs_for_mp
  )
  
  reg_contrib <- fit_reg$contribution
  
  
  # ============================================================
  # 2. Mutation burden QC
  # ============================================================
  
  mut_burden <- colSums(mm)
  
  cat("  Mutation burden summary:\n")
  print(summary(mut_burden))
  
  cat(sprintf(
    "  Samples with <20 mutations: %d/%d\n",
    sum(mut_burden < 20),
    length(mut_burden)
  ))
  
  cat(sprintf(
    "  Samples with <50 mutations: %d/%d\n",
    sum(mut_burden < 50),
    length(mut_burden)
  ))
  
  
  # ============================================================
  # 3. Signature prevalence
  # ============================================================
  
  # 不建议使用 == 0
  zero_frac <- rowMeans(reg_contrib < 1e-6)
  
  signature_prevalence <- 1 - zero_frac
  
  signature_info <- data.frame(
    signature = rownames(reg_contrib),
    prevalence = signature_prevalence,
    median_contribution = apply(
      reg_contrib,
      1,
      median
    )
  )
  
  print(
    signature_info[
      order(-signature_info$prevalence),
    ]
  )
  
  
  # ============================================================
  # 4. Select signatures
  # ============================================================
  
  selected_sigs_names <- signature_info$signature[
    signature_info$prevalence >= 0.10
  ]
  
  # 强制保留 APOBEC signatures
  selected_sigs_names <- union(
    selected_sigs_names,
    c("SBS2", "SBS13")
  )
  
  cat(sprintf(
    "  Selected %d signatures:\n",
    length(selected_sigs_names)
  ))
  
  print(selected_sigs_names)
  
  # ============================================================
  # 5. Remove signatures with zero total weight
  # ============================================================
  
  selected_sigs <- sigs_for_mp[
    ,
    selected_sigs_names,
    drop = FALSE
  ]
  
  valid_sig <- colSums(selected_sigs) > 0
  
  selected_sigs <- selected_sigs[
    ,
    valid_sig,
    drop = FALSE
  ]
  
  selected_sigs_names <- colnames(selected_sigs)
  
  cat("  拟合 COSMIC 签名 (strict)...\n")
  fit_strict <- fit_to_signatures_strict(
    mm, selected_sigs,
    max_delta = 0.004,
    method = "backwards"
  )

  # 保存
  out_dir <- file.path(out_base, "MutationalPatterns")
  write.csv(fit_reg$contribution,        file = file.path(out_dir, sprintf("MutPatterns_%s_regular_weights.csv", label)))
  write.csv(fit_strict$fit_res$contribution, file = file.path(out_dir, sprintf("MutPatterns_%s_strict_weights.csv", label)))
  cat(sprintf("  Regular: 60 个签名全拟合\n"))
  cat(sprintf("  Strict:  %d 个签名输入)\n", length(selected_sigs_names)))
  cat(sprintf("  结果已保存: %s/\n", out_dir))
}

# ============================================================================
# 5. 执行分析 (两种过滤标准)
# ============================================================================
cat("========================================\n")
cat("WES 突变特征信号分析\n")
cat("========================================\n")

# --- Primary (宽松) ---
run_deconstructSigs(primary_vcf_dir, "primary")
run_MutationalPatterns(primary_vcf_dir, "primary")

# --- Strict (严格) ---
run_deconstructSigs(strict_vcf_dir, "strict")
run_MutationalPatterns(strict_vcf_dir, "strict")

cat("\n========================================\n")
cat("R 部分分析完成 (deconstructSigs + MutationalPatterns)\n")
cat("SigProfiler 请运行: python3 2.WES_SigProfiler.py\n")
cat("========================================\n")
