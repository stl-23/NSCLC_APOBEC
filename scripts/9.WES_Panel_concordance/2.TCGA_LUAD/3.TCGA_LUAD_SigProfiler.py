#!/usr/bin/env python3
"""
TCGA-LUAD WES 突变特征信号分析 - SigProfilerAssignment
对过滤后的 VCF 进行 SBS 特征搜索

输入: 过滤后的 VCF 文件目录
输出: Activities (签名权重) 文件

用法: python3 3.TCGA_LUAD_SigProfiler.py
"""

import os
import sys

try:
    from SigProfilerAssignment import Analyzer as Analyze
except ImportError:
    print("请先安装 SigProfilerAssignment:")
    print("  pip3 install SigProfilerAssignment")
    sys.exit(1)

# ============================================================================
# 路径设置
# ============================================================================
base_dir = "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/TCGA/LUAD"
vcf_dir = os.path.join(base_dir, "vcfs")

out_base = os.path.join(base_dir, "WES_signatures", "SigProfiler")
os.makedirs(out_base, exist_ok=True)


def run_sigprofiler(vcf_dir):
    """对 VCF 文件运行 SigProfilerAssignment"""
    if not os.path.isdir(vcf_dir):
        print(f"  VCF 目录不存在: {vcf_dir}")
        return

    vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf')]
    if len(vcf_files) == 0:
        print(f"  VCF 目录为空: {vcf_dir}")
        return

    print(f"\n===== SigProfiler =====")
    print(f"  输入 VCF 目录: {vcf_dir}")
    print(f"  VCF 文件数: {len(vcf_files)}")

    Analyze.cosmic_fit(
        samples=vcf_dir,
        output=out_base + "/",
        input_type="vcf",
        context_type="96",
        collapse_to_SBS96=True,
        cosmic_version=3.4,
        exome=False,
        genome_build="GRCh37",
        signature_database=None,
        exclude_signature_subgroups=None,
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=False,
        verbose=True
    )

    # 输出结果文件位置
    activities_file = os.path.join(
        out_base, "Assignment_Solution", "Activities",
        "Assignment_Solution_Activities.txt"
    )
    print(f"  结果: {activities_file}")


# ============================================================================
# 执行分析
# ============================================================================
if __name__ == "__main__":
    print("=" * 50)
    print("TCGA-LUAD SigProfiler 突变特征信号分析")
    print("=" * 50)

    run_sigprofiler(vcf_dir)

    print("\n" + "=" * 50)
    print("SigProfiler 分析完成!")
    print("=" * 50)
