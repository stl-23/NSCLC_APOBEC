#!/Users/stl/opt/miniconda3/bin/python3
"""
WES 突变特征信号分析 - SigProfilerAssignment
对两种过滤标准 (primary/strict) 的 VCF 进行 SBS 特征搜索

输入: 过滤后的 VCF 文件目录
输出: Activities (签名权重) 文件

用法: python3 2.WES_SigProfiler.py
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
primary_vcf_dir = "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.primary"
strict_vcf_dir  = "/Volumes/Elements/Amoydx/10.NSCLC_APOBEC/WES.filter.final.strict"

out_base = "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/major_revise2026.8.19/2.WES_panel_correlation/WES/SigProfiler"

os.makedirs(os.path.join(out_base, "primary"), exist_ok=True)
os.makedirs(os.path.join(out_base, "strict"),  exist_ok=True)


def run_sigprofiler(vcf_dir, label):
    """对一组 VCF 文件运行 SigProfilerAssignment"""
    print(f"\n===== SigProfiler: {label} =====")
    print(f"  输入 VCF 目录: {vcf_dir}")

    output_dir = os.path.join(out_base, label)

    Analyze.cosmic_fit(
        samples=vcf_dir,
        output=output_dir + "/",
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
        output_dir, "Assignment_Solution", "Activities",
        "Assignment_Solution_Activities.txt"
    )
    print(f"  结果: {activities_file}")


# ============================================================================
# 执行分析
# ============================================================================
if __name__ == "__main__":
    print("=" * 50)
    print("WES SigProfiler 突变特征信号分析")
    print("=" * 50)

    run_sigprofiler(primary_vcf_dir, "primary")
    run_sigprofiler(strict_vcf_dir,  "strict")

    print("\n" + "=" * 50)
    print("SigProfiler 分析完成!")
    print("=" * 50)
