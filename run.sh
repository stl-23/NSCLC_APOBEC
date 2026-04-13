## 1 sigMA sbs analysis
python3 ad_excel2vcf.py --dnafile /Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/ctDNA/ctDNA_snv_filter_c.xlsx --outdir vcfs/
## 2 calculate matrix and cluster
Rscript runsigMA.R

## mutation specturm
#mkdir mutation_spectrum && cd  mutation_spectrum/
#Rscript /Users/stl/Documents/AmoyDX/scripts/draw_plot_sigMA.R ../genome_matrix_96_output_tumortype_ovary_platform_MSK-IMPACT410Genes_cf0.csv
