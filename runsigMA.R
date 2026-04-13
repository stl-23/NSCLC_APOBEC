# short example 2 simulated panels for testing the tool
library(openxlsx)
devtools::load_all("/Users/stl/Documents/pipelines/sigMA/SigMA/")

# data_dir can be replaced by the character containing 
# the directory defined by the user
#data_dir <- system.file("extdata/vcf_examples/", package="SigMA")
#data_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/11.ovarian_cancer_paired_wilcoxon/add_samples/01.analysis/DNA/vcfs"
data_dir <- "/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/ctDNA/05.sigma/vcfs/"

genomes_matrix <- make_matrix(data_dir, file_type = 'vcf', ref_genome_name = 'hg19')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

genome_file = 'genome_matrix_96.csv'

write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)


message(paste0('96-dimensional matrix is saved in ', genome_file))
message('Running SigMA')

output_file_built_in <- run(genome_file, 
    data = "msk", 
    do_assign = T, 
    do_mva = F,
    tumor_type = 'lung', 
    lite_format = F)

m <- read.csv(output_file_built_in)
lite <- lite_df(m)
write.xlsx(m, file = "lite.signature.all.xlsx")
write.xlsx(lite, file = "lite.signature.xlsx")

