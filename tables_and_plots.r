#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript
setwd('~/d/sci/066xten_qc/')
options(stringsAsFactors=FALSE)
require(sqldf)

gcp = read.table("wgs.vs.wes.gq30dp10.molt.all.concordance.proportions",skip=1,header=TRUE)
gcp_table = acast(data=subset(gcp, !(Comp_Genotype %in% c("Mismatching_Alleles","MIXED"))),
                                  formula=Eval_Genotype ~ Comp_Genotype,
                                  value.var="Proportion")
gcp_table

gcp_mat = as.matrix(gcp_table)
gcp_mat_rel = gcp_mat[1:3,1:3]/rowSums(gcp_mat[1:3,1:3])
gcp_mat_rel
