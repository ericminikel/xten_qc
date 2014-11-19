#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript
setwd('~/d/sci/066xten_qc/')
options(stringsAsFactors=FALSE)
require(sqldf)

# color standards
gcolor = '#FFA824' # aureoline yellow
ecolor = '#0D4F8B' # indigo dye

# get midpoints of a vector
midpoints = function(numericvec) {
  midpointvec = numeric(length(numericvec)-1)
  midpointvec = 0.5*(numericvec[1:(length(numericvec)-1)]+
                       numericvec[2:length(numericvec)])
  return(midpointvec)
}

# exome as gold standard, evaluate how good WGS is
gcp = read.table("wgs.vs.wes.gq20dp10.molt.all.concordance.proportions",skip=1,header=TRUE)
gcp_table = acast(data=subset(gcp, !(Comp_Genotype %in% c("Mismatching_Alleles","MIXED"))),
                                  formula=Eval_Genotype ~ Comp_Genotype,
                                  value.var="Proportion")
gcp_table

gcp_mat = as.matrix(gcp_table)
gcp_mat_rel = gcp_mat[1:3,1:3]/rowSums(gcp_mat[1:3,1:3])
gcp_mat_rel

# WGS as gold standard, evaluate how good exome is
gcp = read.table("wes.vs.wgs.gq20dp10.molt.all.concordance.proportions",skip=1,header=TRUE)
gcp_table = acast(data=subset(gcp, !(Comp_Genotype %in% c("Mismatching_Alleles","MIXED"))),
                  formula=Eval_Genotype ~ Comp_Genotype,
                  value.var="Proportion")
gcp_table 

gcp_mat = as.matrix(gcp_table)
gcp_mat_rel = gcp_mat[1:3,1:3]/rowSums(gcp_mat[1:3,1:3])
gcp_mat_rel


### Depth stuff using DiagnoseTargets

wes_depth = read.table("wes.diagnosetargets.table",header=TRUE)
wgs_depth = read.table("wgs.diagnosetargets.table",header=TRUE)

# check that dimensions are same
dim(wes_depth)
dim(wgs_depth)

depth = wes_depth[,c("CHROM","POS","FILTER","GC")]
depth$WESDP = wes_depth$IDP / 24.0 # divide by 24 individuals
depth$WGSDP = wgs_depth$IDP / 24.0 # divide by 24 individuals

# plot depth by chromosomal position
chrbreaks = c(which(!duplicated(depth$CHROM)),dim(depth)[1])
plot(1:dim(depth)[1],depth$WESDP,,pch='.',col=ecolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',ylim=c(0,200),
     ylab='Mean depth',
     cex.lab=1.4,xlab='Interval',
     main='HiSeq 2000 ICE exomes\nmean depth by Gencode CDS interval')
abline(v=chrbreaks,col='black')
axis(side=1,at=midpoints(chrbreaks),labels=unique(depth$CHROM),lty=0,cex.axis=.8)
axis(side=2,at=c(50,100,150,200),labels=c(50,100,150,200),cex.axis=.8)

plot(1:dim(depth)[1],depth$WGSDP,,pch='.',col=gcolor,
     xaxt='n',xaxs='i',yaxs='i',yaxt='n',ylim=c(0,200),
     ylab='Mean depth',
     cex.lab=1.4,xlab='Interval',
     main='X Ten whole genomes\nmean depth by Gencode CDS interval')
abline(v=chrbreaks,col='black')
axis(side=1,at=midpoints(chrbreaks),labels=unique(depth$CHROM),lty=0,cex.axis=.8)
axis(side=2,at=c(50,100,150,200),labels=c(50,100,150,200),cex.axis=.8)


plot(depth$WGSDP, depth$WESDP, pch='.', xlim=c(0,200),ylim=c(0,200),
     xlab='X Ten whole genome depth', ylab='HiSeq 2000 ICE exome depth',
     main='WGS vs. WES depth. Each point is mean coverage\nof a Gencode CDS interval across 24 people',
     cex.main=.8)

# define problem intervals where WGS depth is < 10 and 
wgs_lowcov = depth$WGSDP < 10
cor.test(depth$WGSDP-depth$WESDP, depth$GC)

