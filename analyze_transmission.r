#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

options(stringsAsFactors=FALSE) # as it always should be

args = commandArgs(trailingOnly = TRUE)

if (length(args) > 1) {
    cat("This script only accepts 1 command line argument, the text file outputted by transmission.py")
}

path = args[1]
tx = read.table(path,header=FALSE)
colnames(tx) = c("CHROM","POS","REF","ALT","QUAL","FILTER","INDEL","TRANSMITTED")
tx$INDEL = tx$INDEL==1 # change from int to boolean

# overall transmission rate
overall_tx_rate     = mean(as.numeric(tx$TRANSMITTED))
snp_tx_rate         = mean(as.numeric(tx$TRANSMITTED[!tx$INDEL]))
indel_tx_rate       = mean(as.numeric(tx$TRANSMITTED[tx$INDEL]))
pass_tx_rate        = mean(as.numeric(tx$TRANSMITTED[tx$FILTER=="PASS"]))
filtered_tx_rate    = mean(as.numeric(tx$TRANSMITTED[tx$FILTER!="PASS"]))
