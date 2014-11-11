#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

import vcf
import gzip
import argparse

def read_ped(pedfile):
    '''
    Accepts a path to a PED/FAM file. Returns a dictonary with parents as keys and children as values.
    '''
    d = {}
    with open(pedfile) as f:
        for line in f:
            fid, iid, father, mother, sex, affected = line.strip().split()
            if father != '.':
                d[father] = iid
            if mother != '.':
                d[mother] = iid
    return d

def print_transmission(pedfile,vcfpath,biallelic_only=False,min_gq=20,min_dp=10):
    '''
    Loop through a VCF file and print out info about transmission.
    Requires a PED file, and the VCF and PED need to contain *only*
    trios comprised of two parents and a child. No other family
    structures will give correct results.
    '''
    counter = 0
    children = read_ped(pedfile)
    parents = children.keys()
    if vcfpath[-3:] == ".gz": # open .vcf.gz file with gzip.open, otherwise just use open
        openfunc = gzip.open
    else:
        openfunc = open
    vcf_reader = vcf.Reader(openfunc(vcfpath),filename='ignore',compressed=False,strict_whitespace=True) # split only on tab, allow spaces in ids
    for record in vcf_reader: # iterate over every row of VCF
        if biallelic_only and len(record.ALT) > 1:
            continue
        for alt in record.ALT: # for every alt allele at this site
            this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
            this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
            if record.INFO['AC'][this_alt_allele_index] not in [1,2]: # only care about singletons and doubletons
                continue
            else:
                indel = len(alt) != len(record.REF)
                filterstatus = 'PASS' if len(record.FILTER) < 1 else record.FILTER[0]
                zygosity = {} # dictionary for zygosity of each person for this allele. -1 = no call, 0, 1, 2 = num alleles
                for sample in record.samples:
                    if sample['GT'] is None:
                        zygosity[sample.sample] = -1
                    else:
                        zygosity[sample.sample] = map(int,sample['GT'].split("/")).count(this_alt_allele_number)
                n_het_parents = 0
                het_parent = ''
                for parent in parents:
                    n_het_parents += int(zygosity[parent] == 1)
                    if zygosity[parent] == 1:
                        het_parent = parent # save this one for later
                if n_het_parents != 1:
                    continue
                else:
                    # we have a parent/child pair where the parent has the allele. now check if both meet quality thresholds
                    # and whether child has variant
                    child = children[het_parent]
                    parent_call = record.genotype(het_parent)
                    child_call = record.genotype(child)
                    if parent_call.data.GQ < min_gq or parent_call.data.DP < min_dp or child_call.data.GQ < min_gq or child_call.data.DP < min_dp or zygosity[child] == -1:
                        continue # apply quality filters and make sure child is not no-call
                    else:
                        if zygosity[child] == 1:
                            transmitted = True
                        else:
                            assert zygosity[child] == 0, "Child %s has zygosity %s"%(child,zygosity[child])
                            transmitted = False
                        print record.CHROM, record.POS, record.REF, alt, record.QUAL, filterstatus, str(int(indel)), str(int(transmitted))

def main(args):
    print_transmission(args.pedfile,args.vcfpath,args.biallelic_only,args.min_gq,args.min_dp)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
    parser.add_argument('--pedfile', dest='pedfile', action='store', default='',
                    help='path to PED file', type=str)
    parser.add_argument('--vcfpath', dest='vcfpath', action='store', default='',
                    help='path to VCF (optionally gzipped)', type=str)
    parser.add_argument('--biallelic_only', dest='biallelic_only', action='store_true',
                    help='skip multi-allelic sites')
    parser.add_argument('--min_gq', dest='min_gq', action='store', default=20,
                    help='Minimum genotype quality to include', type=int)
    parser.add_argument('--min_dp', dest='min_dp', action='store', default=10,
                    help='Minimum genotype depth to include', type=int)
    args = parser.parse_args()
    main(args)
