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

def print_transmission(pedfile,vcfpath,biallelic_only=False):
    '''
    Loop through a VCF file and print out info about transmission
    '''
    counter = 0
    ped = read_ped(pedfile)
    founders = ped.keys()
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
                is_het = {}
                for sample in record.samples:
                    if sample['GT'] is None:
                        is_het[sample.sample] = False
                    else:
                        is_het[sample.sample] = map(int,sample['GT'].split("/")).count(this_alt_allele_number) == 1
                n_het_founders = 0
                het_founder = ''
                for founder in founders:
                    n_het_founders += int(is_het[founder])
                    if is_het[founder]:
                        het_founder = founder # save this one for later
                if n_het_founders != 1:
                    continue
                else:
                    if is_het[ped[het_founder]]:
                        transmitted = True
                    else:
                        transmitted = False
                    print record.CHROM, record.POS, record.REF, alt, record.QUAL, filterstatus, str(int(indel)), str(int(transmitted))

def main(args):
    print_transmission(args.pedfile,args.vcfpath,args.biallelic_only)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
    parser.add_argument('--pedfile', dest='pedfile', action='store', default='',
                    help='path to PED file', type=str)
    parser.add_argument('--vcfpath', dest='vcfpath', action='store', default='',
                    help='path to VCF (optionally gzipped)', type=str)
    parser.add_argument('--biallelic_only', dest='biallelic_only', action='store_true',
                    help='skip multi-allelic sites')
    args = parser.parse_args()
    main(args)
