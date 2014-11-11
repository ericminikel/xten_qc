#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

'''
Example usage:
transmission.py --pedfile sample.fam --vcfpath my.vcf.gz --biallelic_only --min_gq 20 --min_dp 10 > transmission.txt

The PED/FAM file and VCF must contain exactly the same set of individuals, and that set must consist
solely of trios with both parents and one child.

Outputs one line per allele that is a singleton among the parents. Each line has this format:
CHROM POS REF ALT QUAL FILTER IS_INDEL TRANSMITTED
'''

import vcf
import gzip
import argparse

def read_ped(pedfile):
    '''
    Accepts a path to a PED/FAM file. Returns a dictonary with parents as keys and children as values.
    Only gives sensible output if there is only one child per parent.
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
    children = read_ped(pedfile) # get dictionary mapping parent to child
    parents = children.keys() # get list of parents
    if vcfpath[-3:] == ".gz": # open .vcf.gz file with gzip.open, otherwise just use open
        openfunc = gzip.open
    else:
        openfunc = open
    vcf_reader = vcf.Reader(openfunc(vcfpath),filename='ignore',compressed=False,strict_whitespace=True) # split only on tab, allow spaces in ids
    for record in vcf_reader: # iterate over every row of VCF
        if biallelic_only and len(record.ALT) > 1: # skip multi-allelics if desired by user
            continue
        for alt in record.ALT: # for every alt allele at this site
            this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
            this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
            if record.INFO['AC'][this_alt_allele_index] not in [1,2]: # only care about singletons and doubletons
                continue
            else:
                indel = len(alt) != len(record.REF) # is this variant an indel?
                filterstatus = 'PASS' if len(record.FILTER) < 1 else record.FILTER[0] # get filtered status
                zygosity = {} # dictionary for zygosity of each person for this allele. -1 = no call, 0, 1, 2 = num alleles
                for sample in record.samples: # loop through and record zygosity of each individual
                    if sample['GT'] is None:
                        zygosity[sample.sample] = -1
                    else:
                        zygosity[sample.sample] = map(int,sample['GT'].split("/")).count(this_alt_allele_number)
                # next we'll count the number of parents that are hets. because we are already limiting to
                # singletons and doubletons according to overall AC, there will 0 to 2 het parents, and
                # we only want cases where there is only 1.
                n_het_parents = 0 # initialize a counter
                het_parent = '' # so we can record the sample ID of the het parent if there is one
                for parent in parents: # loop through parents
                    n_het_parents += int(zygosity[parent] == 1) # add 1 to the counter if the parent is a het
                    if zygosity[parent] == 1: # if we found a het, then...
                        het_parent = parent # ...save this one for later
                if n_het_parents != 1: # if there are 0 or 2 het parents, skip this allele
                    continue
                else: # if there is only 1 het parent, now we can check if variant was transmitted to child.
                    # we have a parent/child pair where the parent has the allele. now check if both meet quality thresholds
                    # and whether child has variant
                    child = children[het_parent] # get sample id of child
                    parent_call = record.genotype(het_parent) # get genotype call details for parent
                    child_call = record.genotype(child) # get genotype call details for child
                    if parent_call.data.GQ < min_gq or parent_call.data.DP < min_dp or child_call.data.GQ < min_gq or child_call.data.DP < min_dp or zygosity[child] == -1:
                        continue # apply quality filters AND make sure child is not a no-call
                    else: # if you reach here, all quality filters were met, we can now ask if variant was transmitted or not.
                        if zygosity[child] == 1: # if child is het, it was transmitted
                            transmitted = True
                        else:
                            # child must not have this allele, so non-transmitted. but we'll just double check with an assertion
                            assert zygosity[child] == 0, "Child %s has zygosity %s"%(child,zygosity[child]) # just double checking
                            transmitted = False
                        # output info about this allele
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
