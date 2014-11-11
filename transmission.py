#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

import vcf

def read_ped(pedfile):
    '''
    Accepts a path to a PED/FAM file. Returns a dictonary with parents as keys and children as values.
    '''
    d = {}
    with open(pedfile) as f:
        for line in f:
            print line
            print len(line.strip().split())
            fid, iid, father, mother, sex, affected = line.strip().split()
            if father != '.':
                d[father] = iid
            if mother != '.':
                d[mother] = iid
    return d
