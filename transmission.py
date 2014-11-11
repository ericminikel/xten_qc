#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

import vcf

class Pedigree:
    def __init__(self,pedfile):
        self.d = {}
        self.founders = []
        with open(pedfile) as f:
            for line in f:
                fid, iid, father, mother, sex, affected = line.strip().split()
                line_dict = {'fid': fid, 'iid': iid, 'father': father, 'mother': mother, 'sex': sex, 'affected': affected}
                self.d[iid] = line_dict
                if father == '.' and mother == '.':
                    self.founders.append(iid)
    def display(self):
        print self.d
    def is_founder(self,id):
        return id in self.founders

