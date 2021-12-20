import time
import os
import pathlib
import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from termcolor import colored

path_to_files = "/mnt/chr12/Data/ZIP/WA2/"
reference_file = "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"

start_time = time.time()
print(colored("\nStarting calculations\n----------------------------------------\n\n", 'red'))


#p = pathlib.Path().absolute().__str__()
#print("path", p)
#os.system('py '+p+'/TADtree.py ' + p +'/firsttry/control_file.txt')

'''Load in the reference arrowhead domain list file'''
f = pd.read_csv(reference_file, sep='\t')
chromosomes_nrs = sorted(f['chr1'].unique())

def load_topdom_output_for_chromosome(filename, chrnr, reference_df=f):
    '''Load the output of the topdom algorithm - a RAWobserved file'''
    
    chr_df = pd.read_csv("%schr%s_100kb.RAWobserved" %(path_to_files, chromosome), sep='\t', header=None, \
                           names = ["contact_start", "contact_end", "raw_observed_score"])
    reference_chr = f[f['chr1']==chromosome]
    reference_chr = reference_chr.rename(columns={'f1':"corner_score", 'f2':'Uvar', 'f3':'Lvar', 'f4':'Usign', 'f5':'Lsign'})
    assert reference_chr['x1'].equals(reference_chr['y1'])
    assert reference_chr['x2'].equals(reference_chr['y2'])
    
    return chr_df, reference_chr

class Domain:
    def __init__(self, start, end, score):
        self.start = start
        self.end = end
        self.score = score
    def __str__(self):
        return ("Domain(start:%d,end:%d)" %(self.start, self.end))
    def concordance(self, otherdomain):
        if otherdomain.end<self.start or otherdomain.start>self.end: return 0
        return min(otherdomain.end, self.end)-max(otherdomain.start, self.start)
    
class DomainList:
    
    def __init__(self, d_list):
        lengths_sum = 0
        self.start_pos = d_list[d_list.columns[0]].iloc[0]
        self.end_pos = d_list[d_list.columns[1]].iloc[0]
        self.domains = [] #list of domains
        
        for idx, row in d_list.iterrows():
            start, end, score = row
            self.domains.append(Domain(start, end, score))
            lengths_sum+=end-start
            if start<self.start_pos: self.start_pos = start
            if end>self.end_pos: self.end_pos = end
        
        self.nr_domains = len(self.domains)
        self.length = self.end_pos-self.start_pos
        self.mean_length = lengths_sum/self.nr_domains
    
    def __str__(self):
        return "DomainList(%d,%d)" %(str(self.start_pos), str(self.end_pos))
    def get_domains(self):
        return self.domains
    def concordance(self, otherlist):
        #sumac = self.domains[0].concordance(otherlist.domains[0])
        sumac=0
        for domain1 in self.domains[:5]:
            for domain2 in otherlist.get_domains():
                sumac += domain1.concordance(domain2)
        return sumac/max(self.length, otherlist.length)

def compare_lists(reflist: DomainList, outlist: DomainList):
    
    nr_diff, nr_proportion = reflist.nr_domains-outlist.nr_domains, reflist.nr_domains/outlist.nr_domains
    print("Number of domains\n\tin reference list: %d\n\tin output list: %d\n\tdifference: %d, %f" \
          %(reflist.nr_domains, outlist.nr_domains, nr_diff, nr_proportion))
    conc = 0#reflist.concordance(outlist)
    print("Mean length of domains:\n\tin reference file: %d\n\tin output file: %d" %(reflist.mean_length, outlist.mean_length))
    print("Overlay in covered distance: %d" %conc)
    return nr_diff, conc


for chromosome in chromosomes_nrs[13:14]:
    chr_df, reference_chr = load_topdom_output_for_chromosome(chromosome, f)
    print("Starting with chromosome %s\n" %chromosome)
    
    reference_domains = DomainList(reference_chr[['x1', 'x2', 'corner_score']])
    output_domains = DomainList(chr_df)
    len_diff, conc = compare_lists(reference_domains, output_domains)
    # print("\n", len_diff, conc)
    
    
    
    


end_time = time.time()
print(colored("\n\n----------------------------------------\nDuration of execution: %f seconds\n" %(end_time-start_time), 'red'))