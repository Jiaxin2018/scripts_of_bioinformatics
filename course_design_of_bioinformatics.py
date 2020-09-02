# -*- coding: utf-8 -*-
"""
Created on Tue Sep 1 18:18:35 2020

@author: Jiaxin2018
"""
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from Bio.SeqUtils import CodonUsage as CU
import os
import os.path
import re
#from CAI import RSCU
from matplotlib import pyplot as plt
plt.switch_backend('Agg')
import csv
import warnings
warnings.filterwarnings('ignore')

"""define function & class"""
def Generate_Codon_Table():
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
    standard_table['TAA'] = '*'
    standard_table['TAG'] = '*'
    standard_table['TGA'] = '*'
    return standard_table

class RSCU:
    def __init__(self, cds, translation_table):
        self.cds = cds
        for window_cusp in Sliding_windows(cds,window_size=3,step_size=3):
            continue
        self.cds_codon_num  = window_cusp.window_index/3
        self.codon_fre = window_cusp.frequency
        
        self.aa_codon_dict = {}
        self.aa_fre = {}
        for (key,value) in translation_table.items():
            if value not in self.aa_codon_dict.keys():
               self.aa_codon_dict[value] = []
            self.aa_codon_dict[value].append(key)
        for codon in (self.codon_fre.keys()):
            codon_aa=translation_table[codon] if codon in translation_table.keys() else 'other'
            if codon_aa not in self.aa_fre.keys():
                self.aa_fre[codon_aa] = 0
            self.aa_fre[codon_aa] += self.codon_fre[codon]
   
    def fraction (self):
        RSCU_dict = {
                'fraction':{},
                'frequency':{},
                'RSCU':{},
                }        
        for aa in (self.aa_codon_dict.keys()):
            for codon in self.aa_codon_dict[aa]:
                if codon not in self.codon_fre.keys():
                    self.codon_fre[codon]=0
                    RSCU_dict['fraction'][codon] = 0
                else:
                    RSCU_dict['fraction'][codon] = self.codon_fre[codon]/self.aa_fre[aa]
                RSCU_dict['frequency'][codon] = self.codon_fre[codon]/self.cds_codon_num*1000                    
                RSCU_dict['RSCU'][codon] = RSCU_dict['fraction'][codon]*(len(self.aa_codon_dict[aa]))
        return RSCU_dict
    
    def gc_letter (self):
        gc_dict = {0:0,1:0,2:0}
        for codon in self.codon_fre:
            for i in range(3):
                if re.match("[gc]",codon[i],re.I):
                    gc_dict[i] += self.codon_fre[codon]/self.cds_codon_num
        return gc_dict
    

class Sliding_windows:
    def __init__(self,seq,window_size=100,base='',step_size=1):
        self.seq = seq
        self.base = base
        self.window_size = window_size
        self.step_size = step_size
        self.seq_length = len(seq)
        self.window_index = 0
        self.ratio = 0
        self.frequency={}
        if window_size > self.seq_length:
            print('window_size must be smaller than length of the sequence!')
        
    def __iter__(self):
        return self        
    
    def __next__(self):
        seq = self.seq
        window_index = self.window_index
        window_size = self.window_size
        if (window_index+window_size)<=self.seq_length:
            base = self.base
            window_seq = seq[window_index:window_index+window_size]
            if window_seq not in self.frequency.keys():
                self.frequency[window_seq] = 0
            self.frequency[window_seq]+=1
            ratio = 0
            for letter in base:
                ratio += window_seq.count(letter)/window_size
            self.ratio = ratio
            self.window_index += self.step_size
            return self
        else:
            raise StopIteration
            
"""Get the input gbk file && define the output files"""
gbk_file=''
while not os.path.exists(gbk_file):
    gbk_file=input('Please enter your sequence file in gbk format?\n')
out_dir=input('Please provide a valid output directory[default:./]\n')
if not os.path.exists(out_dir):
    out_dir='./'
out_dir+='/'
basename = os.path.basename(gbk_file)
basename = re.sub(r'.gbk?$', "", basename)
fna_file = out_dir+basename+".CDS.fna"
faa_file = out_dir+basename+".CDS.faa"
geneNC = open(fna_file, "w")
geneAA = open(faa_file, "w")

"""Write the cds, protein and the codon usage"""
for seq in SeqIO.parse(gbk_file, "genbank"):
    print ("\n\n###Dealing with GenBank file of %s, \nOutput: \ngene nc: %s \ngene aa: %s\n\n" % (
            seq.id,
            fna_file,
            faa_file))
    for seq_feature in seq.features :
        if seq_feature.type == "CDS" :
            geneSeq = seq_feature.extract(seq.seq)
            NC_record = SeqRecord(geneSeq, 
                                  id=seq_feature.qualifiers['protein_id'][0]+'_CDS_'+seq.name+str(seq_feature.location), 
                                  description=seq.description)
            AA_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0]), 
                                  id=seq_feature.qualifiers['protein_id'][0]+'_CDS', 
                                  description=seq.description)
            SeqIO.write(NC_record, geneNC,"fasta")
            SeqIO.write(AA_record, geneAA,"fasta")

            rscu_file = out_dir+seq_feature.qualifiers['protein_id'][0]+".codon.csv"
            print("###1st.Caculate RSCU with sliding windows, write to %s\n" % rscu_file)
            rscu_csvfile = open (rscu_file, 'w',newline ='') 
            spamwriter=csv.writer(rscu_csvfile,dialect='excel')
            standard_table = Generate_Codon_Table()
            RSCU_class = RSCU(geneSeq, standard_table)
            gc_letter = RSCU_class.gc_letter()
            cds_gc = (gc_letter[0]+gc_letter[1]+gc_letter[2])/3
            spamwriter.writerow(['#Coding GC',cds_gc])             
            spamwriter.writerow(['#1st letter GC',gc_letter[0]])             
            spamwriter.writerow(['#2nd letter GC',gc_letter[1]])             
            spamwriter.writerow(['#3rd letter GC',gc_letter[2]])
            spamwriter.writerow(['#Codon','AA','Fraction','Frequency','Number'])
            RSCU_dict = RSCU_class.fraction()
            aa_codon_dict = RSCU_class.aa_codon_dict
            for aa in sorted(aa_codon_dict.keys()):
                for codon in sorted( aa_codon_dict[aa], key=lambda x:[['A','C','G','T'].index(y) for y in x] ):
                    spamwriter.writerow([codon,aa,RSCU_dict['fraction'][codon],
                                         RSCU_dict['frequency'][codon],RSCU_class.codon_fre[codon]]) 
            rscu_csvfile.close()
            
            rscu_fig = out_dir+seq_feature.qualifiers['protein_id'][0]+".RSCU.png"            
            print("###2nd.Draw the RSCU graph %s\n" % rscu_fig)
            fig,ax=plt.subplots(dpi=128,figsize=(10,6))
            plt.bar(list(RSCU_dict['RSCU'].keys()),list(RSCU_dict['RSCU'].values()))
            plt.title('RSCU of '+seq_feature.qualifiers['protein_id'][0],fontsize=20) 
            plt.ylabel('RSCU') 
            plt.xlabel('Codons') 
            for label in ax.xaxis.get_ticklabels():                
                label.set_rotation(90)
                label.set_fontsize(5)
            plt.savefig(rscu_fig)

    fre_dict={}
    fre_file = out_dir+seq.id+".compseq.csv"
    fre_csvfile = open (fre_file, 'w',newline ='') 
    print("###3rd.Caculate 1/2 base composition and write to %s\n" % fre_file)
    spamwriter=csv.writer(fre_csvfile,dialect='excel')                        
    for i in [1,2]:
        exp_frq = 1/4**i
        for window_compseq in Sliding_windows(seq.seq,i):
            continue
        sliding_num = window_compseq.window_index
        spamwriter.writerow(['Word size',i])
        spamwriter.writerow(['Total count',sliding_num])
        spamwriter.writerow(['# Word','Obs Count','Obs Frequency','Exp Frequency','Obs/Exp Frequency'])    

        fre = window_compseq.frequency
        for base in sorted( fre.keys(), key=lambda x:[['A','C','G','T'].index(y) for y in x] ):
            fre_num = window_compseq.frequency[base]
            spamwriter.writerow([base,fre_num,fre_num/sliding_num,exp_frq,fre_num/sliding_num/exp_frq])  
    fre_csvfile.close()
            
    cs_dict={}
    print("###4th.Caculate ATGC contents with sliding windows:\n")    
    for base in ['A','C','G','T','AT','GC']:
        cs_file = out_dir+seq.id+"_"+base+".content.csv"
        print("###Caculate %s contents with sliding windows and write to %s\n" % (base,cs_file))    
        cs_dict[base] = []
        with open(cs_file, 'w',newline ='') as csvfile:
            spamwriter=csv.writer(csvfile,dialect='excel') 
            for ratio_window in Sliding_windows(seq.seq,base=base):        
                cs_dict[base].append(ratio_window.ratio)
                spamwriter.writerow([ratio_window.ratio])
    
    ATGC_fig = out_dir+seq.id+".ATGC.png"
    AT_GC_fig = out_dir+seq.id+".AT_GC.png"
    print("###5th.Final part: draw the ATGC content pictures: %s and %s\n" % (ATGC_fig,AT_GC_fig))
    fig=plt.figure(dpi=128,figsize=(10,6))#plot the density graph
    plt.plot(cs_dict['A'],'black',label="A")
    plt.plot(cs_dict['C'],'red',label="C")
    plt.plot(cs_dict['G'],'green',label="G")
    plt.plot(cs_dict['T'],'blue',label="T")
    plt.xlabel('Position',fontsize=16)
    plt.ylabel('A/T/G/C % of '+seq.id,fontsize=16)
    plt.title('Contents of bases of '+seq.id,fontsize=20)
    plt.legend(loc='upper left')
    plt.savefig(ATGC_fig)
    
    fig=plt.figure(dpi=128,figsize=(10,6))
    plt.plot(cs_dict['GC'],'m',label="GC")
    plt.plot(cs_dict['AT'],'yellow',label="AT")
    plt.xlabel('Position',fontsize=16)
    plt.ylabel('AT/GC % of '+seq.id,fontsize=16)
    plt.title('Contents of AT&GC in DNA from '+seq.id,fontsize=20)
    plt.legend(loc='upper left')
    plt.savefig(AT_GC_fig)
    print("###The analysis of %s is finished. Enjoy it!\n\n\n\n" % seq.id)
            
geneNC.close()
geneAA.close()
"""You can also use CodonUsage to caculate the codon usage, however the result seems not correct"""
#myIndex = CU.CodonAdaptationIndex()
#myIndex.generate_index(fna_file)
#myIndex.print_index()
"""You can also use CAI.RSCU to caculate the codon usage"""
#geneseq_list = []
#geneseq_list.append(geneSeq)
#geneseq_codon = RSCU(geneseq_list)#caculate codon usage