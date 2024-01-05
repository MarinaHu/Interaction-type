import os
import pandas as pd  
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

intypes = ['wc','cc','bo','so','hb','all']

def heatmap(filename): #read in res_intfreq_dict.dat
    filename = filename.replace('.',' ')
    filename = filename.split()
    filename = filename[0] 

    data = pd.read_fwf(filename + '.dat', sep=' ', header=None, names=intypes)
    data.to_csv(filename + '.csv')


    a = sns.heatmap(data,annot=True,fmt='.0f',yticklabels=True,cmap="crest",vmin=0)
    plt.gcf().set_size_inches(10, 15)
    plt.savefig(filename + '.png',dpi=300)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Heatmap')
    parser.add_argument('-hm', dest='intdat', default = 'None', help='like res_intfreq_dict.dat')
    args = parser.parse_args()

    intdat = args.intdat
    
    heatmap(intdat)
#For the seed and solvent, there are only one row in the heatmap
#For the residues, there are two rows: mc (top) and sc (bottom)

#python3 heatmap.py -hm res_intfreq_dict.dat
