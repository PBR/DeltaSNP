#!/usr/bin/python
from Bio import SeqIO
import sys
from Bio.SeqIO import FastaIO
import vcf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functools
from matplotlib import colors as mcolors
import seaborn as sns
#pd.set_option('display.max_rows', 2000)

def conjunction(*conditions):
    return functools.reduce(np.logical_and, conditions)

def getMedian(list):
    med=list.median()
    return med

def parseVCF(infile, windowsize):   
    ## chromosome, minPos, maxPos, windowsize):
    d=[]
    
    vcf_reader = vcf.Reader(open(infile, 'r'))
    #print dir(vcf_reader)
    #print vcf_reader.contigs.keys()
    #for contig in vcf_reader.contigs.keys():
        
    samples=vcf_reader.samples
    for sample in samples:
        print sample
    for record in vcf_reader:
#         #print record.is_snp
#       	#if record.is_snp:
       	    if len(record.ALT) == 1 :
                 if (isinstance(record.genotype(samples[0]).data.AO, int)
                     and isinstance(record.genotype(samples[1]).data.AO, int)):
#                     #print "Test: {0}".format(float(record.genotype(samples[0]).data.DP))
       	    	    d.append({'CHROM': record.CHROM, \
                         'POS': int(record.POS), \
                         'REF': record.REF, \
                         'ALT': record.ALT[0], \
                         "DP_{0}".format(samples[0]) : record.genotype(samples[0]).data.DP, \
                         "DP_{0}".format(samples[1]) : record.genotype(samples[1]).data.DP, \
                         "AO_{0}".format(samples[0]) : record.genotype(samples[0]).data.AO, \
                         "AO_{0}".format(samples[1]) : record.genotype(samples[1]).data.AO, \
                         "RATIO_{0}".format(samples[0]) : (record.genotype(samples[0]).data.AO/float(record.genotype(samples[0]).data.DP)), \
                         "RATIO_{0}".format(samples[1]) : (record.genotype(samples[1]).data.AO/float(record.genotype(samples[1]).data.DP)), \
                         "dSNP" : (record.genotype(samples[1]).data.AO/float(record.genotype(samples[1]).data.DP))-(record.genotype(samples[0]).data.AO/float(record.genotype(samples[0]).data.DP))})
# #      	    		
    df=pd.DataFrame(d)
    for contig in df.CHROM.unique():
        sns.set_style("whitegrid")
        colorlist=sns.color_palette("husl", 12)
        chromosome=contig
        df_temp=df.loc[df["CHROM"] == contig]
        # print dfls -
        header = ["CHROM", "POS", "REF", "ALT", "DP1", "DP2", "AO1", "AO2", "RATIO1", "RATIO2", "dSNP"]
        median1=getMedian(df_temp["DP_{0}".format(samples[0])])
        median2=getMedian(df_temp["DP_{0}".format(samples[1])])
        maxThreshold1=3*median1
        maxThreshold2=3*median2
        print "Max Threshold sample {0}={1}".format(samples[0], maxThreshold1)
        print "Max Threshold sample {0}={1}".format(samples[1], maxThreshold2)
        #print df
        ###df_filteredTemp=df[df["DP_{0}".format(samples[0])] < maxThreshold1]
        ###df_filtered=df_filteredTemp[df_filteredTemp["DP_{1}".format(samples[1])] < maxThreshold2]
        #p1 = df["POS"] > int(minPos)
        #p2 = df["POS"] < int(maxPos)
        #p3 = df["REF"] == "SL2.40ch04"
        #df_sub = df[conjunction(p1, p2, p3)]
        df_sub = df_temp #[conjunction(p1, p2)]
        #print "subdataframe:"
        #print df_sub
        c1 = df_sub["DP_{0}".format(samples[0])] < maxThreshold1
        c2 = df_sub["DP_{0}".format(samples[1])] < maxThreshold2
        df_filtered = df_sub[conjunction(c1,c2)]
        #df_filtered.to_csv("{0}_output.csv".format(infile), columns = header, index = None, sep = '\t')
        
        maxThreshold=float(0.6)
        minThreshold=float(-0.6)

        c3 = df_sub["dSNP"] > maxThreshold 
        c4 = df_sub["dSNP"] < minThreshold
        df_filteredMax = df_sub[conjunction(c1,c2,c3)]
        df_filteredMin = df_sub[conjunction(c1,c2,c4)]
        

        df_filtered['SW_mean'] = pd.rolling_mean(df_filtered["dSNP"], int(windowsize), center=True)
        df_filtered['SW_median'] = pd.rolling_median(df_filtered["dSNP"], int(windowsize), center=True)
        t1 = abs(df_filtered['SW_median']) > 0.75
        df_test= df_filtered[conjunction(t1)]
        if len(df_test)>-1:
            df_filtered.to_csv("{0}.SW{2}_{1}.csv".format(infile, chromosome, windowsize), columns = header, index = None, sep = '\t')
            print df_test
            plt.figure(figsize=(20,12))
            plt.subplots_adjust(hspace=0.4)
            plt.subplot(311)
            plt.ylim(-0.1,1.1)
            #plt.xlim(int(minPos), int(maxPos))
            plt.grid(b=True, which='both', color='0.65',linestyle='-')
            plt.axhline(y=0, linewidth=2, linestyle='-', color ='k', alpha=0.5)
            plt.axhline(y=1, linewidth=2, linestyle='-', color ='k', alpha=0.5)
            plt.plot(df_filtered["POS"], df_filtered["RATIO_{0}".format(samples[1])], linestyle='', color=colorlist[1], alpha=0.5, marker='o', label="Ratio {0}".format(samples[1]))
            plt.ylabel("Ratio ALT/DP sample {0}".format(samples[1]))
            plt.title("ALT ratio plot of both pools & delta SNP ratio : {0} ----- reference:{1}, windowsize={2}".format(infile, chromosome, windowsize))
       
            plt.subplot(312)
            plt.ylim(-0.1,1.1)
            #plt.xlim(int(minPos), int(maxPos))
            plt.grid(b=True, which='both', color='0.65',linestyle='-')
            plt.axhline(y=0, linewidth=2, linestyle='-', color ='k', alpha=0.5)
            plt.axhline(y=1, linewidth=2, linestyle='-', color ='k', alpha=0.5)
            plt.plot(df_filtered["POS"], df_filtered["RATIO_{0}".format(samples[0])], linestyle='', color=colorlist[4], alpha=0.5, marker='o', label="Ratio {0}".format(samples[0]))
            plt.ylabel("Ratio ALT/DP sample {0}".format(samples[0]))
        
            maxThreshold=float(0.6)
            minThreshold=float(-0.6)
            plt.subplot(313)
            plt.ylim(-1.1,1.1)
            #plt.xlim(int(minPos), int(maxPos))
            plt.grid(b=True, which='both', color='0.65',linestyle='-')
            plt.axhline(y=maxThreshold, linewidth=1, linestyle='-', color ='#9b9b9b')
            plt.axhline(y=minThreshold, linewidth=1, linestyle='-', color ='#9b9b9b')
            plt.axhspan(maxThreshold, minThreshold, color='#9b9b9b', alpha=0.5, lw=0)
            # #plt.axvspan(42000, 44000, color='c', alpha=0.3, lw=0)
            # #plt.axvspan(0, 13200, color='c', alpha=0.3, lw=0)

        
            plt.plot(df_filtered["POS"], df_filtered["dSNP"], linestyle='', color='#9b9b9b', marker='o', alpha=0.6,)
            plt.plot(df_filteredMax["POS"], df_filteredMax["dSNP"], linestyle='', color=colorlist[6], marker='o', alpha=0.5,)
            plt.plot(df_filteredMin["POS"], df_filteredMin["dSNP"], linestyle='', color=colorlist[6], marker='o', alpha=0.5,)
            plt.plot(df_filtered["POS"], df_filtered["SW_mean"], linestyle='-', color="#e74c3c", marker='', alpha=0.8, linewidth=2.0, label="Mean dSNP, windowsize={0}".format(windowsize))
            plt.plot(df_filtered["POS"], df_filtered["SW_median"], linestyle='-', color="#34495e", marker='', alpha=0.9, linewidth=3.0, label="Median dSNP, windowsize={0}".format(windowsize))
            #plt.plot(df_filtered["POS"], df_filtered["SW_mean"], linestyle='-', color='#990033', marker='', alpha=0.5, linewidth=2.0, label="Mean dSNP, windowsize={0}".format(windowsize))
            #plt.plot(df_filtered["POS"], df_filtered["SW_median"], linestyle='-', color='#ff1aff', marker='', alpha=0.8, linewidth=3.0, label="Median dSNP, windowsize={0}".format(windowsize))
            plt.xlabel('position')
            #plt.ylim(-1,1)
            plt.ylabel("dSNP (Ratio_{0} - Ratio_{1})".format(samples[1], samples[0]))
            #plt.legend()
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=3, borderaxespad=0.)
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, borderaxespad=0.)
            #    plt.legend(handles=[mean, median])
            #plt.annotate("Mean_SlidWin{0}".format(windowsize), (0,0), (0, -20), xycoords='axes fraction', textcoords='offset points', va='top')
            #plt.figtext(0.99, 0.01, "text", horizontalalignment='right')
        
     
            plt.savefig("{0}.SW{2}_{1}.png".format(infile, chromosome, windowsize), format='png')

if __name__ == "__main__":
    parseVCF(sys.argv[1], sys.argv[2])  #, sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
