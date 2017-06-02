#!/usr/bin/env python
'''
This script takes a popoolation2 .sync file and runs a
QuasiBinomial (default) or Binomial GLM for each SNP.

The script is modeled after cmh.pl from the Popoolation2 package
(Kofler et al., 2011 Bioinformatics 27: 3435 - 3436).

A correction to add a "1" count to each allele frequency if there are any
frequencies of "0" is implemented.

.sync file count part format needs to be:

pop1            pop-2           pop-3          pop-4          ...\n
1:0:6:155:0:0   2:0:9:239:0:0   1:0:3:127:0:0  0:0:5:155:0:0  ...\n

The script also requires a separate data file for some number
of additional variables. All variables will be in tested simultaneously
in the model. If you want to test variables separately you need to split
the data file.

datafile format should be a comma [,] separated table in
the format (including headers):
_________________
Head | pop,var1,var2
Line1| pop1,35.2,1
Line2| pop2,51.8,-2
Line3| pop3,78.0,3
Line4| pop4,45.1,5
  ...| ...

The first column should be population names or ID
The actual variable names in the header don't matter.

Output format:

QBGLM: [SNP RAW DATA] var1_pval var1_effect var2_pval var2_effect...

'''

import sys
import argparse
from argparse import RawTextHelpFormatter
import os
import csv
import operator

def GetMajorAlleles(cnts,mincnt,minc):
    # checks if a SNP is truly biallelic across all pops
    # returns the major and minor alleles across all pops
    # also checks if there are too many INDELS
    Acnt = 0
    Tcnt = 0
    Ccnt = 0
    Gcnt = 0
    Dcnt = 0
    allele_i = {"A":0, "T":1, "C":2, "G":3,"INDEL":5}
    all_cnt = 0
    for i in cnts:
        alleles = i.split(":")
        Acnt = Acnt + int(alleles[allele_i['A']])
        Tcnt = Tcnt + int(alleles[allele_i['T']])
        Ccnt = Ccnt + int(alleles[allele_i['C']])
        Gcnt = Gcnt + int(alleles[allele_i['G']])
        Dcnt = Dcnt + int(alleles[allele_i['INDEL']])
    # only count an allele if the number of reads >= mincnt
    if Acnt >= mincnt:
        all_cnt = all_cnt + 1
    if Tcnt >= mincnt:
        all_cnt = all_cnt + 1
    if Ccnt >= mincnt:
        all_cnt = all_cnt + 1
    if Gcnt >= mincnt:
        all_cnt = all_cnt + 1
    # if there are not enough alleles, the position is not a SNP
    #print Acnt, Tcnt, Ccnt, Gcnt, all_cnt
    if all_cnt <= 1:
        return "not SNP" #not enough alleles
    counts = [("A",Acnt), ("T",Tcnt), ("C",Ccnt), ("G",Gcnt)]
    counts.sort(key=operator.itemgetter(1))
    # last two alleles are the major and minor alleles
    major_alleles = counts[-2:]
    #print major_alleles, mincnt
    # if Dcnt has >= mincnt then allele is "tainted" by INDELS
    if Dcnt >= mincnt:
        return "not SNP"# SNP tainted by INDELS
    #print Dcnt
    # if the *third* most common allele has count >= mincnt then the SNP
    # is not truly biallelic. There are more than 2 alleles.
    #print counts
    if counts[1][1] > 0:
        #print counts[1][1]
        return "not SNP"# SNP is not biallelic
    # if the *second* most common allele has count < mincnt then the
    # allele is fixed
    #print major_alleles
    for i in major_alleles:
        if i[1] < mincnt:
            return "not SNP"# major allele considered fixed
    return major_alleles


def checkSNP(line,maxc,minc,mincnt):
    # checks the SNP line for coverage
    # 1) each population has at least minc
    # 2) no population has > maxc
    # 3) no indels
    line = [i.replace("\n","") for i in line.split("\t")]
    cnts = line[3:]
    #print cnts #script tester line
    #print maxc, minc #script tester line
    for i in range(0,len(cnts)):
        pop = [int(j) for j in cnts[i].split(":")]
        #print pop, sum(pop), sum(pop) > maxc, sum(pop) < minc #script tester line
        # If the coverage within a population is > maxc or < minc
        # then the SNP is considered invalid
        # coverage is counted across As,Ts,Cs and Gs,
        # Ns and INDELs *not* counted
 	if sum(pop[:-2]) > maxc or sum(pop[:-2]) < minc:
            return "not SNP"#: coverage too high or too low
        
    return cnts

def printRlines(cnts,major_alleles,npops,n,\
                mincnt,line,rescale,scale,zeroes,datfile,test):
    # prints a line that will run a glm in R and a line that will
    # print the results of that glm
    counts = []
    if rescale == 'nr':
    # don't rescale the data
    # don't change zeroes to 1s
        for p in range(0,len(cnts)):
            pop=cnts[p]
            #print pop #script tester line
            allele_i = {"A":0, "T":1, "C":2, "G":3}
            alleles = pop.split(":")
            #print alleles #script tester line
            counts.append(alleles[allele_i[major_alleles[1][0]]])
            counts.append(alleles[allele_i[major_alleles[0][0]]])
            
    if rescale == 'r':
    # rescale the data
        scale=float(scale)
        #print scale
        for p in range(0,len(cnts)):
            pop=cnts[p]
            #print pop #script tester line
            allele_i = {"A":0, "T":1, "C":2, "G":3}
            alleles = pop.split(":")
            asum=int(alleles[allele_i[major_alleles[1][0]]])+\
                  int(alleles[allele_i[major_alleles[0][0]]])
            afreq=float(alleles[allele_i[major_alleles[1][0]]])/float(asum)
            #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
            ac=afreq*scale
            counts.append(str(int(round(ac))))        
            asum=int(alleles[allele_i[major_alleles[1][0]]])+\
                  int(alleles[allele_i[major_alleles[0][0]]])
            afreq=float(alleles[allele_i[major_alleles[0][0]]])/float(asum)
            #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
            ac=afreq*scale
            counts.append(str(int(round(ac))))
    if rescale == 'neff':
    # rescale the data
        #print scale
        for p in range(0,len(cnts)):
            pop=cnts[p]
            #print pop #script tester line
            samples = []
            allele_i = {"A":0, "T":1, "C":2, "G":3}
            #print samples #script tester line
            # Get alleles counts from each sample
            # Deal with Major Allele
            #print sam #script tester line
            alleles = pop.split(":")
            # Calculate the coverage at the SNP (asum)
            asum=int(alleles[allele_i[major_alleles[1][0]]])+\
                  int(alleles[allele_i[major_alleles[0][0]]])
            afreq=float(alleles[allele_i[major_alleles[1][0]]])/float(asum)
            #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
            # Calculate neff and rescale the counts 
            neff=(asum*(n*2)-1)/((n*2)+asum)
            ac=afreq*neff
            counts.append(str(int(round(ac))))        
            # Deal with Minor Allele
            asum=int(alleles[allele_i[major_alleles[1][0]]])+\
                  int(alleles[allele_i[major_alleles[0][0]]])
            afreq=float(alleles[allele_i[major_alleles[0][0]]])/float(asum)
            #print "ASUM: ", asum, "SCALE: ", scale, "ALLELE FREQUENCY: ", afreq # script tester line
            neff=(asum*(n*2)-1)/((n*2)+asum)
            ac=afreq*neff
            counts.append(str(int(round(ac))))        
            #print "NEFF: ", neff, "AC: ", ac, "ASUM: ", asum, "N(n*2): ", n, "(",n*2,")"

    #print counts #script tester line
    counts=','.join(counts)

    # Get additional variables from datafile
    datfilreader=csv.reader(open(datfile,'r'))
    header = datfilreader.next()
    variables = '+'.join(header[1:])
    popnames = []
    for var in header[1:]:
        exec "%s = []" %(var)
    for dat in datfilreader:
        for var in header[1:]:
            indx = header.index(var)
            exec "%s.append(dat[indx])" %(var)
            #print dat[indx]
        # Population names should always be in the first column
        popnames.append(dat[0])
    #print ",".join(popnames)
        
    #print '#Variables: '+variables
    #for var in header[1:]:
    #    exec "print '#'+var+': ',','.join(%s)"%(var)
    print 'matrix<-array(c('+counts+'),'+\
          'dim=c(1,2,'+str(npops)+'),'+\
          'dimnames = list(c("1"),c("A","a"),c("'+\
          '","'.join(popnames)+'")))'
    #print line #script tester line
    #####
    # without warning message:
    # to turn on warnings, comment out the below lines.
    #####
    print 'dat<-get_dat(matrix,zeroes='+zeroes+')'
    for var in header[1:]:
        exec "print 'dat$'+var+'<-c('+','.join(%s)+')'" %(var)
    #print 'print(dat)'
    if test == "qbglm":
        print 'res<-glm(cbind(A_Cnt,Tot_Cnt-A_Cnt)~'+variables+\
              ',family="quasibinomial",data=dat)'
    elif test == "bglm":
        print 'res<-glm(cbind(A_Cnt,Tot_Cnt-A_Cnt)~'+variables+\
              ',family="binomial",data=dat)'
    #print 'print(summary(res))'
    # The treatment and additionaly variable results start after npops rows of the summary
    print 'n_rows<-nrow(summary(res)$coefficients)'
    print 'varrows<-nrow(summary(res)$coefficients)'
    print 'cat(c("'+line.replace('\n','')+'"'+\
          ',summary(res)$coefficients[seq(2,varrows),4]'+\
          ',sep="\\t","\\n"))'
    #####
    # with warning message:
    # to turn on warnings, uncomment the below lines.
    #####

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'],formatter_class=RawTextHelpFormatter)

#MAIN CONTROLS
    parser.add_argument('-filename',
                       help = '.sync file name')

    parser.add_argument('-datfile',
                       help = 'tab separated table of data for each population/line')
    
    parser.add_argument('-npops',
                        help = 'Number of populations.')

    parser.add_argument('-n',
                        help = 'Number of individuals in the pool.')

    parser.add_argument('-mincnt', default = 16,
                        help = 'Minimum count needed to consider a SNP.'+\
                        ' default = 16')

    parser.add_argument('-minc', default = 30,
                        help = 'Minimum coverage needed to consider a SNP.'+\
                        ' default = 30')

    parser.add_argument('-maxc', default = 400,
                        help = 'Maximum coverage needed to consider a SNP.'+\
                        ' default = 400')

    parser.add_argument('-rescale', default = 'nr',
                        help = 'r = rescale to -scale, '+\
                        'nr = no rescaling, neff = rescale to neff')

    parser.add_argument('-scale',default = 'XX',
                        help = 'The new total count for rescaling')

    parser.add_argument('-zeroes', default = 1,
                        help = 'Add (1) or not (0) 1 to each cell if any'+\
                        'cell has count 0')

    parser.add_argument('-test', default = 'qbglm',
                        help = 'qbglm - Quasibinomial GLM, bglm - Binomial GLM')

    args = vars(parser.parse_args())

    try:
        filnam = args['filename']
        npops = int(args['npops'])
        n = int(args['n'])
        mincnt = int(args['mincnt'])
        minc = int(args['minc'])
        maxc = int(args['maxc'])
        rescale = args['rescale']
        scale = args['scale']
        zeroes=str(args['zeroes'])
        datfile=args['datfile']
        test=args['test']

        lines = open(filnam, 'rb')
        
        print 'suppressWarnings(library(methods))'
        print 'source("~/Desktop/Data/RData/RScripts/poolAlFreqTest.R")'
        print '#Parameters: ',"npops =",npops,\
              "mincnt =",mincnt,"min coverege =",minc,\
              "max coverage =",maxc,"rescale =",rescale,\
              "scale =",scale,"zeroes =",zeroes, "test = ",test

        for line in lines:
            SNP = checkSNP(line,maxc,minc,mincnt)
            #print SNP #script tester line
            if SNP != "not SNP":
                #print SNP #script tester line
                major_alleles = GetMajorAlleles(SNP,mincnt,minc)
                if major_alleles != "not SNP":
                    #print major_alleles
                    printRlines(SNP,major_alleles,npops,n,mincnt,\
                                line,rescale,scale,zeroes,datfile,test)
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
