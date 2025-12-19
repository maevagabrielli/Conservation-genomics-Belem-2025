#!/usr/bin/python

import sys
import argparse
import textwrap
import gzip
import re

def parse_lines(fin, fou, OUT1, OUT2, SMALL, LARGE):
    totalcount = 0
    for line in fin:
        if line[0] != '#':
            totalcount += 1
            linetab = line.strip().split('\t')
            linesplit=[linetab[x].split(':')[0] for x in range(0,len(linetab))]
            # get maf in outgroup individuals
            #outgroup1
            genoo1=[]
            for i in OUT1:
                geno=linesplit[int(i)]
                if geno=="0/0" or geno=="0|0":
                    genoo1.append(0)
                elif geno=="0/1" or geno=="0|1" or geno=="1/0" or geno=="1|0":
                    genoo1.append(1)
                elif geno=="1/1" or geno=="1|1":
                    genoo1.append(2)
            if len(genoo1)!=0:
                mafOUT1=sum(genoo1)/(2.0*len(genoo1))
            else:
                mafOUT1="NA"
            #outgroup2
            genoo2=[]
            for i in OUT2:
                geno=linesplit[int(i)]
                if geno=="0/0" or geno=="0|0":
                    genoo2.append(0)
                elif geno=="0/1" or geno=="0|1" or geno=="1/0" or geno=="1|0":
                    genoo2.append(1)
                elif geno=="1/1" or geno=="1|1":
                    genoo2.append(2)
            if len(genoo2)!=0:
                mafOUT2=sum(genoo2)/(2.0*len(genoo2))
            else:
                mafOUT2="NA"
            #tot outgroups
            if len(genoo1)!=0 and len(genoo2)!=0:
                OUT=OUT1+OUT2
                genoo=genoo1+genoo2
                mafOUT=sum(genoo)/(2.0*len(genoo))
            else:
                mafOUT="NA"
            # get maf in individuals from the small population
            genos=[]
            for i in SMALL:
                geno=linesplit[int(i)]
                if geno=="0/0" or geno=="0|0":
                    genos.append(0)
                elif geno=="0/1" or geno=="0|1" or geno=="1/0" or geno=="1|0":
                    genos.append(1)
                elif geno=="1/1" or geno=="1|1":
                    genos.append(2)
            if len(genos)!=0:
                mafSMALL=sum(genos)/(2.0*len(genos))
            else:
                mafSMALL="NA"
            # get maf in individuals from the large population
            genol=[]
            for j in LARGE:
                geno=linesplit[int(j)]
                if geno=="0/0" or geno=="0|0":
                    genol.append(0)
                elif geno=="0/1" or geno=="0|1" or geno=="1/0" or geno=="1|0":
                    genol.append(1)
                elif geno=="1/1" or geno=="1|1":
                    genol.append(2)
            if len(genol)!=0:
                mafLARGE=sum(genol)/(2.0*len(genol))
            else:
                mafLARGE="NA"
            ###### NO MISSING DATA: number of genotypes used for the maf calculation = number of genotypes in the individuals of each group (OUT, SMALL, LARGE)
            if len(OUT)==len(genoo) and len(SMALL)==len(genos) and len(LARGE)==len(genol):
                ### ancestral state (of the outgroups) is 0 in the vcf
                if mafOUT==0:
                    if mafSMALL==1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFBXXXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==0 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFOLGXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFOSMXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOBOXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOLGSD'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==0 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOLGSA'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOSMLD'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOSMLA'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    else:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'INV'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                ### ancestral state (of the outgroups) is 1 in the vcf
                elif mafOUT==1:
                    if mafSMALL==0 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFBXXXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFOLGXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==0 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNFOSMXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOBOXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==0 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOLGSD'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOLGSA'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOSMLD'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMNPOSMLA'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    else:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'INV'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                else:
                    fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'NOPOL'+"\t"+"-"+"\t"+"-\n")
            ###### MISSING DATA: polarize when data is present in at least two thirds of the individuals in each group (OUT1, OUT2, SMALL, LARGE)
            elif len(genoo1)>=len(OUT1)*2/3 and len(genoo2)>=len(OUT2)*2/3 and len(genos)>=len(SMALL)*2/3 and len(genol)>=len(LARGE)*2/3:
                ### ancestral state (of the outgroups) is 0 in the vcf
                if mafOUT==0:
                    if mafSMALL==1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFBXXXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==0 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFOLGXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFOSMXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOBOXX'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOLGSD'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL==0 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOLGSA'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOSMLD'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOSMLA'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                    else:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'INV'+"\t"+"0"+"\t"+linesplit[3]+"\n")
                ### ancestral state (of the outgroups) is 1 in the vcf
                elif mafOUT==1:
                    if mafSMALL==0 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFBXXXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFOLGXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==0 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYFOSMXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOBOXX'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==0 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOLGSD'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL==1 and mafLARGE!=0 and mafLARGE!=1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOLGSA'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==0:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOSMLD'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    elif mafSMALL!=0 and mafSMALL!=1 and mafLARGE==1:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'OFMYPOSMLA'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                    else:
                        fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'INV'+"\t"+"1"+"\t"+linesplit[4]+"\n")
                else:
                    fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'NOPOL'+"\t"+"-"+"\t"+"-\n")
            else:
                fou.write(linesplit[0]+"\t"+linesplit[1]+"\t"+'NOPOL'+"\t"+"-"+"\t"+"-\n")

def main():
    parser = argparse.ArgumentParser(prog='polarizeVCFbyOutgroup', description='Write a tag of polarization', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-vcf', help=textwrap.dedent('''specify vcf input file '''))
    parser.add_argument('-output', help='specify output file')
    parser.add_argument('-OUT1', nargs='+', help='vector of individual idx of the OUTGROUP 1')
    parser.add_argument('-OUT2', nargs='+', help='vector of individual idx of the OUTGROUP 2')
    parser.add_argument('-SMALL', nargs='+', help='vector of individual idx of the SMALL pop')
    parser.add_argument('-LARGE', nargs='+', help='vector of individual idx of the LARGE pop')
    args = parser.parse_args()
    print(args)
    if args.vcf is None:
        parser.print_help()
        sys.exit('Please specify vcf input file')
    if args.output is None:
        parser.print_help()
        sys.exit('Please specify output file')
    if args.output.endswith('gz'):
        with gzip.open(args.output, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.OUT1, args.OUT2, args.SMALL, args.LARGE)
            if not args.vcf.endswith('gz'):
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.OUT1, args.OUT2, args.SMALL, args.LARGE)
    if not args.output.endswith('gz'):
        with open(args.output, 'wt') as fou:
            if args.vcf.endswith('gz'):
                with gzip.open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.OUT1, args.OUT2, args.SMALL, args.LARGE)
            if not args.vcf.endswith('gz'):
                with open(args.vcf, 'rt') as fin:
                    parse_lines(fin, fou, args.OUT1, args.OUT2, args.SMALL, args.LARGE)

if __name__ == '__main__':
    main() 
