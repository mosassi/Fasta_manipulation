#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:35:56 2015

@author: Mohamed
"""
from optparse import OptionParser
from Bio import SeqIO


def dic_length(DBfasta):
    dic={}
    for cur_record in SeqIO.parse(DBfasta, "fasta") :
        dic[cur_record.name]=0 
    return dic
    
def make_compare(DBfasta,blastFile,outputFile):
    blast = open(blastFile,"r")
    out = open(outputFile,"w")
    dic = dic_length(DBfasta)
    line = blast.readline()
    while line : 
        column = line.split("\t")
        if column[1] in dic.keys() : 
            dic[column[1]]=1
        line = blast.readline()
    blast.close()
    count = 0
    for key in dic.keys() :
        out.write(key+"\t" + str(dic[key])+"\n")
        count = count + dic[key]
    out.write("Total\t"+str(count))
    
def main() :

    parser = OptionParser(usage="usage: ./%prog -b the_blast_file -d SRD_database -o output output_file ",
                          version="%prog 1.0")
    parser.add_option("-d", "--database",
                      action="store",
                      dest="database",
                      default=False,
                      help="the SRD database (fasta format)")
    parser.add_option("-o", "--output",
                      dest="output",
                      default=False,
                      help="the name of the outputFile")
    parser.add_option("-b", "--blast",
                      dest="blast",
                      default=False,
                      help="Blast of your proteins against the COG database")

                      
                      
    (options, args) = parser.parse_args()
    

    if options.blast==False or  options.output==False or options.database==False :
        parser.error("You didn't put all the required arguments\n")
        
    make_compare(options.database,options.blast,options.output)

    
main()