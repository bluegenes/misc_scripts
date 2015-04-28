# Tessa Pierce
# 4.24.2015

"""Script to rename FastQ reads that were mistakently named "@No name" & split the interleaved files into 
separate paired files. Grabs barcode from filename.
run via "python rename_fastq.py
"""

import sys, re, os, glob
from optparse import OptionParser

desc = """ Rename Matt's weird RAD-Seq files"""
parser = OptionParser(description = desc)
parser.add_option("--inDir", "--INDirectory", help = "path to all input files" , action="store", type="string", dest="inDir", default = '/Volumes/HARD_DRIVE/ddRAD1/')

(opts, args) = parser.parse_args()

os.chdir(opts.inDir)
files = glob.glob("*.fastq") #list all fastq files in directory

outPath = opts.inDir + '/paired'
if not os.path.exists(outPath):
    os.makedirs(outPath)

for file in files:
    baseName= file.split('.')[0]
    reMatch = re.search('-(BC\d*)', baseName)
    if reMatch is not None:
        barcode = reMatch.groups()[0]
    out1 = os.path.join(outPath, baseName + '_1.fastq')
    out2 = os.path.join(outPath, baseName + '_2.fastq')
    out_fq1 = open(out1, 'w')
    out_fq2 = open(out2, 'w')
    with open(file, 'r') as f:
        counter = 1
        for line in f.readlines():
            if line.startswith('@No name'):
                counter = counter + 1
                if counter % 2 == 0:
                    name = '@' + str(barcode) + '_' + str(counter - 1) + '_' + str(1)
                    out_fq1.write(name + '\n')
                else:
                    name = '@' + str(barcode) + '_' + str(counter - 1) + '_' + str(2)
                    out_fq2.write(name + '\n')
            else:
                if counter % 2 == 0:
                    out_fq1.write(line)
                else:
                    out_fq2.write(line)
        f.close()
    out_fq1.close()
    out_fq2.close()


