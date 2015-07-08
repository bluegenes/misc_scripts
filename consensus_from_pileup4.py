#Tessa Pierce
#orig. 11.4.2013. This version finished 7.7.2015 to reduce RAM requirements & deal with contigs of 0 coverage & names with spaces (alignment files usually split on the space/provide just the first part of the name in the mpileup file)

# make a new consensus reference based on mpileup, generated via: samtools mpileup -f FASTA BAM --output MPILEUP
import os, sys, re
#from memory_profiler import memory_usage
from optparse import OptionParser
desc = """ A script that reads in a reference file and ONE samtools pileup file from reads that were mapped to that reference, and prints out a new ('alternate') consensus reference using the most commonly found base at each position in the pileup file. By default, uses 'n' if there is insufficient coverage (use -b option to instead use the reference base).  ONLY TAKES A SINGLE PILEUP FILE, so if you have many bam files, merge them first with samtools merge, then make a pileup file (samtools mpileup -f FASTA BAM --output MPILEUP). This version of the script is MUCH more memory efficient, at the slight expense of time. NOTE: it relies on the pileup entries being in the same order as the fasta entries, which will be the case for any true samtools pileup file + the reference that you used to create it.
"""
parser = OptionParser(description = desc)
### Fasta Reference ### 
parser.add_option("-r", "--ref", "--InReference", help = "fasta reference file", action="store", type="string", dest="reference")
#### Piluep Files. Use -p for each file #### 
parser.add_option("-p", "--pileup", "--InPileup", help = "name of the Input pileup file" , action="store", type="string", dest="pileup")
### Input Min Coverege Depth to call a Base ###
parser.add_option("-d", "--depth", "--DepthCoverage", help = "Minimum pileup coverage to call a base" , action="store", type="int", dest="covDepth", default= 8)
### In case of insufficient coverage, use reference base? True/False (Default = False; use 'n')
parser.add_option("-b", "--refBase", "--UseReferenceBase", help = "In cases of insufficient coverage, use reference base instead of n" , action="store_true", default=False, dest="useRefBase")
### Output Files ###
parser.add_option("-o", "--outRef", "--OutAltReference", help = "name of the Output Reference File" , action="store", type="string", dest="out")
parser.add_option("-s", "--outSnp", "--OutSnpTable", help = "name of the Output Snp Table" , action="store", type="string", dest="snp")
# parse args
(opts, args) = parser.parse_args()
####################
# parse bases from pileupLine --> dictionary 'A': #, 'G': #, etc
def parseBases(refAllele, baseInf):
    i = 0
    alleleList = []
    while i < len(baseInf):
	base = baseInf[i]
	if base == '.' or base == ',':# match on forward(.) or reverse(,) strand
	    alleleList.append(refAllele)
	elif base == 'A' or base == 'C' or base == 'G' or base == 'T' or base =='N' or base == 'a' or base == 'c' or base == 'g' or base == 't' or base =='n':
	#base.isalpha(), but then might catch the'M' .. maybe should do: isalpha() and not == 'M'
            alleleList.append(base.upper()) # add this alternate allele (in  uppercase)
	elif base == '+' or base == '-': # insertion or deletion...
	    try: 
		numBases = int(line[i+1])
	        endIndex = i + 1 + numBases
	        alleleList.append(line[i:endIndex].upper()) # need to turn bases into uppercase
	        i = endIndex # advance i by the number of items in the indel
	    except:
		error = True # there's at least one line where - appears and isn't followed by a number or any bp ... not correct pileup format = ignore for now.
	i = i + 1
    alleleD ={}
    for b in alleleList:
        alleleD[b] = alleleList.count(b)
    return alleleD 

def refToAlt_byContig(contigName, refBases, outAltReference, outSNP, contigDt, depth, refBase):    
    refLength = len(refBases)
    newBases = ['n'] * refLength
    refBases = list(refBases)
    insertions = False
    deleteIndices = []
    for coordinate, coordInfo in contigDt.items(): 
	writeSNP_byContig(contigName, coordinate, coordInfo, outSNP) #write snp table; let coordinates be 1-based.
	coordinate = coordinate -1 # python list is 0-indexed, coords from pileup are 1-indexed
	#consensus fasta:
	if int(coordInfo[0]) >= depth: #user-specified depth coverage
	    newAllele = coordInfo[2]
	    if newAllele[0] =='-':  # indel provision
		endDeletion = coordinate + int(newAllele[1])
		deleteIndices.append((coordinate, endDeletion))
		newAllele = ''
	    elif newAllele[0] =='+':
		try:
		    newAllele = newAllele[2:] 
		except:
		    issue = True #just need something here...
		    #import pdb; pdb.set_trace()
		insertions = True
	elif refBase:
	    newAllele = refBases[coordinate]
	else:
	    newAllele = 'n'
	newBases[coordinate] = newAllele 
    for deletion in deleteIndices: 
	i = deletion[0]
	while i < int(deletion[1]):
	    newBases[i] = ''
	    i = i + 1
    bases = ''.join(newBases)
    outAltReference.write('>' + contigName + '_len:' + str(len(bases)) + '\n' + bases + '\n')

# SNP Frequency Table
# contig \t coordinate \t coverage \t refAllele \t mostFreq \t alteredRef? \t bases : freq
def writeSNP_byContig(contigName, coord, coordInfo, outSNP):
    alteredRef = False
    coverage,refAllele,mostFreqBase,freq,baseDt = coordInfo
    if refAllele != mostFreqBase:
	alteredRef = True # can later sort and only print 'True' Items if we want ...
    baseList = []
    for base, hits in baseDt.items():
        baseList.append(base + ': ' + str(hits))
    outSNP.write(contigName + '\t' + str(coord) + '\t' + str(coverage) + '\t' + refAllele + '\t' + mostFreqBase + '\t' + str(alteredRef) + '\t' + '\t'.join(baseList) + '\n')

#main
outAltReference = open(opts.out, 'w')
outSnpTable = open(opts.snp, 'w')
outSnpTable.write("contig" + '\t' + "coordinate" + "\t" + "coverage" + "\t" + "refAllele" + "\t" + "mostFreq" + "\t" + "alteredRef?" + "\t" + "base:freq" + '\n')
contigD = {}
bp = ''
with open(opts.pileup, 'r') as p:
    with open(opts.reference, 'r') as f:
        currentContig = f.readline().strip()[1:].split(' ')[0] # get first contigName
        for line in p:
	    line = line.strip().split('\t')
	    if len(line) == 6: # locations with no coverage only have four items in a pileup line
                cName,coordinate,refAllele,coverage,baseInfo,baseQuals = line
	        baseDt = parseBases(refAllele, baseInfo)
	        mostFreqBase, freq = max(baseDt.items(), key=lambda x:x[1])
	        if cName == currentContig: 
	            contigD[int(coordinate)] = [int(coverage), refAllele, mostFreqBase, freq, baseDt]
	        else:
	            for line in f:
	                if line.startswith('>'):
	   	            refToAlt_byContig(currentContig, bp, outAltReference, outSnpTable, contigD, opts.covDepth, opts.useRefBase) # write newFasta, snptable entries
		            currentContig = line.strip()[1:] 
		            bp = ''
		            contigD = {} #resetdictionary
			    contigD[int(coordinate)] = [int(coverage), refAllele, mostFreqBase, freq, baseDt] # catch first base
		            break
		        else:    
		            bp = bp + line.strip()
	for line in f: #catch bp of last contig
	    bp = bp + line.rstrip()
	refToAlt_byContig(currentContig, bp, outAltReference, outSnpTable, contigD, opts.covDepth, opts.useRefBase) # catch last contig
	
	outAltReference.close()
	outSnpTable.close()
        f.close()
    p.close()
