#Tessa Pierce
#started 10.24.2013, finished 11.4.2013

# make a new consensus reference based on pileup

import os, sys, re
from optparse import OptionParser

desc = """ This is a script that reads in a reference file and samtools pileup file from reads that were mapped to that reference, and prints out a new ('alternate') reference using the most commonly found base at each position (from the pileup file).  
"""
parser = OptionParser(description = desc)

### Fasta Reference ### 
parser.add_option("-r", "--ref", "--InReference", help = "fasta reference file", action="store", type="string", dest="reference")

#### Piluep Files. Use -p for each file #### 
parser.add_option("-p", "--pileup", "--InPileup", help = "name of the Input pileup files" , action="append", type="string", dest="pileup", default = [])

### Input Min Coverege Depth to call a Base ###
parser.add_option("-d", "--depth", "--DepthCoverage", help = "Minimum pileup coverage to call a base" , action="store", type="int", dest="covDepth", default= 8)

### Output Files ###
parser.add_option("-o", "--outRef", "--OutAltReference", help = "name of the Output Reference File" , action="store", type="string", dest="out")
parser.add_option("-s", "--outSnp", "--OutSnpTable", help = "name of the Output Snp Table" , action="store", type="string", dest="snp")

(opts, args) = parser.parse_args()

outAltReference = open(opts.out, 'w')
outSnpTable = open(opts.snp, 'w')
contig_name = re.compile('^>(.*)')

# function to turn a fasta file into a dictionary of contigName : bases
def fastaToDict (filename):
    inReference = open(filename, 'r')
    fastaLines = [ x.strip() for x in inReference.readlines() ]
    fastaDt = {}
    bp = ''
    for line in fastaLines:
        if line.startswith('>'):
            if len(bp) > 0:
                fastaDt[name] = [len(bp), bp]
                bp = ''
            name = re.match(contig_name, line).groups()[0] # don't really need the regex .. just want whole line after '>'
        else:
            bp = bp + line
    fastaDt[name] = [len(bp), bp]
    inReference.close()
    return fastaDt

# read in the Pileup file
def readPileup(filename):
    pileupLines =  [ x.strip().split('\t') for x in filename.readlines() ]
    return pileupLines

# parse bases from pileupLine --> dictionary 'A': #, 'G': #, etc
def parseBases(refAllele, line):
    Error = False
    i = 0
    alleleList = []
    while i < len(line):
        base = line[i]
        if base == '.' or base == ',':
            alleleList.append(refAllele)
        elif base == 'A' or base == 'C' or base == 'G' or base == 'T' or base == 'a' or base == 'c' or base == 'g' or base == 't':
            alleleList.append(base.upper()) # put all bases in uppercase
        elif base == '+' or base == '-': # deal with INDELS here
            try:
                numBases = int(line[i+1])
                indelBases = line[i + 1 + numBases].upper()
                endIndex = i + 1 + numBases
                alleleList.append(line[i:endIndex]) # need to turn bases into uppercase
                i = endIndex # advance i by the number of items in the indel
            except:
                error = True # there's at least one line where - appears and isn't followed by a number or any bp ... not correct pileup format, don't have solution = just ignore for now
        i = i + 1
    #alleleDt = {x:alleleList.count(x) for x in alleleList} # create frequency dictionary, same as next 3 lines, works in ipython, but syntax prob here (py 2.6 ->2.5 version issue?)
    alleleD ={}
    for b in alleleList:
        alleleD[b] = alleleList.count(b)
    return alleleD

# function to do the actual work of adding a pileup file to the dictionary
# add to dict of geneName: index of snp, (ref allele, count of ref allele), (alt allele, count of alt allele), etc)
def addPileupToDict(fastaRef, pileupFile, pileupD):
    pileupLines = readPileup(pileupFile)
    for line in pileupLines:
        contigD = {}
        contig_name = line[0]
        coordinate = int(line[1])
        refAllele = line[2]
        coverage = int(line[3])
        baseDt = parseBases(refAllele, line[4])
        if contig_name in pileupD:
            contigD = pileupD.get(contig_name)
            if coordinate in contigD:
                coordinateInfo = contigD.get(coordinate)
                coverage = int(coordinateInfo[0]) + coverage
                oldBaseDt = coordinateInfo[4]
                mergedBaseDt = dict( (n, oldBaseDt.get(n, 0)+ baseDt.get(n, 0)) for n in set(oldBaseDt)|set(baseDt) ) # sum the values for keys that are identical
                baseDt = mergedBaseDt # prob could do this in one step above ...
        mostFreqBase, freq = max(baseDt.items(), key=lambda x:x[1])
        contigD[coordinate] = [coverage, refAllele, mostFreqBase, freq, baseDt]
        pileupD[contig_name] = contigD
    return pileupD


# function to manage addition of each pileup to the dictionary
def parsePileups(fastaRef, pileupList, pileupDt):
    i = 0
    for file in pileupList:
        inF = open(file, 'r')
        addPileupToDict (fastaRef, inF, pileupDt)
        inF.close()
    return pileupDt

#verify INDEL functionality
def pileupDtToFasta(fastaRef, pileupDt):
    newFastaD = {}
    newBases = []
    for contig,info in fastaRef.items():
        refLength = info[0]
        newBases = ['n'] * refLength
        snpInfo = pileupDt.get(contig)
        insertions = False
        deleteIndices = []
        for coordinate, coordInfo in snpInfo.items(): # snpInfo being the dict of coordinates:info for each contig
            coordinate = int(coordinate)
            if int(coordInfo[0]) >= opts.covDepth: #user-specified depth coverage
                newAllele = coordInfo[2]
                if newAllele[0] =='-':  # indel provision
                    endDeletion = coordinate + int(newAllele[1])
                    deleteIndices.append((coordinate, endDeletion))
                    newAllele = ''
                elif newAllele[0] =='+':
                    try:
                        newAllele = newAllele[2:]
                    except:
                        import pdb; pdb.set_trace()
                    insertions = True
            else:
                newAllele = 'n'
            newBases[coordinate - 1] = newAllele # python list is 0-indexed, coords from pileup are 1-indexed
        for deletion in deleteIndices:
            i = deletion[0]
            while i < int(deletion[1]):
                newBases[i] = ''
                i = i + 1
        newFastaD[contig] = ''.join(newBases)
        newBases =[]
    return newFastaD

#main
ref = fastaToDict(opts.reference)
pileupDict = {}
pileupDict = parsePileups(ref, opts.pileup, pileupDict)
newFasta = pileupDtToFasta(ref, pileupDict)

#write out fasta file
for contigName, bases in newFasta.items():
    outAltReference.write('>' + contigName + '_len:' + str(len(bases)) + '\n' + bases + '\n')

# SNP Frequency Table
# contig \t coordinate \t coverage \t refAllele \t mostFreq \t alteredRef? \t bases : freq
for contig, coordDt in pileupDict.items():
    for coord, coordInfo in coordDt.items():
        alteredRef = False
        coverage,refAllele,mostFreqBase,freq,baseDt = coordInfo
        if refAllele != mostFreqBase:
           alteredRef = True # can later sort and only print 'True' Items if we want ...
        baseList = []
        for base, hits in baseDt.items():
            baseList.append(base + ': ' + str(hits))
        outSnpTable.write(contig + '\t' + str(coord) + '\t' + str(coverage) + '\t' + refAllele + '\t' + mostFreqBase + '\t' + str(alteredRef) + '\t' + '\t'.join(baseList) + '\n')

outAltReference.close()
outSnpTable.close()
