#!/usr/bin/env python
# encoding: utf-8
#
# Robin Kageyama - Ansel Lab 2015
#

"""
v1.1 added clip_tools.bedheader
v1.2 added clip_tools.getfasta and .reversecomplement
v1.3 changed BetweenDict. No longer lookup with __get__:  bd[chrom,strand,position]. 
     have to explicitly call bd.lookup(chrom,strand,position). Much faster.

Set of useful tools for CLIP analysis

CLASSES
   
    BetweenDict : Lookups for geneID by location
        bD = BetweenDict("location/to/bedfile.bed")
         call with bD("chr1","+",100000,200000)
    BedFile
        b = BedFile("bedfile.bed")
    BedLine
        bl = BedLine([chrom,start,stop,name,0,strand]): bl.line bl.geneid bl.location
    Depth : Contains depths for a region in a *stranded* depth dictionary {"+":dd}
        d = Depth(depthDictionary,chr,start,stop,strand)

FUNCTIONS

    MIRNAs:

    parse_TS_miRNA_names(miRNAstring) - > [List of miRNAs]
        Converts lengthy targetscan famlily names into single miRNA names
    miRbase_reads(miNRA_string) 
        Returns the number of reads from miRbase recorded.
    getSeed(mature_miRNA_sequence)
        Returns a 8mer 7mer1a, 7merm8 and 6mer TARGET seed sequences.

    Gene Info:
    
    bedheader(**kargs) - > "String"
        writes a header for a bedfile
    depths_to_dictionary(depth_file)
        converts a depthfile into a depth dictionary. not strand specific.
    getfasta(chrom, start, stop, strand='+', genome='/Users/DarthRNA/Documents/Adam/genomes/mm10.fa')
        returns a sequence for a given location
    max_depths(bed_file, depthDictionary, conditions=[''])
        Finds the maximum depths from the regions of a bedfile.
    reversecomplement(sequence)
        returns uppercase RC of sequence string.
    easylocation(location)
        converts an igv location "chr1:100-1,000" to chrom,start,stop (integers)


"""

__version__ = "v.1.2"

defaultgenome = "Data/genomes/mm10.fa"
defaultUTR = "Data/genomes/refseq_3utr.bed"

import clip_tools
import csv
import subprocess
from subprocess import Popen, PIPE, STDOUT


class Depth():
    '''
    Makes a depth object
    This object is designed to take in a depthdictionary formatted: 
        (x = {"+":cliptools.depths_to_dictionary(watson_depths),
             "-":cliptools.depths_to_dictionary(crick_depths)}
            )
    AND, a bedline or chr,start,stop,strand

    USAGE: d = Depth(depthDictionary,chr,start,stop,strand)
    OR   : d = Depth(depthDictionary,bedline=[chrom,start,stop,geneid,0,strand])

    HOLDS: depths between start and stop

    To Get info out: d[0]  - Returns -  <int depth at position 0>
                     d.locate(genomicposition) - Returns - <int depth at position 0>
                     d.max() - Returns - max depth in region.

    Iterate through d to get all the depths. 

    '''

    def __init__(self, depthfile, *args, **kargs):
        """Takes either (self,depthfile,bedline=<bedline>) or the args:
        chrom,start,stop,strand
        Depth Dictionary must be in format:
        {"+":{"Chr1":{Position:Depth}}}
        
        Can call max, sum, etc...
        
        """
        from collections import Counter
        if 'bedline' in kargs:  # Read in only a bedline
            self.chrom, start, stop, self.geneid, self.zero, self.strand = kargs['bedline'][0:6]
        if len(args) == 4:  # If chrom,start,stop,strand format.
            self.chrom, start, stop, self.strand = args
        if 'geneid' in kargs:  # If a genename was given
            self.geneid = kargs['geneid']

        self.start = str(start)
        self.stop = str(stop)
        self.int_start = int(self.start)
        self.int_stop = int(self.stop)

        self.depths = []
        for n in xrange(self.int_start, self.int_stop):
            try:
                self.depths.append(float(depthfile[self.strand][self.chrom][str(n)]))
            except KeyError:
                self.depths.append(0)

    def __getitem__(self, key):
        return self.depths[key]

    def __len__(self):
        return len(self.depths)

    def locate(self, key):
        '''
        Query a chromosomal position, returns the depth
        example: depth.locate(234342389)
        '''
        return self.depths[int(key) - self.int_start]


class BedLine():
    '''A single line of a bed file, seperated into components.
        Requires a strict order 
        [chrom,start,stop,name,0,strand]

        bl = BedLine([chrom,start,stop,name,0,strand])
        bl.line - returns - [chrom,start,stop,name,0,strand]
        bl.location - returns "chrX:100-200"

    '''

    def __init__(self, line, extra_features=list()):
        self.bed_length = 6
        self.line = line[0:6]
        self.extra_features = extra_features
        if self.bed_length > 6:
            self.extra_features = line[6:]
        self.parse(self.line)

    def igv(self):
        """returns 'chr:1-100' for easy copy paste to igv"""
        return "%s:%s-%s" % (self.chrom, str(self.start), str(self.stop))

    def scan(self):
        '''Return a generator yielding the numbers 
        beginning with start and ending before stop'''
        for n in range(self.start, self.stop):
            yield n

    def parse(self, bedline):
        '''Returns an index of features in the bedfile
        Takes in either a list of features, or a literal line
        from a bedfile (including tabs).
        '''
        if "\t" in bedline:  # Convert Tabbed to List
            print("Tabbed!")
            bedline = bedline.split("\t")
        if self.check_bed_line(bedline):
            self.chrom = bedline[0]
            self.start = int(bedline[1])
            self.stop = int(bedline[2])
            self.geneid = bedline[3]
            self.zero = bedline[4]
            self.strand = bedline[5]
            self.location = "%s:%s-%s" % (self.chrom, str(self.start), str(self.stop))

    def check_bed_line(self, bedline):
        '''Raises warnings if the bedfile is not what is expected.
        Checks for length(fatal, the 0 column (nonfatal) and strand column (nonfatal)
        '''
        if len(bedline) < self.bed_length:
            raise ValueError("Incomplete bedfile. Recieved %s lines instead of %s." % (len(bedline), self.bed_length))
            return False
        if len(bedline) > self.bed_length:
            if not len(bedline) == len(self.extra_features) + 6:
                raise ValueError('Error: Abnormal Length. Specify self.extra_features if necessary')
        if not int(bedline[4]) == 0:
            print("Warning: Column 5 is not 0. Make sure to check format.")
            return True
        if not (str(bedline[5]) == '+' or str(bedline[5]) == '-'):
            print("Warning: Column 6 is not strand. Make sure to check format")
            return True
        return True


class BedFile():
    ''' 
    Holds a list of Bedline objects in self.bedlines  
    initialize with: 
    
    b = BedFile("bedfile.bed") OR
    b = BedFile([chr,start,stop,geneid,score,strand]) OR
    b = BedFile([[chr,start,stop,geneid,score,strand],[]])

    Can also specify header=True if the first line of the file is a header.

    Call with:
    b[0] > returns first BedObject of the file.
        can specify b[0].chrom or b[0].start to get individual properties
    b[0].line > returns first line of the file: ['chr1', '3214499', '3214505', 'miR-155', '0', '-']
    b.getgene("Icos") > returns a list of BedLines that have that geneid
        specify r=True to get the literal bedlines instead of objects.
    b.save("filename") > to save the object as a bedfile.
        specify header=True if you want to write a header
        specify overwrite = True to ignore overwrite prompt

    '''

    def __init__(self, bed_file_in="", lines=[], header=True):
        '''pass it either,
        > a bedfile or
        > some number of lines [[chr,start,stop,geneid,score,strand],[]]
        which it will convert into BedLine objects and store in a list (bedlines)'''
        import os
        self.bedlines = []
        if len(lines) > 1:
            for line in lines:
                self.bedlines.append(BedLine(line))
        elif len(lines) == 1:
            self.bedlines.append(BedLine(lines))
        elif len(bed_file_in) > 0:
            import csv
            self.bed_file_in = bed_file_in
            if os.path.isfile(bed_file_in):
                with open(bed_file_in, "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    if header == True:
                        self.header = reader.next()
                    for line in reader:
                        self.bedlines.append(BedLine(line))
            else:
                raise "Error: bed file not found"
        else:
            raise "InputError: requires (bedfile) or ([bedline,bedline])"
        self.indexbygene()

    def __len__(self):
        return len(self.bedlines)

    def __getitem__(self, key):
        return self.bedlines[key]

    def addline(self, line):
        '''Adds a line as a BedLine'''
        self.bedlines.append(BedLine(line))

    def removeline(self, position_to_cut=-1, silent=False):
        '''removes BedLine at specified position, or the last line if none
            provided. Returns the removed Bedline '''
        if silent:
            self.bedlines.pop(position_to_cut)
        else:
            return self.bedlines.pop(position_to_cut)

    def getgene(self, geneid, r=False):
        '''Return a list of bedlines that match the given geneid.
        Rebuilds the geneindex if new Bedlines have been added.
        Set r=True to return literal bed lines instead of Bedline objects
        '''
        lines = []
        if self.indexflag != len(self.bedlines):  # check to make sure the index is updated
            print("rebuilding index")
            self.indexbygene()
        if r:
            templine = []
            for i in self.geneindex[geneid]:
                templine.append(i.line)
            return templine
        return self.geneindex[geneid]

    def indexbygene(self):
        '''creates a dict of lines that can be looked up by geneid'''
        from collections import defaultdict
        self.indexflag = 0  # Sets how many lines are indexed
        self.geneindex = defaultdict(list)
        for line in self.bedlines:
            self.geneindex[line.geneid].append(line)
            self.indexflag += 1

    def save(self, savefile, overwrite=False, header=True):
        import os
        import csv
        if (
                        overwrite == False and
                    os.path.isfile(savefile)
        ):
            p = str(raw_input("The file %s exists already. Overwrite? (y/n)" % savefile))
            if p == "y":
                print(p)
            elif p == "n":
                raise ValueError("Terminated by user to prevent overwrite of %s" % savefile)
            else:
                raise ValueError("Input Error")
        with open(savefile, 'w') as f:
            writer = csv.writer(f, delimiter="\t")
            if header == True:
                writer.writerow(self.header)
            for line in self:
                writer.writerow(line.line)

class BetweenDict():
    '''Creates a dictionary object that contains GENEids that look up a list of tuples for
    the region that their 3'UTR start/stops. 

    FIRST:
    Best to initialize with a bed file. 
    bD = BetweenDict("location/to/bedfile.bed")

    You can feed the dict a location and it will return GENEids that match that location (or range)

    bD.lookup("chr1","+",100000)  will look up a singe position
    bD.lookup("chr1","+",100000,200000) will look up a range and find any overlapping geneIDs. 

    returntable: instead of returning a geneid list, return the whole dataframe of results.
        Default(False)

    '''

    def __init__(self,bedfile):
        import pandas as pd
        self.bedfile=bedfile
        self.myTable = pd.read_table(self.bedfile,names=["chrom","start","stop","geneid","misc","strand"],
                                    usecols=[0,1,2,3,4,5],
                                    index_col=None)
    
    def lookup(self,chrom,strand,start,stop=None,returntable=False):
        tempTable = self.myTable[(self.myTable.chrom==chrom)]
        tempTable = tempTable[tempTable.strand==strand]
        if stop:
            tempTable = tempTable[tempTable.stop>int(start)]
            tempTable = tempTable[tempTable.start<int(stop)]
        else:
            tempTable = tempTable[tempTable.start<int(start)]
            tempTable = tempTable[tempTable.stop>int(start)]
        if returntable:
            return tempTable
        return set(tempTable.geneid)

# OLD VERSION (v1.2)
# class BetweenDict(dict):
#     '''Creates a dictionary object that contains GENEids that look up a list of tuples for
#     the region that their 3'UTR start/stops. 
    
#     FIRST:
#     Best to initialize with a bed file. 
#     bD = BetweenDict("location/to/bedfile.bed")

#     You can feed the dict a location and it will return GENEids that match that location (or range)

#     bD("chr1","+",100000)  will look up a singe position
#     bD("chr1","+",100000,200000) will look up a range and find any overlapping geneIDs. 

#     '''

#     def __init__(self, d={}, caps=False, header=False):
#         '''
#         Need to send a dict containing a {(start,stop,chrom):value}
#         with the tuple being a range where the value is true
#         '''
#         if type(d) == dict:
#             for k, v in d.items():
#                 self[k] = v
#         elif type(d) == str:
#             import csv
#             with open(d, "r") as f:
#                 reader = csv.reader(f, delimiter="\t")
#                 if header: reader.next()
#                 for line in reader:
#                     chrom, start, stop, geneid, zero, strand = line[:6]  # Define line elements
#                     if caps: geneid = geneid.upper()
#                     self[(chrom, strand, int(start), int(stop))] = geneid
#         else:
#             raise TypeError('Requires a dictionary {(start,stop,chr):GeneId} or a bedfile')

#     def __getitem__(self, key, debug=False):
#         '''Requires a tuple (chrom,strand,position)'''
#         s = set()
#         if type(key) == str:  # Called if only a chromosome is requested
#             return dict.__getitem__(self, key)
#         if len(key) == 3:
#             # if a single position is given (chrom,strand,position)
#             chrom, strand, position = key
#             for k, v in self[chrom + strand].items():
#                 if k[0] <= position < k[1]:
#                     for i in v:
#                         s.add(i)
#             return s
#         if len(key) == 4:
#             # if a range of positions (chrom,strand,start,stop)
#             chrom, strand, start, stop = key
#             start = int(start)
#             stop = int(stop)
#             try:
#                 for k, v in self[chrom + strand].items():
#                     if max(0, min(stop, k[1]) - max(start, k[0])) > 0:
#                         for i in v:
#                             s.add(i)
#             except KeyError:
#                 s = [None]
#             return s

#         raise KeyError("Key '%s' is not between any values in the BetweenDict" % key)

#     def __setitem__(self, key, value):
#         from collections import defaultdict
#         if len(key) == 4:
#             chrom, strand, start, stop = key
#             if not chrom + strand in self:
#                 dict.__setitem__(self, chrom + strand, defaultdict(list))
#             if start < stop:
#                 # dict.__setitem__(self[chrom+strand], (start, stop), value)
#                 dict.__getitem__(self, chrom + strand)[(start, stop)].append(value)
#             else:
#                 raise RuntimeError('First element of a BetweenDict key '
#                                    'must be strictly less than the '
#                                    'second element')
#         else:
#             raise ValueError('Key of a BetweenDict must be an iterable '
#                              'with length two')

#     def __contains__(self, key):
#         try:
#             return bool(self[key]) or True
#         except KeyError:
#             return False






# class GeneTranslator():
#     '''DEPRECIATERD
#         Super slow'''

#     def __init__(self,bed,caps=False):
#         import csv
#         from collections import defaultdict
#         self.bed = bed
#         self.genedict = defaultdict(list)
#         with open(self.bed,"r") as f:
#             reader = csv.reader(f,delimiter="\t")
#             for line in reader:
#                 chrom,start,stop,geneid,zero,strand=line[:6] # Define line elements
#                 if caps: geneid = geneid.upper()  # Indexes gene in all caps if caps=True
#                 self.genedict[geneid].append((chrom,int(start),int(stop),strand))

#     def bypos(self,chrom,position,strand):
#         ''' Looks at EVERY geneid to see if your value is
#         between the start and stop values.
#         '''
#         genes=[]
#         n=0
#         for key,value in self.genedict.iteritems():
#             for lookupchrom,start,stop,lookupstrand in value:
#                 if (start<position<stop) and \
#                         (strand==lookupstrand) and \
#                         (lookupchrom==chrom) :
#                     genes.append(key)
#         if len(genes) > 0:
#             return genes
#         else:
#             return None

#     def byrange(self,chrom,start,stop,strand):
#         '''looks up by a range instead of a position'''
#         genes=[]
#         n=0
#         if start==stop:
#             print self.bypos(chrom,start,strand)
#         for key,value in self.genedict.iteritems():
#             for lookupchrom,lookupstart,lookupstop,lookupstrand in value:
#                 if (self.checkoverlap(start,stop,lookupstart,lookupstop)) and \
#                         (strand==lookupstrand) and \
#                         (lookupchrom==chrom) :
#                     genes.append(key)
#         if len(genes) > 0:
#             return genes
#         else:
#             return None

#     def checkoverlap(self,st1,sp1,st2,sp2):
#         if max(0, min(sp1, sp2) - max(st1, st2))>0:
#             return True
#         else:
#             return False


def depths_to_dictionary(depth_file, bedgraph = False):
    '''
    Takes in a depths file. 
    Create with: 
    samtools view -b -F 16 [alignment.bam] | samtools depth -b ~/Documents/Robin/genomes/refseq_3utr.bed - > watson_depths.depth
    Makes a nested dictionary with the first
        level labeled as the chromosome, and the second level 
    Returns dictionary: `{chr:{location:depth,location:depth...}}

    Can use bedgraph=True if working with a Bedtools Coverage file. 
    '''
    print("making depth_dictionary from ", depth_file)
    with open(depth_file, 'r') as depth_in:
        # print depth_in
        depth_reader = csv.reader(depth_in, delimiter="\t")
        depthdict = {}
        chromosome = ""
        for i in depth_reader:
            position,depth = i[1],i[2]
            # i = [chr,position,depth]
            
            if bedgraph:
                start,stop,depth = i[1],i[2],i[3]
                for position in xrange(int(start),int(stop)):
                    if i[0] == chromosome:
                        depthdict[chromosome][str(position)] = float(depth)
                    else:
                        chromosome = i[0]
                        depthdict[chromosome] = {}
                        depthdict[chromosome][str(position)] = float(depth)

            else:
                if i[0] == chromosome:
                    depthdict[chromosome][position] = depth
                else:
                    chromosome = i[0]
                    depthdict[chromosome] = {}
                    depthdict[chromosome][position] = depth
        return depthdict


def max_depths(bed_file, depthDictionary, conditions=[""]):
    '''
    Finds the maximum depths from the regions of a bedfile.
    Requires depths_to_dictionary created depth dictionary.
    Returns the original bedlines with the expression added on. 
    If multiple conditions, Return:
        bedlines    cond1   cond2 ....
    '''
    import csv
    with open(bed_file, 'r') as f:
        results = []
        peak_reader = csv.reader(f, delimiter="\t")
        for n, peak in enumerate(peak_reader):
            # if n == 500: break
            chrom, start, stop, geneid, score, strand = peak
            # For each condition, get the depths
            expression = []
            for n, condition in enumerate(conditions):
                depths = get_depths(chrom, start, stop, strand, conditions[n], depthDictionary)
                if len(depths) > 0:
                    expression.append(max(depths))
                else:
                    expression.append(0)
            results.append([chrom, start, stop, geneid, score] + expression)
    return results


# def find_geneid(chromosome,location,strand=False,annotation_bed="//DarthRNA/Documents/Robin/genomes/refseq_3utr.bed"):
#     '''
#     DEPRECIATED: Super Slow Garbage
#     Take a chromosome and a location and a reference (default is refseq_3utr)
#     Return a list of the genes that intersect that location
#     ex: clip_tools.find_geneid("chr3",133466156)
#         Returns ['Tet2']
#     '''
#     import csv
#     genes=[]
#     with open(annotation_bed,'r') as f:
#         reader = csv.reader(f,delimiter="\t")
#         for line in reader:
#             chrom,start,stop,geneid,zero,linestrand = line
#             if (int(start)<location<int(stop)) and (chrom==chromosome):
#                     if strand:
#                         if strand==linestrand:
#                             pass
#                         else:
#                             continue
#                     if geneid in genes:
#                         pass
#                     else:
#                         genes.append(geneid)
#     return genes

def bedheader(**kargs):
    """Returns a string for writing to a bed file trackline:
    name="My Bed File" description="test.py v1.5 cliptools v.12.0"
    
    Takes arguments: 
    filename=""    # usually send __file__
    version=""     # usually send __version__
    name=""        # easy to understand name of track
    color=[Bolean] # set True if using a color track
    
    """
    filename = ""
    version = ""
    trackname = ""
    color = ""
    cliptoolsinfo = __version__
    if "filename" in kargs:
        filename = kargs["filename"]
    if "version" in kargs:
        version = kargs["version"]
    if "name" in kargs:
        name = kargs["name"]
    if "color" in kargs:
        if kargs["color"] == True:
            color = " itemRgb=\"\"On\"\""
    print(version)
    header = "name=\"%s\" description=\"%s %s cliptools %s\"%s" % (name, filename, version, cliptoolsinfo, color)
    return header


def getfasta(chrom, start, stop, strand="+", genome=defaultgenome):
    '''
    Takes a chromosome,start,stop,strand,genome
    and converts to strand specific sequence.
    '''
    region = chrom + ":" + str(start) + "-" + str(stop)
    x = subprocess.check_output(["samtools", "faidx", genome, region])
    x = "".join(x.split("\n")[1:])
    if strand == "-":
        x = clip_tools.reversecomplement(x)
    return x


def reversecomplement(sequence):
    '''
    Feed a sequence of ACTG to get the reverse complement.
    Currently only returns uppercase, and will convert
    lowercase to uppercase
    '''
    sequence = sequence.upper()
    reference = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse = sequence[::-1]
    rc = ""
    for nt in reverse:
        rc = rc + reference[nt]
    return rc


def parse_TS_miRNA_names(miRNAstring):
    ''' Takes names like 'miR-129-5p/129ab-5p' from Targetscan,
    and converts them to 
    [miR-192-5p,miR192ab-5p].
    Also accepts let(s)
    
    Example: parse_TS_miRNA_names('miR-124/124ab/506')
    returns: ['miR-124', 'miR-124ab', 'miR-506']
    '''
    import re
    f1 = re.compile("(let-|miR-|/)\d+([a-z]+)?-?(3p|5p)?")
    f2 = re.compile("let-|miR-")
    f3 = re.compile("(let-|miR-|/)(\d+([a-z]+)?-?(3p|5p)?)")
    x = f1.finditer(miRNAstring)

    miRNAs = []

    for i in x:
        istring = i.group()
        firstgroup = f2.match(istring)
        if firstgroup:
            header = firstgroup.group()
        number = f3.match(istring).group(2)
        miRNAs.append(header + number)
    return miRNAs

def miRbase_reads(miRNAstring,quiet=False, kind="readcount"):
    '''
    Takes in a MATURE miRNA sequence and returns the counts on miRBase. 
    If you arent sure which one you want, you can be ambiguous
    ex miR-29a instead of miR-29a-3p
    Helps to determine which is the star strand. 
    '''
    import pandas as pd
    mirna_db = pd.read_csv("Data/genomes/mature_with_counts.csv",index_col=0)
    results = mirna_db[mirna_db.index.str.contains(miRNAstring)]
    if results.ReadCount.count()>1 and not(miRNAstring in results.index) :
        if not quiet: 
            print("Found multiple hits. Returning all:")
        return results
    else:
        idx = results.index[0]
        if kind == "readcount":
            return int(results.ReadCount)
        if kind == "sequence":
            return str(results.loc[idx,"sequence"])
        if kind == "ascension":
            return str(results.loc[idx,"ascension"])
        if kind =="seed":
            return clip_tools.getSeed(str(results.loc[idx,"sequence"]))
    raise KeyError


def getSeed(mature_miRNA):
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna
    '''takes in a mature miRNA sequence, returns 8mer

    mRNA    | | N N N N N N N N | | |
    miRNA       8 7 6 5 4 3 2 1 
    8mer    | | N N N N N N N A | | |
    7mer-m8 | | N N N N N N N | | | |
    7mer-A1 | | | N N N N N N A | | |
    6mer    | | | N N N N N N | | | |
    '''

    my_rna = Seq(mature_miRNA, generic_rna) #sequence is RNA
    mRNA = str(my_rna.back_transcribe().reverse_complement()[-8:])
    seed = mRNA[:-1] # 2-7
    #print "found seed " + seed
    mer6 = seed[1:7] # 2-6
    mer8 = seed+"A"
    mer7a = seed[1:7]+"A"
    mer7m8 = str(seed)
    seeds = {"8mer":mer8,"7mer-m8":mer7m8,"7mer1a":mer7a,"6mer":mer6}
    return seeds


if __name__ == "__main__":
    print(__doc__)
    print(dir(clip_tools))
    print("test")

def easylocation(location,to_easy=False):
    ''' Send a location from IGV:
    'chrX:42,521,693-42,522,443'
    and return
    chrom, start,stop (as integers)
    '''
    if to_easy==False:
        chrom,remainder = location.split(":")
        start = remainder.split("-")[0].replace(",","")
        stop = remainder.split("-")[1].replace(",","")
        return chrom,int(start),int(stop)

