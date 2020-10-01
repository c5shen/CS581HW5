from Bio import SeqIO
import sys

class Alignment(object):
    def __init__(self, seq):
        self.seq = seq

def writePhylip(alignment, filePath, taxa = None):
    maxChars = 0
    lines = []
    for tag in alignment:
        if taxa is None or tag in taxa:
            lines.append("{} {}\n".format(tag, alignment[tag].seq))
            maxChars = max(maxChars, len(alignment[tag].seq))
    
    with open(filePath, 'w') as textFile:
        textFile.write("{} {}\n".format(len(lines), maxChars))
        for line in lines:
            textFile.write(line)

# original way of converting it
#input = sys.argv[1]
#output = '.'.join(input.split('.')[:-1]) + ".phylip"
#records = SeqIO.parse(input, "fasta")
#count = SeqIO.write(records, output, "phylip")
#print("Converted %i records" % count)


# vlad's way of converting fasta to phylip
input = sys.argv[1]
output = '.'.join(input.split('.')[:-1]) + ".phylip"
data = open(input, 'r').read().split('\n')[:-1]
taxa = {}
cur_taxon = ''
cur_seq = ''
for i in range(len(data)):
    if data[i] == '':
        taxa[cur_taxon] = Alignment(cur_seq)
        cur_taxon = ''
        cur_seq = ''
        continue

    # check for header
    if data[i][0] == '>':
        if cur_taxon != '':
            taxa[cur_taxon] = Alignment(cur_seq)
        cur_taxon = data[i][1:]
        cur_seq = ''
        continue
    cur_seq += data[i]
if cur_taxon != '':
    taxa[cur_taxon] = Alignment(cur_seq)

#print(taxa)
writePhylip(taxa, output)
print("Converted {} records".format(len(taxa)))
