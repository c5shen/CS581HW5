import re
import sys

def convertToUnalign(f):
    path = '/'.join(f.split('/')[:-1])
    out = open(path + '/rose.unaln.true.fasta', 'w')

    data = open(f, 'r').read().split('\n')[:-1]
    for line in data:
        if line == '':
            out.write('\n')
            continue
        if line[0] == '>':
            out.write(line + '\n')
        else:
            out.write(re.sub(r'[-]*', '', line) + '\n')
    out.close()

# To conversion true alignment files back to unaligned fasta files
def main():
    # data are from 1000M1 and 1000M4
    pre = '../../'
    post = 'rose.aln.true.fasta'
    datasets = ['1000M1', '1000M4']
    
    # choose the first 10 replicates for the task (for both datasets)
    for d in datasets:
        for i in range(10):
            f = pre + d + '/' + d + '/R' + str(i) + '/' + post
            print("Processing {}".format(f))
            convertToUnalign(f)

if __name__ == "__main__":
    main()
