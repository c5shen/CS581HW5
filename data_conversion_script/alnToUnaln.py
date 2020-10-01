import re
import sys

def convertToUnalign(f):
    path = '/'.join(f.split('/')[:-1])
    out = open(path + '/rose.unaln.true.fasta', 'w')

    data = open(f, 'r').read().split('\n')[:-1]
    combined = []
    for line in data:
        if line == '':
            #out.write('\n')
            continue
        if line[0] == '>':
            if len(combined) > 0:
                out.write(''.join(combined) + '\n')
                combined.clear()
            out.write(line + '\n')
        else:
            combined.append(re.sub(r'[-]*', '', line))
            #out.write(re.sub(r'[-]*', '', line) + '\n')
    if len(combined) > 0:
        out.write(''.join(combined) + '\n')
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

    # also parse for 16S.M
    f = pre + '16S.M/R0/cleaned.alignment.fasta'
    data = open(f, 'r').read().split('\n')[:-1]
    out = open(pre + '16S.M/R0/cleaned.unaln.fasta', 'w')
    print("Processing {}".format(f))
    combined = []
    for line in data:
        if line == '':
            #out.write('\n')
            continue
        if line[0] == '>':
            if len(combined) > 0:
                out.write(''.join(combined) + '\n')
                combined.clear()
            out.write(line + '\n')
        else:
            combined.append(re.sub(r'[-]*', '', line))
            #out.write(re.sub(r'[-]*', '', line) + '\n')
    if len(combined) > 0:
        out.write(''.join(combined) + '\n')
    out.close()
if __name__ == "__main__":
    main()
