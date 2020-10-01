from Bio import SeqIO
import subprocess
import os, shutil, sys
import pexpect
from dendropy.calculate import treecompare
from dendropy.datamodel.treemodel import Tree
import dendropy
tns = dendropy.TaxonNamespace()
binary = '../../FastTree/FastTreeMP'

# compare two trees and return the FP and FN (t2 is reference)
def compareTrees(t1, t2):
    return treecompare.false_positives_and_negatives(t1, t2)

# evaluate results for true_aln, est_aln and pasta_tree on gt_tree
def evaluation(gt_tree, true_aln_tree, est_aln_tree, pasta_tree):
    # 2.0 Ground truth tree encoding
    gt_tree_in = Tree.get_from_path(gt_tree, "newick", taxon_namespace=tns)
    gt_tree_in.encode_bipartitions()
    
    # 2.1 True Alignment inferred tree
    if (not os.path.isfile(true_aln_tree) or
            os.stat(true_aln_tree).st_size == 0):
        true_aln_fpfn = ['-', '-']
    else:
        tree2 = Tree.get_from_path(true_aln_tree, "newick", taxon_namespace=tns)
        tree2.encode_bipartitions()
        true_aln_fpfn = compareTrees(tree2, gt_tree_in)

    # 2.2 Estimated Alignment inferred tree
    if (not os.path.isfile(est_aln_tree) or
            os.stat(est_aln_tree).st_size == 0):
        est_aln_fpfn = ['-', '-']
    else:
        tree2 = Tree.get_from_path(est_aln_tree, "newick", taxon_namespace=tns)
        tree2.encode_bipartitions()
        est_aln_fpfn = compareTrees(tree2, gt_tree_in)

    # 2.3 Maybe put in PASTA generated tree too
    tree3 = Tree.get_from_path(pasta_tree,
                "newick", taxon_namespace=tns)
    tree3.encode_bipartitions()
    pasta_tre_fpfn = compareTrees(tree3, gt_tree_in)

    return true_aln_fpfn, est_aln_fpfn, pasta_tre_fpfn

def main():
    result_file = open('result_pasta.txt', 'w')
    result_file.write('Dataset,Repetition,True Aln FP,True Aln FN,Est Aln FP,Est Aln FN,PASTA FN,PASTA FP\n')

    #16S.M
    #true_aln = '../../16S.M/R0/cleaned.alignment.phylip'
    #est_aln = '16S.M/pastajob.marker001.cleaned.unaln.phylip'
    #true_aln_tree = '16S.M/true_aln_tree.nwk'
    #est_aln_tree = '16S.M/est_aln_tree.nwk'
    #gt_tree = '../../16S.M/16S.M.reference.nwk'
    #pasta_tree = '16S.M/pastajob.tre'
    #print("running FastTree on 16S.M on R0, method GTR")
    #command = binary + ' -nt -gamma -gtr '
    #os.system(command + true_aln + ' > ' + true_aln_tree)
    #os.system(command + est_aln + ' > ' + est_aln_tree)
    #true_aln_fpfn, est_aln_fpfn, pasta_tre_fpfn = evaluation(
    #        gt_tree, true_aln_tree, est_aln_tree, pasta_tree)
    #result_file.write('16S.M'+',R0,'+
    #                ','.join([str(x) for x in true_aln_fpfn])+
    #                ','+','.join([str(x) for x in est_aln_fpfn])+
    #                ','+','.join([str(x) for x in pasta_tre_fpfn])+
    #                '\n')

    #1000M1 and 1000M4
    targets = ['1000M1', '1000M4']
    for target in targets:
        for i in range(10):
            # input file info
            true_aln = '../../'+target+'/'+target+'/R'+str(i)+'/rose.aln.true.phylip'
            est_aln = target+'/R'+str(i)+'/pastajob.marker001.rose.unaln.true.phylip'

            # output file info
            outname = target+'/R'+str(i)
            true_aln_tree = outname + '/true_aln_tree.nwk'
            est_aln_tree = outname + '/est_aln_tree.nwk'
            
            # ground truth tree info
            gt_tree = '../../'+target+'/'+target+'/R'+str(i)+'/rose.tt'

            method = 'GTR'

            print("running FastTree on {} dataset on R{}, method {}".format(
                target, str(i), method))

            m = '-gtr'

            # 1.1 tree build based on true MSA
            command = binary + ' -nt -gamma {} '.format(m)
            os.system(command + true_aln + ' > ' + true_aln_tree)

            # 1.2 tree build based on estimated MSA
            os.system(command + est_aln + ' > ' + est_aln_tree)

            # 2. evaluate tree correctness
            pasta_tree = outname + '/pastajob.tre'
            true_aln_fpfn, est_aln_fpfn, pasta_tre_fpfn = evaluation(
                    gt_tree, true_aln_tree, est_aln_tree, pasta_tree)

            # 3. Write to result file
            # columns are: dataset, repetition, true alignment FP, true alignment FN,
            # est alignment FP, est alignment FN, pasta FP, pasta FN
            result_file.write(target+',R'+str(i)+','+
                    ','.join([str(x) for x in true_aln_fpfn])+
                    ','+','.join([str(x) for x in est_aln_fpfn])+
                    ','+','.join([str(x) for x in pasta_tre_fpfn])+
                    '\n')

    result_file.close()

if __name__ == "__main__":
    main()
