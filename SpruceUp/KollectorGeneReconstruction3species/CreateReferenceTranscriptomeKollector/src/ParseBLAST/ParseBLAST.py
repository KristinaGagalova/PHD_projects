#!/usr/bin/env python

#######################################
######Kristina Gagalova################
#######################################
###########09 Feb 2017#################
#######################################

'''
Description: parse the output from tabular BLAST output (-outfmt "6 std qcov qcovhsp")
example of the input - no header:
k60_k38_30027   k72_k38_3058026 84.21   323     46      4       270     588     347     668     1e-82   309     20      20
k60_k38_166418  k62_k40_9638430 84.34   1245    178     11      1       1231    710     1951    0.0     1203    76      76
k60_k38_208364  k62_k72_4644991 79.94   952     174     11      735     1674    1423    2369    0.0     684     56      56
k60_k38_208364  k60_k56_1056760 100.00  55      0       0       1620    1674    1       55      3e-20   102     3       3
k60_k38_234344  k72_k58_3283932 76.36   846     175     17      470     1301    611     1445    3e-119  431     44      44

Example of the usage:
ParseBLAST.py -f my_file.tsv -i 25 -c 80
'''

# import libraries
import csv
import sys
import getopt


# parse command line
def usage():
    print (
        'usage: ParseBLAST.py -f <file_name> -o <out_dir> -i <identity_th%> -c <qcov_th%> \n'
        'Insert file to be processed, output directory, identity threshold (%) and query coverage theshold (5) \n'
        'Input file should be tab delimited, blast output -outfmt "6 std qcov qcovhsp"')
    sys.exit(2)


def main(argv):
    in_file = ''
    out_dir = ''
    identity = ''
    coverage = ''

    try:
        opts, args = getopt.getopt(argv, "hf:o:i:c:", ["input_file=", "output_dir=", "identity%=", "coverage%="])

    except getopt.GetoptError:
        usage()

    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ("-f", "--input_file"):
            in_file = arg
        elif opt in ("-o", "--output_dir"):
            out_dir = arg
        elif opt in ("-i", "--identity%"):
            identity = arg
        elif opt in ("-c", "--coverage%"):
            coverage = arg

    if in_file == '':
        usage()
    if (int(identity) > 100) or (int(coverage) > 100):
        usage()
    if not in_file.endswith('.tsv'):
        usage()
    return [in_file, out_dir, identity, coverage]


def parsed_blast(in_file, ident_th, cov_th):
    '''
    given the blast output table, create a dictionary of the parsed values
    '''
    dict_vals = dict()
    with open(in_file, 'rb') as f:
        for line in f:
            qseqid, sseqid, pident, evalue, qcov = map(lambda x: line.rstrip().split('\t')[x], [0, 1, 2, 10, 12])
            if float(pident) > ident_th and float(qcov) > cov_th:
                if qseqid not in dict_vals:
                    dict_vals[qseqid] = [[sseqid, pident, evalue, qcov]]
                else:
                    dict_vals[qseqid].append([sseqid, pident, evalue, qcov])
        return dict_vals


def create_out(in_file, out_dir, ident_th, cov_th):
    '''
    given the dictionaries with the parsed data , create the output data frame
    '''
    parsed_vals = parsed_blast(in_file, ident_th, cov_th)
    out_file = in_file.replace(".tsv", "").split("/")[-1]  # create out file name
    w_parsed = csv.writer(open(out_dir + "/parsed_" + out_file + str(int(ident_th)) + "_" + str(int(cov_th)) + ".out", "w"), delimiter='\t')
    w_parsed.writerow(['qseqid', 'sseqid', 'pident', 'evalue', 'qcov'])  # create header
    for qseq in parsed_vals:
        reord_list = list()
        for value in [list(t) for t in zip(*parsed_vals.get(qseq))]:  # reorder values
            reord_list.append(",".join(value))
        w_parsed.writerow([qseq] + reord_list)


if __name__ == "__main__":
    in_file, out_dir, identity, coverage = main(sys.argv[1:])
    create_out(in_file, out_dir, float(identity), float(coverage))
