from Bio import SeqIO
import khmer
import matplotlib.pylab as plt
import pandas as pd

kmers_file = "/projects/btl_scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/KmerCounting/DSK_countKmers/res-k75/dsk_75.txt"
fasta_file = "/projects/btl_scratch/kgagalova/PipelineCNV/ProofOfConcept/CheckCopies/HumanReference/CheckDifferentCopies/KmerCounting/ENSG00000005889.fa"

#use odd numbers so that the window is symmetric
#KMER_SIZE needs to be the same as the kmer counts
KMER_SIZE = 75

def load_kmers(textfile):
    #load the kmers counts from an external txt file
    with open(textfile) as f:
        kmers_counts = dict()
        for line in f:
            cols = line.rstrip().split(' ')
            kmers_counts[cols[0]] = int(cols[1])
        return kmers_counts

def get_myfastaSequence(fasta_file):
    #get the sequence of the first file, more info here: http://biopython.org/wiki/SeqIO
    record_dict = list(SeqIO.parse(open(fasta_file), 'fasta'))
    return record_dict[0].seq

def assign_kmers(fasta_file, kmer_file, kmer):
    fasta = get_myfastaSequence(fasta_file)
    kmer_counts = load_kmers(kmer_file)
    half_kmer = int((kmer -1) / 2)
    kmer_pos = dict()
    for window in range(half_kmer - 1, len(fasta) - half_kmer):
        fasta_kmer= str(fasta[window - half_kmer - 1 : window + half_kmer])
        if fasta_kmer in kmer_counts.keys():
            kmer_pos[window] = kmer_counts[fasta_kmer]
        else:
           kmer_pos[window] = 0
    return kmer_pos

#def write_to_df(dictionary):
#    df = [ pd.DataFrame.from_dict(dictionary) ]
#   return df

#print(load_kmers(kmers_file))

#plot the coverage values
coord_counts = assign_kmers(fasta_file,kmers_file,KMER_SIZE)
#print(test.items())
df = pd.Series(coord_counts)
df.to_csv("OutputCounts.bedgraps", sep='\t',header=False)
#DataFrame.to_csv(df)
#write_to_df(test)
lists = sorted(df.items())
x, y = zip(*lists)
plt.plot(x, y)
plt.axvspan(552, 663, facecolor='g', alpha=0.5)
plt.axvspan(23543, 23628, facecolor='g', alpha=0.5)
plt.axvspan(24435, 24551, facecolor='g', alpha=0.5)
plt.axvspan(30011, 30598, facecolor='g', alpha=0.5)
plt.axvspan(58154, 58303, facecolor='g', alpha=0.5)
plt.axvspan(58540, 58683, facecolor='g', alpha=0.5)
plt.axvspan(59046, 59198, facecolor='g', alpha=0.5)
plt.axvspan(59728, 59868, facecolor='g', alpha=0.5)
plt.axvspan(61021, 67083, facecolor='g', alpha=0.5)
plt.show()