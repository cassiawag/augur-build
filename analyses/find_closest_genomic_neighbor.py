import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


def hamming(array1, array2):
    return np.sum(array1 != array2)

def seattle_strains(metadata_file):
    '''
    creates list from metadata of all Seattle strains
    '''
    seattle_strains_list = []
    for mname in metadata_file:
        with open(mname) as mfile:
            meta_file = pd.read_csv(mfile, sep='\t')
            meta_strains_seattle = meta_file[meta_file["region"] == 'seattle']
            seattle_strains_list = meta_strains_seattle.loc[:,'strain'].tolist()
    return seattle_strains_list


def mapping_sequences(files):
    '''
    returns dictionary mapping of all strains across input files
    '''
    mapping = {}
    for fname in files:
        with open(fname) as fasta_file:
            seq_file = SeqIO.parse(fasta_file, "fasta")
            for record in seq_file:
                mapping[record.name] = np.array([])
                sequ = str(record.seq)
                sequ = sequ.replace("A", "0")
                sequ = sequ.replace("T", "1")
                sequ = sequ.replace("G", "2")
                sequ = sequ.replace("C", "3")
                array = np.asarray(list(sequ), dtype = int)
                mapping[record.name] = np.concatenate((mapping[record.name], array), axis=0).astype(int)
    return mapping

def genetic_distance (files, mapping, strain_list):
    '''
    returns dictionary of the lowest genetic distance to the closest sample (from the entire dataset) for each Seattle strain
    '''
    genetic_distance = {}
    for fname in files:
        with open(fname) as fasta_file:
            seq_list = list(SeqIO.parse(fasta_file, "fasta"))
            for strainA in seq_list: #strainA is the strain of intrest and should be in region = Seattle
                distance_hold = None #resets the counters for every strainA
                strain_neighbor= None
                for strainB in seq_list:
                    if (strainA.name in strain_list) and (strainA.name != strainB.name):
                        distance = hamming(mapping[strainA.name], mapping[strainB.name])
                        if (distance_hold is None) or (distance_hold > distance):
                            distance_hold = distance
                            strain_neighbor = strainB.name
                        genetic_distance[strainA.name] = distance_hold
    return genetic_distance

def write_histo(gen_dis_dict, histo = "analyses/genetic_distance_seattle.png", table = "analyses/genetic_distance_h3n2.tsv"):
    if table != None:
        mydict = [gen_dis_dict]
        with open(table, 'w') as fh:
            pd.DataFrame(mydict).to_csv(fh, sep='\t')
    if histo != None:
        fig, ax1 = plt.subplots(figsize = (15,15))
        ax1.hist(gen_dis_dict.values(),bins = len(gen_dis_dict.values()),width = 0.6, align = 'left', color='g')
        plt.xlabel('distance in # of (nucleotides)', size = '20')
        plt.ylabel('count', size = '20')
        plt.title('Genetic Distance to Closest Sample for Seattle H3N2 Flu Strains', size = '25')
        fig.savefig(histo)
        plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="For the Seattle strains, finds the closest genetic distance when compared to all other viruses in the dataset",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', nargs='+', type=str, required=True, help="metadata per lineage based on HA segment")
    parser.add_argument('--alignment', nargs='+', type=str, required=True, help= "aligned sequences")
    parser.add_argument('--output-figure', type = str, help="name of the file to write figure to")
    parser.add_argument('--output-table', type = str, help="name of the file to write TSV to")
    args = parser.parse_args()

    strains_in_seattle = seattle_strains(args.metadata)

    mapping = mapping_sequences(args.alignment)

    lowest_genetic_distance = genetic_distance(args.alignment, mapping, strains_in_seattle)

    write_histo(lowest_genetic_distance, histo = args.output_figure, table = args.output_table)
