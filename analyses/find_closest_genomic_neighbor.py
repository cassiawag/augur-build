import os
import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt

def hamming(array1, array2):
    return np.sum(array1 != array2)

def seattle_strains(metadata_file):
    '''
    creates list from metadata of all Seattle strains
    '''
    seattle_strains_list = []
    with open(metadata_file) as mfile:
        meta_file = pd.read_csv(mfile, sep='\t')
        meta_strains_seattle = meta_file[meta_file["region"] == 'seattle']
        seattle_strains_list = meta_strains_seattle.loc[:,'strain'].tolist()
    return seattle_strains_list


def mapping_sequences(alignment_file):
    '''
    returns dictionary mapping of all strains across input files
    '''
    mapping = {}
    with open(alignment_file) as fasta_file:
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

def genetic_distance(alignment_file, mapping, strain_list):
    '''
    returns dictionary of the lowest genetic distance to the closest sample (from the entire dataset) for each Seattle strain
    '''

    counter = 0
    interval = 100
    length = len(mapping)
    print("progress")

    genetic_distance = {}
    with open(alignment_file) as fasta_file:
        seq_list = list(SeqIO.parse(fasta_file, "fasta"))
        for strainA in seq_list: # strainA is the strain of interest and should be in region = Seattle

            if counter % interval == 0:
                print("[", end = '')
                for x in range(int(counter/interval)):
                    print("-", end = '')
                for x in range(int(length/interval) - int(counter/interval)):
                    print(" ", end = '')
                print("]")

            if strainA.name in strain_list:
                distance_hold = None # resets the counters for every strainA
                for strainB in seq_list:
                    if strainA.name != strainB.name:
                        distance = hamming(mapping[strainA.name], mapping[strainB.name])
                        if (distance_hold is None) or (distance_hold > distance):
                            distance_hold = distance
                        genetic_distance[strainA.name] = distance_hold

            counter += 1


    return genetic_distance

def write_histo(gen_dis_dict, histo = "analyses/genetic_distance_h3n2_seattle.png", table = "analyses/genetic_distance_h3n2.tsv"):
    if table != None:
        table_dir = os.path.dirname(table)
        if not os.path.exists(table_dir):
            os.makedirs(table_dir)
        mydict = [gen_dis_dict]
        with open(table, 'w') as fh:
            pd.DataFrame(mydict).to_csv(fh, sep='\t')
    if histo != None:
        histo_dir = os.path.dirname(histo)
        if not os.path.exists(histo_dir):
            os.makedirs(histo_dir)
        fig, ax1 = plt.subplots(figsize = (15, 10))
        ax1.hist(gen_dis_dict.values(), bins = len(gen_dis_dict.values()), width = 0.6, align = 'left', color="#4C90C0")
        ax1.set_xlim(0, 25)
        plt.xlabel('distance (# of nucleotides)', size = '32')
        plt.ylabel('count', size = '32')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.tick_params(axis='both', which='major', labelsize=28)
        ax1.tick_params(axis='both', which='minor', labelsize=28)
        plt.tight_layout()
        fig.savefig(histo)
        plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="For the Seattle strains, finds the closest genetic distance when compared to all other viruses in the dataset",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, default="data/metadata_h3n2_ha.tsv", help="metadata per lineage based on HA segment")
    parser.add_argument('--alignment', type=str, default="results/aggregated/aligned_h3n2_genome_1y.fasta", help= "aligned sequences")
    parser.add_argument('--output-table', type=str, default="analyses/tables/closest_genomic_neighbor_h3n2.tsv", help="name of the file to write TSV to")
    parser.add_argument('--output-figure', type=str, default="analyses/figures/closest_genomic_neighbor_h3n2.png", help="name of the file to write figure to")
    args = parser.parse_args()

    strains_in_seattle = seattle_strains(args.metadata)

    mapping = mapping_sequences(args.alignment)

    lowest_genetic_distance = genetic_distance(args.alignment, mapping, strains_in_seattle)

    write_histo(lowest_genetic_distance, histo = args.output_figure, table = args.output_table)
