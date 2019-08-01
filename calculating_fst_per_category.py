import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np

def extracting_categories(metadata_file, label):
    '''
    creates a dictionary of all strains where the keys are the label of interest
    '''
    cat_dict = {}
    for mname in metadata_file:
        with open(mname) as mfile:
            meta_file = pd.read_csv(mfile, sep='\t')
            for category, df_cat in meta_file.groupby(label):
                cat_dict[label] =  df_cat.loc[:,'strain'].tolist()
    return cat_dict

def mapping_sequences(files):
    '''
    returns dictionary mapping of all strains across input files. All ambiguous letter codes are treated equally.
    '''
    mapping = {}
    ambig = ['N', 'R', 'M', 'S', 'Y', 'W', 'K']
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
                for letter in ambig:
                    sequ = sequ.replace(letter, "4")
                array = np.asarray(list(sequ), dtype = int)
                mapping[record.name] = np.concatenate((mapping[record.name], array), axis=0).astype(int)
    return mapping

def hamming(array1, array2):
    '''
    calculates the number of nucleotide differences per site for each pair
    '''
    array1_mask = np.ma.masked_where(array1 > 3, array1)
    array2_mask = np.ma.masked_where(array2 > 3, array2)
    count = min(np.ma.count(array1_mask), np.ma.count(array2_mask))
    return (np.sum(array1_mask != array2_mask))/count

def genetic_distance (cat_dict, mapping):
    '''
    creates a dictionary for within-category genetic distance and all-category genetic distance
    '''
    all_cat_dict = {}
    within_cat_dict = {}
    for catA, strain_list in cat_dict.items():
        for strainA in strain_list:
            for strainB in strain_list:
                if (strainA != strainB) and (strainA in mapping.keys()) and (strainB in mapping.keys()):
                    within_cat_distance = hamming(mapping[strainA], mapping[strainB])
                    within_cat_dict.setdefault(catA, []).append(within_cat_distance)

            for catC, strainC_list in cat_dict.items():
                for strainC in strainC_list:
                    if (strainA != strainC) and (strainA in mapping.keys()) and (strainC in mapping.keys()):
                        '''if you want to make it within vs. outside, can just add and (catA != catC)'''
                        all_cat_distance = hamming(mapping[strainA], mapping[strainC])
                        all_cat_dict.setdefault(catA, {}).setdefault(catC,[]).append(all_cat_distance)
    return within_cat_dict, all_cat_dict

def calculating_fst(within_cat_dict, all_cat_dict):
    '''
    calculates Fst analogue comparing within category distance with all-category distance
    outputs a TSV file with the results
    '''
    ave_dict = {}
    for catA, distancesA in within_cat_dict.items():
        all_cat_distances = []
        cat_average = sum(distancesA)/len(distancesA)
        for catB, distancesB in all_cat_dict[catA].items():
            all_cat_distances.extend(distancesB)
        all_cat_average = sum(all_cat_distances)/len(all_cat_distances)
        ave_dict[catA] = (cat_average, all_cat_average)

    final_dict = {}
    for cat, averages in ave_dict.items():
        fst = ((averages[1] - averages[0])/averages[1])
        final_dict[cat] = fst
    final = [final_dict]

    with open("fst_calculations.tsv", 'w') as fh:
        pd.DataFrame(final).to_csv(fh, sep='\t', index = False)

    return final_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculates an Fst analogue based on input label.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', nargs='+', type=str, required=True, help="metadata for desired pathogen")
    parser.add_argument('--alignment', nargs='+', type=str, required=True, help= "aligned sequences")
    parser.add_argument('--label', type = str, help="name of the file to write figure to")
    args = parser.parse_args()


categories_dic = extracting_categories(args.metadata, args.label)

mapping = mapping_sequences(args.alignment)

within_cat_dict, all_cat_dit = genetic_distance(categories_dic, mapping)

final_fst_dict = calculating_fst(within_cat_dict, all_cat_dit)
