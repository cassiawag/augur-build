'''
Assign selection priorities for each strain based on distance to Seattle strains
Priority is the distance from a strain to closest Seattle strain, including Ns
This naturally penalizes incomplete genomes
'''
import argparse
import re
import numpy as np
import pandas as pd
import Bio.SeqIO


'''
Calculate score between strains
This is the total length of overlapping coverage minus Hamming distance
Or alternatively the total length of overlapping matching nucleotides
'''
def sequence_score(seq_array1, seq_array2, coverage_array1, coverage_array2):
    overlap_array = np.bitwise_and(coverage_array1, coverage_array2)
    difference_array = overlap_array * seq_array1 != overlap_array * seq_array2
    return int(np.sum(overlap_array) - np.sum(difference_array))

'''
Return tuple of Seattle and global strain name lists
'''
def collect_strains(alignment, metadata):
    records = [record for record in Bio.SeqIO.parse(alignment, 'fasta')]
    strains = []
    for record in records:
        strains.append(record.name)

    seattle_strains = []
    global_strains = []

    with open(metadata) as metadatafile:
        meta = pd.read_csv(metadatafile, sep='\t')
        seattle_df = meta[meta["region"] == str("Seattle")]
        seattle_strains_meta = seattle_df['strain'].tolist()
        global_df = meta[meta["region"] != str("Seattle")]
        global_strains_meta = global_df['strain'].tolist()
        for strain in strains:
            if strain in seattle_strains_meta:
                seattle_strains.append(strain)
            if strain in global_strains_meta:
                global_strains.append(strain)

    return (seattle_strains, global_strains)

'''
Return dictionary mapping of strain to seq array
'''
def sequence_mapping(alignment, strains):
    mapping = {}
    records = [record for record in Bio.SeqIO.parse(alignment, 'fasta')]

    for strain in strains:
        match = next((record for record in records if record.name == strain), None)
        if match:
            seq = str(match.seq)
            seq = seq.replace("A", "0")
            seq = seq.replace("T", "1")
            seq = seq.replace("G", "2")
            seq = seq.replace("C", "3")
            seq = re.sub(r'[A-Z]', '4', seq)
            array = np.asarray(list(seq), dtype = int)
            mapping[strain] = array

    return mapping

'''
Return dictionary mapping of strain to coverage array
'''
def coverage_mapping(alignment, strains):
    mapping = {}
    records = [record for record in Bio.SeqIO.parse(alignment, 'fasta')]

    for strain in strains:
        match = next((record for record in records if record.name == strain), None)
        if match:
            seq = str(match.seq)
            seq = seq.replace("A", "1")
            seq = seq.replace("T", "1")
            seq = seq.replace("G", "1")
            seq = seq.replace("C", "1")
            seq = re.sub(r'[A-Z]', '0', seq)
            array = np.asarray(list(seq), dtype = int)
            mapping[strain] = array

    return mapping

'''
Collect priorities
'''
def priority_mapping(strains, seattle_strains, sequence, coverage):
    mapping = {}
    counter = 0
    interval = 100
    length = len(strains)
    print("progress")
    for strain in strains:
        if counter % interval == 0:
            print("[", end = '')
            for x in range(int(counter/interval)):
                print("-", end = '')
            for x in range(int(length/interval) - int(counter/interval)):
                print(" ", end = '')
            print("]")
        score = 0
        for sstrain in seattle_strains:
            temp_score = sequence_score(sequence[strain], sequence[sstrain], coverage[strain], coverage[sstrain])
            if temp_score > score:
                score = temp_score
        mapping[strain] = score
        counter += 1
    return mapping

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign selection priorities for each strain based on genome completeness and distance to Seattle strains",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type = str, required = True, help = "genome alignment FASTA")
    parser.add_argument('--metadata', type = str, required = True, help = "metadata TSV")
    parser.add_argument('--output', type = str, required = True, help = "output TSV mapping of strain ")
    args = parser.parse_args()

    # collect Seattle strains and global strains separately and together
    seattle_strains, global_strains = collect_strains(args.alignment, args.metadata)
    strains = list(set(seattle_strains + global_strains))

    print(str(len(global_strains)) + " global strains")
    print(str(len(seattle_strains)) + " Seattle strains")
    print(str(len(strains)) + " total strains")        

    # mapping of strains to sequence
    sequence = sequence_mapping(args.alignment, strains)

    # mapping of strains to coverage
    coverage = coverage_mapping(args.alignment, strains)

    # mapping of strains to priority
    strain_to_priority = priority_mapping(strains, seattle_strains, sequence, coverage)

    with open (args.output, "w") as outfile:
        for strain, priority in strain_to_priority.items():
            outfile.write(strain + '\t' + str(priority) + '\n')
