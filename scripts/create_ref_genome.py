'''
Concatenates segment reference genbank files into a genome reference genbank file.
'''
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

'''
Creates full-genome Bio Seq object
'''
def concat(files, reference):
    for fname in files:
        seq = SeqIO.read(fname, 'gb')
        if 'ha' in fname:
            flu_ref = seq
        else:
            flu_ref += seq

    #Genbank forces names to be 16 characters or less.
    if 'California' in reference:
        flu_ref.name = reference.replace('California', 'Cali')
    else:
        flu_ref.name = reference

    flu_ref.description = reference
    flu_ref.id = reference

    #Adds source feature to Genbank file for the full-length of genome.
    source_feature = flu_ref.features[0]
    full_genome_feature = SeqFeature(FeatureLocation(0, len(flu_ref), strand = 1), type = "source", qualifiers = source_feature.qualifiers)
    flu_ref.features.insert(0, full_genome_feature)

    return flu_ref

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create full-genome reference genbank file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--references', nargs = '+', type = str, required = True, help = "List of reference gb files")
    parser.add_argument('--output', type = str, required = True, help = "output location")
    parser.add_argument('--ref-strain', type = str, required = True, help = "Name of reference strain")
    args = parser.parse_args()

# Writes out reference genbank file for full influenza genome
SeqIO.write(concat(args.references, args.ref_strain), args.output, 'gb')
