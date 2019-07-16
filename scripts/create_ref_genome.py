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
    model_source_feature = flu_ref.features[0]
    full_genome_feature = SeqFeature(FeatureLocation(0, len(flu_ref), strand = 1), type = "source", qualifiers = model_source_feature.qualifiers)
    #Removes qualifiers that don't apply to the entire sequence.
    if "segment" in full_genome_feature.qualifiers.keys():
        del full_genome_feature.qualifiers["segment"]
    if "lab_host" in full_genome_feature.qualifiers.keys():
        del full_genome_feature.qualifiers["lab_host"]
    flu_ref.features.insert(0, full_genome_feature)

    #Removes features of type source that don't encompass whole genome.
    source_features = []
    for feature in flu_ref.features:
        if feature.type == "source" and flu_ref.features.index(feature) != 0:
            source_features.append(flu_ref.features.index(feature))
    source_features.reverse()
    for index in source_features:
        flu_ref.features.pop(index)

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
