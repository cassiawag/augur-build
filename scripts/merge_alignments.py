'''
Creates full-genome alignment FASTA file by merging aligned segment FASTAs
'''
import argparse
import Bio.SeqIO

'''
Return comprehensive list of strain names, representing union across alignments
'''
def collect_strains(records_mapping):
    strains = []
    for alignment, records in records_mapping.items():
        for record in records:
            strains.append(record.name)
    return list(set(strains))

'''
Return list of alignment lengths
'''
def collect_lengths(records_mapping):
    mapping = {}
    for alignment, records in records_mapping.items():
        mapping[alignment] = len(records[0].seq)
    return mapping

'''
Construct full-genome alignment
'''
def write_alignment(strains, records_mapping, length_mapping, output):
    with open (output, "w") as outfile:
        for strain in strains:
            outfile.write('>' + strain + '\n')
            for alignment, records in records_mapping.items():
                match = next((record for record in records if record.name == strain), None)
                if match:
                    outfile.write(str(match.seq))
                else:
                    outfile.write('N' * length_mapping[alignment])
            outfile.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates full-genome alignment FASTA file by merging aligned segment FASTAs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignments', nargs = '+', type = str, required = True, help = "list of aligned segment FASTAs")
    parser.add_argument('--output', type = str, required = True, help = "output file for full-genome alignment FASTA")
    args = parser.parse_args()

    records_mapping = {}
    for alignment in args.alignments:
        records = [record for record in Bio.SeqIO.parse(alignment, 'fasta')]
        records_mapping[alignment] = records

    # collect comprehensive list of strain names
    strains = collect_strains(records_mapping)

    # collect segment-specific length mapping
    length_mapping = collect_lengths(records_mapping)

    # construct full-genome alignment
    write_alignment(strains, records_mapping, length_mapping, args.output)
