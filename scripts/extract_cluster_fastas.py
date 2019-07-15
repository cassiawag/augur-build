'''
Creates full-genome fasta files for each clusters.
'''
import os
import json
import argparse
import pandas as pd

'''
Load clusters into usable python dictionary
'''
def clusters(file, min_size):
    diction = {}
    with open(file) as jfile:
        json_data = json.load(jfile)
        for strain, cluster_dict in json_data["nodes"].items():
            cluster_number = cluster_dict['cluster']
            if cluster_number in diction:
                diction[cluster_number][strain] = []
            else:
                diction[cluster_number] = {strain : []}
    dict_scan = dict(diction)
    for cluster_id, cluster_strains in dict_scan.items():
        if len(cluster_strains) < min_size:
            del diction[cluster_id]
    return diction


'''
Adds sequences into python dictionary by cluster
'''
def sequences(clusters, files):
    for fname in files:
        with open(fname) as jfile:
            json_data = json.load(jfile)
            for cluster_id, cluster_strains in clusters.items():
                for strain in cluster_strains.keys():
                    if strain in json_data["nodes"].keys():
                        seq = json_data["nodes"][strain]["sequence"]
                        clusters[cluster_id][strain].append(seq)
    for cluster_id, cluster_strains in clusters.items():
        for strain, sequence in cluster_strains.items():
            clusters[cluster_id][strain] = ''.join(sequence)
    return clusters

'''
Filters clusters to only those that contain at least one strain from selected regions
if region not selected, then it just returns the above clusters
'''
def filter_to_region(clusters, metadata_file = None, region = None):
    if region is None:
        return clusters
    else:
        strains_region_list = []
        cluster_region = []
        final_clusters = {}
        for mfname in metadata_file:
            with open(mfname) as mfile:
                meta_file = pd.read_csv(mfile, sep='\t')
                meta_strains_region = meta_file[meta_file["region"] == str(region)]
                strains_region_list = meta_strains_region.loc[:,'strain'].tolist()

        for cluster, strain_att in clusters.items():
            check = any(strain in strain_att.keys() for strain in strains_region_list)
            if check:
                cluster_region.append(cluster)

        final_clusters = clusters.copy()

        for cluster in clusters.keys():
            if cluster not in cluster_region:
                del final_clusters[cluster]
        return final_clusters

'''
Outputs python dictionary sequences into fasta files
'''
def export_genomes(seq_dict, output_dir):
    for cluster, cluster_strains in seq_dict.items():
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        filename = output_dir + '/cluster' + str(cluster) + '.fasta'
        with open (filename, "w") as fasta:
            for strain, seq in cluster_strains.items():
                fasta.write('>' + strain + '\n' + seq + '\n')
    return fasta

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create full-genome fasta files of clusters, output one FASTA per cluster",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--clusters', type = str, required = True, help = "cluster JSON files")
    parser.add_argument('--nt-muts', nargs = '+', type = str, required = True, help = "list of nt-muts JSON files")
    parser.add_argument('--metadata', nargs = '+', type = str, help = "metadata of strains, one for each segment")
    parser.add_argument('--filter-to-region', type = str, help = 'filters to clusters only containing strains from that region')
    parser.add_argument('--min-size', type = int, default = 2, help = "Minimum number of strains in cluster. Default is 2.")
    parser.add_argument('--output-dir', type = str, required = True, help = "output directory")
    args = parser.parse_args()

    # Create clusters dictionary
    cluster_dict = clusters(args.clusters, args.min_size)

    # Makes sequences dictionary
    cluster_seq_dict = sequences(cluster_dict, args.nt_muts)

    #filters for clusters with region strains
    region_clus = filter_to_region( cluster_seq_dict, metadata_file = args.metadata, region = args.filter_to_region)

    # Outputs fasta files
    export_genomes(region_clus, args.output_dir)
