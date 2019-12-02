'''
This script plots genetic distance from closest viral strain vs. geographic distance from closest viral strain.
The inputs are --metadata, --alignment, and --lat-longs. The output is a PNG.
'''

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopy.distance
from scipy import stats
import statsmodels.api as sm

def strains_mapping(alignment, metadata):
    '''
    Returns dictionary mapping strains to location and sequence.
    Only includes strains with region = Seattle.
    Removes strains with location = '?'.
    '''
    diction = {}
    with open(alignment) as align:
        for index, record in enumerate(SeqIO.parse(align, 'fasta')):
            seq = str(record.seq)
            seq = seq.replace('A', '0')
            seq = seq.replace('T', '1')
            seq = seq.replace('G', '2')
            seq = seq.replace('C', '3')
            seq_array = np.asarray(list(seq), dtype = int)
            diction[record.id] = {'seq' : seq_array, 'location' : ''}
    with open(metadata) as mfile:
        meta_file = pd.read_csv(mfile, sep = '\t')
        for row in meta_file.itertuples():
            if row.region == 'seattle':
                if row.strain in diction:
                    diction[row.strain]['location'] = row.location
    mapping = {strains: attributes for strains, attributes in diction.items()
               if attributes['location'] != '' and attributes['location'] != '?'
               and attributes['location'].startswith('530')}
    #Some of the strains in metadata_h3n2_ha.tsv with region = seattle have
    #locations that are in other states.
    return mapping

def hamming(array1, array2):
    return np.sum(array1 != array2)

def genetic_distance(mapping):
    '''
    Returns dictionary containing strain, closest strain, and genetic distance (computed as hamming distance) to closest strain.
    '''
    gen_distance_dict = {}
    for strainA in mapping:
        gen_distance_low = None
        for strainB in mapping:
            if strainA != strainB:
                gen_distance = hamming(mapping[strainA]['seq'], mapping[strainB]['seq'])
                if gen_distance_low is None or gen_distance < gen_distance_low:
                    gen_distance_low = gen_distance
                    strain_name = strainB
        gen_distance_dict[strainA] = {'gen_distance' : gen_distance_low, 'closest_strain' : strain_name}
    return gen_distance_dict

def location(lat_longs, mapping):
    '''
    Returns dictionary mapping strain to latitude and longitude.
    '''
    location_dict = {}
    with open(lat_longs) as llfile:
        coordinates = pd.read_csv(llfile, sep = '\t')
        for strain in mapping:
            for row in coordinates.index:
                if coordinates.iloc[row, 1] == mapping[strain]['location']:
                    location_dict[strain] = {'latitude' : coordinates.iloc[row, 2], 'longitude' : coordinates.iloc[row, 3]}
    return location_dict

def geographic_distance(gen_distance_dict, location_dict):
    '''
    Returns dictionary mapping strain to geographic distance and genetic distance from closest strain.
    '''
    distance_dict = {}
    for strain in gen_distance_dict:
        closest_strain = gen_distance_dict[strain]['closest_strain']
        strain_location = (location_dict[strain]['latitude'], location_dict[strain]['longitude'])
        closest_strain_location = (location_dict[closest_strain]['latitude'], location_dict[closest_strain]['longitude'])
        geo_distance = geopy.distance.distance(strain_location, closest_strain_location).km
        distance_dict[strain] = {'gen_distance' : gen_distance_dict[strain]['gen_distance'], 'geo_distance' : geo_distance}
    return distance_dict

def scatter_plot(distance_dict, table, figure):
    '''
    Saves a table and scatterplot of genetic distance and geographic distance from each strain to its closest strain.
    '''
    distance_df = pd.DataFrame.from_dict(distance_dict, orient = 'index')

    with open(table, 'w') as tsv:
        distance_df.to_csv(tsv, sep = '\t', index_label = 'strain')

    corr, pvalue = stats.spearmanr(distance_df.gen_distance, distance_df.geo_distance)
    lowess = sm.nonparametric.lowess(distance_df.geo_distance, distance_df.gen_distance, frac=1./3, it=2)
    print(lowess)
    transposed_lowess = lowess.transpose()
    lx, ly = np.vsplit(transposed_lowess, 2)
    plt.figure(figsize=(10,8), facecolor = 'white')
    ax = plt.axes()
    ax.axis(xmin = -1, xmax = np.percentile(distance_df.gen_distance, 99))
    plt.title('Genetic distance vs. geographic distance to closest strain', fontsize=18)
    plt.ylabel('Geographic distance (km)', fontsize=12)
    plt.xlabel('Genetic distance (bp)', fontsize=12)
    ax.text(32.6, 165, "ρ: {0:.{1}f} \nP: {2:.{3}f}".format(corr, 4, pvalue, 6), fontsize=11)
    ax.scatter(distance_df.gen_distance, distance_df.geo_distance, color='#4C90C0', label = '_nolegend_')
    ax.plot(lx[0], ly[0], c='black', linewidth=1, label = 'LOESS')
    ax.legend(loc=1, frameon = False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(figure)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot genetic distance vs. geographic distance from closest strain",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, default='data/metadata_h3n2_ha.tsv', help="name of metadata file")
    parser.add_argument('--alignment', type=str, default='results/aggregated/aligned_h3n2_genome_1y.fasta', help="name of alignment file")
    parser.add_argument('--lat-longs', type=str, default='config/lat_longs.tsv', help="name of lat-longs file")
    parser.add_argument('--output-table', type=str, default='analyses/tables/h3n2_gen_geo_distance.tsv', help="name of output table")
    parser.add_argument('--output-figure', type=str, default='analyses/figures/h3n2_gen_geo_distance.png', help="name of output figure")
    args = parser.parse_args()

    #Make dictionary mapping strain to location and sequence
    strains = strains_mapping(args.alignment, args.metadata)

    #Compute genetic distance for each strain
    gen_distance = genetic_distance(strains)

    #Make dictionary mapping strain to latitude and longitude
    locations = location(args.lat_longs, strains)

    #Compute geographic distance
    geo_distance = geographic_distance(gen_distance, locations)

    #Plot
    scatter_plot(geo_distance, args.output_table, args.output_figure)
