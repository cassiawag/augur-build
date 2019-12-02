'''
This script finds number of influenza introductions into Seattle. It parses JSONs produced for auspice,
determines the number of distinct introductions into Seattle from other regions (as defined by a change from trait state of an internal node from a non-Seattle region to Seattle),
and determines the number of children attributed to each introduction (e.g. how much transmission occurred post introduction)
and also the frequency of introductions coming from different source regions.

Inputs are: --tree, --output-table, --lineage, --output-size-figure, --colors, and --output-region-figure.

Outputs are:
(1) TSV file with the number of descendants from an introduction and the source of that introduction.
(2) PNG showing introductions by size.
(3) PNG showing introductions by region.
'''

import argparse
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import baltic as bt

def assign_region_trait_to_stub_nodes(baltic_tree):
    """ The tree JSONs made for SFS are in fact collated subtrees.
        Each subtree starts with a 'stub'. These nodes don't have any region trait information.
        This throws KeyErrors when traversing the tree looking for geographic jumps.
        So to avoid this issue, I'm assigning a region value of 'undefined' to the stubs.

        Note that this function takes a baltic tree object specifically."""
    new_tree = baltic_tree
    for k in new_tree.Objects:
        if isinstance(k,bt.node):
            try:
                k.traits['region']
            except KeyError:
                k.traits['region'] = 'undefined'
    return new_tree

def extract_seattle_subtrees(baltic_tree):
    """
    Function takes in a baltic tree object and returns a list of baltic tree objects.
    Each returned tree object represents a point in the tree in which the tree moved from
    a non-seattle deme into seattle. Note that this will return sub-sub trees, e.g.
    where the clade moves into seattle, leaves seattle, and then strains within that clade
    return to seattle (this scenario would make 2 separate subtrees, even though they are nested).

    If that behavior is not desired, the returned list must be pruned out of any trees that have an
    intersection of internal nodes.
    """
    subtrees_list = []
    for k in sorted(baltic_tree.Objects, key = lambda x:x.height):
        if k.traits['region'] == 'ancestor':
            continue
        else:
            parent = k.parent
            k_region = k.traits['region']
            parent_region = k.parent.traits['region']

            if k_region != parent_region and k_region == 'seattle':
                seattle_subtree = baltic_tree.subtree(k)
                subtrees_list.append(seattle_subtree)
    return subtrees_list

def get_basal_subtrees(subtrees_list):
    """
    Given a list of subtrees identified from a baltic_tree, this function removes, sub-subtrees and returns a list of basal subtrees.
    In other words, if the input list contained trees that were subtrees of other trees -- termed basal trees -- in the input list,
    this function only returns the basal trees.
    """
    subtree_leaves = []
    for tree in subtrees_list:
        leaf_list = []
        for k in tree.Objects:
            if isinstance(k,bt.leaf):
                leaf_list.append(k.name)
        subtree_leaves.append(leaf_list)

    non_basal_subtrees = []
    for indexA, leaf_list_A in enumerate(subtree_leaves):
        for indexB, leaf_list_B in enumerate(subtree_leaves):
            if indexA != indexB:
                for leaf in leaf_list_A:
                    if leaf in leaf_list_B:
                        if len(leaf_list_A) < len(leaf_list_B):
                            non_basal_subtrees.append(subtrees_list[indexA])
                        elif len(leaf_list_A) > len(leaf_list_B):
                            non_basal_subtrees.append(subtrees_list[indexB])
                        break

    basal_subtrees = [tree for tree in subtrees_list if tree not in non_basal_subtrees]
    return basal_subtrees

def count_number_of_seattle_leaves(baltic_tree):
    """
    Given a baltic tree or subtree, this function counts how many of the leaves in the tree were sampled from Seattle.
    """
    n_seattle_leaves = 0
    for k in sorted(baltic_tree.Objects, key = lambda x:x.height):
        if isinstance(k,bt.leaf) and k.traits['region'] == 'seattle':
            n_seattle_leaves +=1
    return n_seattle_leaves


def find_region_of_introduction(baltic_tree):
    """
    This function returns the 'region' trait on the root of a subtree
    """
    sorted_tree = sorted(baltic_tree.Objects, key = lambda x:x.height)
    region = sorted_tree[0].parent.traits['region'] # root fill have smallest height, and therefore be 0th element
    return region

def find_clade_of_introduction(baltic_tree):
    """
    This function returns the 'clade_membership' trait of members of a subtree
    """
    clade = "unknown"
    for node in baltic_tree.Objects:
        if 'clade_membership' in node.traits:
            clade = node.traits['clade_membership']
    return clade

def find_date_of_introduction(baltic_tree):
    """
    This function returns the 'num_date' trait on the root of a subtree
    """
    sorted_tree = sorted(baltic_tree.Objects, key = lambda x:x.height)
    date = sorted_tree[0].traits['num_date'] # root fill have smallest height, and therefore be 0th element
    return date

def create_dataframe(basal_subtrees, table):
    """
    Creates dataframe from basal subtrees for origin of introduction and the number of leaves in Seattle.
    Writes out dataframe as TSV to table.
    """
    data = []
    for tree in basal_subtrees:
        data.append([count_number_of_seattle_leaves(tree),
                     find_region_of_introduction(tree),
                     find_clade_of_introduction(tree),
                     find_date_of_introduction(tree)])
    seattle_introductions_df = pd.DataFrame(data)
    seattle_introductions_df.columns = ["size", "region", "clade", "root_date"]
    with open(table, 'w') as tsv:
        seattle_introductions_df.to_csv(tsv, sep='\t', index=False)
    return seattle_introductions_df

def plot_intro_by_size(dataframe, lineage, output):
    """
    Plots histogram of introductions by size.
    """
    mpl.rcParams['font.weight']=110
    mpl.rcParams['axes.labelweight']=110
    mpl.rcParams['font.size']=14

    # make histogram of introductions by size
    fig, ax = plt.subplots(figsize=(10,8), facecolor='white')
    plt.hist(dataframe['size'], bins = 50, color = '#4292c6')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_aspect(0.8)
    ax.set_xlabel("Number of sequenced cases in Seattle")
    ax.set_ylabel("Number of introductions")
    ax.set_title(lineage)
    # ax.set_xlim(1,190)
    # ax.set_ylim(0,100)
    plt.tight_layout()

    plt.savefig(output, dpi=250)

def plot_intro_by_region(dataframe, lineage, colors, output):
    """
    Plots histogram of introductions by region.
    """
    with open(colors) as cfile:
        color_df= pd.read_csv(cfile, sep = '\t', header = None, )
    color_df.columns = [ 'classification', 'region', 'color']
    color_df = color_df[color_df.classification == 'region']
    color_df['label'] = ['China', 'Southeast Asia', 'South Asia', 'Japan/Korea',
                         'Oceania', 'West Asia', 'Africa', 'Europe', 'South America',
                         'Central America', 'Northeast USA', 'Midwest USA', 'South USA',
                         'West USA', 'Canada', 'Seattle']
    style = color_df.set_index('region')
    introductions_by_region = dataframe['region'].value_counts()
    introduction_df = introductions_by_region.to_frame(name = 'num_introductions')
    plot_df = style.join(introduction_df, how='inner')

    mpl.rcParams['font.weight'] = 110
    mpl.rcParams['axes.labelweight'] = 110
    mpl.rcParams['font.size'] = 14

    fig, ax = plt.subplots(figsize=(10,8), facecolor='white')
    ax.bar(plot_df.index, plot_df['num_introductions'], color = plot_df['color'])#, color=colors)
    ax.set_xticklabels(labels = plot_df['label'], rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Region")
    ax.set_ylabel("Number of introductions")
    ax.set_title(lineage)
    plt.tight_layout()

    plt.savefig(output, dpi=250)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot genetic distance vs. geographic distance from closest strain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, default="auspice/seattleflu_flu_seasonal_h3n2_genome_1y_tree.json", help="name of JSON tree")
    parser.add_argument('--output-table', type=str, default="analyses/tables/introductions_h3n2.tsv", help = "name of output TSV")
    parser.add_argument('--lineage', type=str, default="h3n2", help = "lineage of tree")
    parser.add_argument('--colors', type=str, default="config/colors.tsv", help = "name of colors TSV")
    parser.add_argument('--output-size-figure', type=str, default="analyses/figures/introductions_size_h3n2.png", help = "name of output for introductions by size figure")
    parser.add_argument('--output-region-figure', type=str, default="analyses/figures/introductions_region_h3n2.png", help = "name of output for introductions by region figure")
    args = parser.parse_args()

# LoadsJSON tree
tree = bt.loadJSON(args.tree, stats=False)

# Fixes region naming
stub_fixed_tree = assign_region_trait_to_stub_nodes(tree)
stub_fixed_tree.root.traits['region'] = 'ancestor'

# Extracts subtrees_list
subtrees = extract_seattle_subtrees(stub_fixed_tree)

# Extracts basal subtrees list
basal_subtrees = get_basal_subtrees(subtrees)

# Makes pandas introduction dateframe and writes out as TSV
introduction_df = create_dataframe(basal_subtrees, args.output_table)

# Writes out histogram of introductions by size
plot_intro_by_size(introduction_df, args.lineage, args.output_size_figure)

# Writes out histogram of introductions by region
plot_intro_by_region(introduction_df, args.lineage, args.colors, args.output_region_figure)
