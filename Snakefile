path_to_fauna = '../fauna'
if os.environ.get('FAUNA_PATH'):
    path_to_fauna = os.environ.get('FAUNA_PATH')

segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']
lineages = ['h3n2', 'h1n1pdm']
resolutions = ['2y']

def reference_strain(wildcards):
    references = {
        'h3n2': "A/Beijing/32/1992",
        'h1n1pdm': "A/California/07/2009",
        'vic': "B/HongKong/02/1993",
        'yam': "B/Singapore/11/1994"
    }
    return references[wildcards.lineage]

rule all:
    input:
        auspice_tree = expand("auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_tree.json", lineage=lineages, segment=segments, resolution=resolutions),
        auspice_meta = expand("auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_meta.json", lineage=lineages, segment=segments, resolution=resolutions),
        aggregated = expand("results/aggregated/tree-raw_{lineage}_genome_{resolution}.nwk", lineage=lineages, resolution=resolutions)

rule files:
    params:
        seattle_metadata = "metadata/seattle_metadata.tsv",
        outliers = "config/outliers_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        reference = "config/reference_{lineage}_{segment}.gb",
        colors = "config/colors.tsv",
        lat_longs = "config/lat_longs.tsv",
        auspice_config = "config/auspice_config_{lineage}.json",

files = rules.files.params

rule download_background_seqmeta:
    message:
        """
        download_background_seqmeta: Downloading background sequences from fauna
        {wildcards.lineage} {wildcards.segment}
        """
    output:
        seqmeta = "data/background_seqmeta_{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem background_seqmeta_{wildcards.lineage}_{wildcards.segment}
        """

rule parse_background_seqmeta:
    message:
        """
        parse_background_seqmeta: Parsing fasta into sequences and metadata
        {wildcards.lineage} {wildcards.segment}
        """
    input:
        seqmeta = rules.download_background_seqmeta.output.seqmeta
    output:
        sequences = "data/background_sequences_{lineage}_{segment}.fasta",
        metadata = "data/background_metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields = "strain virus isolate_id date region country division location passage authors age sex"
    shell:
        """
        augur parse \
            --sequences {input.seqmeta} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule download_seattle_sequences:
    message:
        """
        download_seattle_sequences: Downloading Seattle sequences from fauna
        {wildcards.lineage} {wildcards.segment}
        """
    output:
        sequences = "data/seattle_sequences_{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus seattle \
            --fasta_fields {params.fasta_fields} \
            --select segment:{wildcards.segment} type:{wildcards.lineage} \
            --path data \
            --fstem seattle_sequences_{wildcards.lineage}_{wildcards.segment}
        """

rule concat_sequences:
    message:
        """
        concat_sequences: Concatenating background and Seattle sequences
        {wildcards.lineage} {wildcards.segment}
        """
    input:
        background_sequences = rules.parse_background_seqmeta.output.sequences,
        seattle_sequences = rules.download_seattle_sequences.output.sequences
    output:
        sequences = "data/sequences_{lineage}_{segment}.fasta"
    shell:
        """
        cat {input.background_sequences} {input.seattle_sequences} > {output.sequences}
        """

rule concat_metadata:
    message:
        """
        concat_metadata: Concatenating background and Seattle metadata
        {wildcards.lineage} {wildcards.segment}
        """
    input:
        background_metadata = rules.parse_background_seqmeta.output.metadata,
        seattle_metadata = files.seattle_metadata
    output:
        metadata = "data/metadata_{lineage}_{segment}.tsv"
    shell:
        """
        python3 scripts/concat_metadata.py \
            --files {input.background_metadata} {input.seattle_metadata} \
            --mergeby strain \
            --fields date region country division residence_census_tract site site_type site_category flu_shot age age_category sex \
            > {output.metadata}
        """

rule filter:
    message:
        """
        filter: Filtering {wildcards.lineage} {wildcards.segment} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region metadata
          - egg-passaged samples
        {wildcards.lineage} {wildcards.segment}
        """
    input:
        sequences = rules.concat_sequences.output.sequences,
        metadata = rules.concat_metadata.output.metadata,
        exclude = files.outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = 800
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude {input.exclude} \
            --exclude-where region=? passage=egg \
            --output {output}
        """

rule select_strains:
    message:
        """
        select_strains: Subsampling background strains, but include all strains with region=seattle
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        sequences = expand("results/filtered_{{lineage}}_{segment}.fasta", segment=segments),
        metadata = expand("data/metadata_{{lineage}}_{segment}.tsv", segment=segments),
        include = files.references
    output:
        strains = "results/strains_{lineage}_{resolution}.txt",
    params:
        viruses_per_month = 80
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {segments} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --resolution {wildcards.resolution} \
            --viruses_per_month {params.viruses_per_month} \
            --output {output.strains}
        """

rule extract:
    message:
        """
        extract: Extract sequences from a given FASTA file that match the given list of sample names
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/extracted_{lineage}_{segment}_{resolution}.fasta'
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """

rule align:
    message:
        """
        align: Aligning sequences to {input.reference}
          - filling gaps with N
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{lineage}_{segment}_{resolution}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads 1
        """

rule tree:
    message:
        """
        tree: Building tree
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{lineage}_{segment}_{resolution}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 1
        """

def clock_rate(w):
    rate = {
        ('h3n2', 'ha'): 0.00356, ('h3n2', 'na'): 0.00298, ('h3n2', 'mp'): 0.00082,
        ('h3n2', 'np'): 0.00140, ('h3n2', 'ns'): 0.00193, ('h3n2', 'pa'): 0.00226,
        ('h3n2', 'pb1'): 0.00177, ('h3n2', 'pb2'): 0.00230,
        ('h1n1pdm', 'ha'): 0.00329, ('h1n1pdm', 'na'): 0.00342, ('h1n1pdm', 'mp'): 0.00209,
        ('h1n1pdm', 'np'): 0.00196, ('h1n1pdm', 'ns'): 0.00278, ('h1n1pdm', 'pa'): 0.00235,
        ('h1n1pdm', 'pb1'): 0.00188, ('h1n1pdm', 'pb2'): 0.00224,
        ('vic', 'ha'): 0.0024, ('vic', 'na'): 0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'): 0.0013
    }
    return rate[(w.lineage, w.segment)]

rule refine:
    message:
        """
        refine: Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.concat_metadata.output.metadata
    output:
        tree = "results/tree_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/branch-lengths_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message:
        """
        ancestral: Reconstructing ancestral sequences and mutations
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{lineage}_{segment}_{resolution}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message:
        """
        translate: Translating amino acid sequences
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{lineage}_{segment}_{resolution}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits:
    message:
        """
        traits: Inferring traits for {params.columns!s}
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_metadata.output.metadata
    output:
        node_data = "results/traits_{lineage}_{segment}_{resolution}.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def _get_clades_file_for_wildcards(wildcards):
    if wildcards.segment == "ha":
        return "config/clades_%s_ha.tsv"%(wildcards.lineage)
    else:
        return "results/clades_%s_ha_%s.json"%(wildcards.lineage, wildcards.resolution)

rule clades:
    message:
        """
        clades: Annotating clades
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = "results/tree_{lineage}_ha_{resolution}.nwk",
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clades = _get_clades_file_for_wildcards
    output:
        node_data = "results/clades_{lineage}_{segment}_{resolution}.json"
    run:
        if wildcards.segment == 'ha':
            shell("""
                augur clades \
                    --tree {input.tree} \
                    --mutations {input.nt_muts} {input.aa_muts} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)
        else:
            shell("""
                python3 scripts/import_tip_clades.py \
                    --tree {input.tree} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)

rule lbi:
    message:
        """
        lbi: Calculating segment LBI
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = 0.3,
        window = 0.5,
        names = "lbi"
    output:
        node_data = "results/lbi_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window} \
            --output {output.node_data}
        """

rule hamming_distance:
    message:
        """
        hamming_distance: Calculating segment hamming distance
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        aligned = rules.align.output.alignment
    output:
        data = "results/distance_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        python3 scripts/reassort/hamming.py --fasta {input.aligned} --output {output.data}
        """

rule clustering:
    message:
        """
        clustering: Identifying clusters based on connected components
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        nt_muts = expand("results/nt-muts_{{lineage}}_{segment}_{{resolution}}.json", segment=segments),
    params:
        cutoff = 30
    output:
        node_data = "results/clustering_{lineage}_{resolution}.json"
    shell:
        """
        python3 scripts/connected_components.py \
            --nt-muts {input.nt_muts} \
            --cutoff {params.cutoff} \
            --output {output.node_data}
        """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades.output.node_data,
        rules.traits.output.node_data,
        rules.lbi.output.node_data,
        rules.clustering.output.node_data
    ]

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export:
    message:
        """
        export: Exporting Ausice JSONs
        {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_metadata.output.metadata,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_tree = "auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_meta.json"
    wildcard_constraints:
        segment = "\S\S"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

checkpoint clusters_fasta:
    message:
        """
        clusters_fasta: Creating directory of fasta files of full-genome clusters
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        clusters = rules.clustering.output.node_data,
        nt_muts = expand("results/nt-muts_{{lineage}}_{segment}_{{resolution}}.json", segment=segments),
        metadata = expand("data/metadata_{{lineage}}_{segment}.tsv", segment = segments)
    params:
        min_size = 2
    output:
        dir = directory("results/clusters/pre_{lineage}_genome_{resolution}")
    shell:
        """
        python3 scripts/extract_cluster_fastas.py \
            --clusters {input.clusters} \
            --nt-muts {input.nt_muts} \
            --metadata {input.metadata} \
            --min-size {params.min_size} \
            --output-dir {output.dir}
        """

rule clusters_intermediate:
    message:
        """
        clusters_intermediate: Copying cluster fastas
        {wildcards.lineage} {wildcards.resolution} {wildcards.cluster}
        """
    input:
        "results/clusters/pre_{lineage}_genome_{resolution}/{cluster}.fasta"
    output:
        "results/clusters/post_{lineage}_genome_{resolution}/{cluster}.fasta"
    shell:
        "cp {input} {output}"

rule reference_genome:
    message:
        """
        reference_genome: Creating full-genome reference genbank file
        {wildcards.lineage}
        """
    input:
        references = expand("config/reference_{{lineage}}_{segment}.gb", segment=segments)
    params:
        ref_strain = reference_strain
    output:
        ref_genome = "config/reference_{lineage}_genome.gb"
    shell:
        """
        python3 scripts/create_ref_genome.py \
            --references {input.references} \
            --ref-strain {params.ref_strain} \
            --output {output.ref_genome}
        """

rule align_clusters:
    message:
        """
        align_clusters: Aligning sequences, keep reference
        {wildcards.lineage} {wildcards.resolution} {wildcards.cluster}
        """
    input:
        sequences = rules.clusters_intermediate.output,
        reference = rules.reference_genome.output.ref_genome
    output:
        alignment = "results/clusters/aligned_{lineage}_genome_{resolution}/{cluster}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads 1
        """

rule tree_clusters:
    message:
        """
        tree_clusters: Building tree for each cluster
        {wildcards.lineage} {wildcards.resolution} {wildcards.cluster}
        """
    input:
        alignment = rules.align_clusters.output.alignment
    output:
        tree = "results/clusters/tree-raw_{lineage}_genome_{resolution}/{cluster}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 1
        """

rule refine_clusters:
    message:
        """
        refine_clusters: Refining tree, reroot to reference, don't produce timetree
        {wildcards.lineage} {wildcards.resolution} {wildcards.cluster}
        """
    input:
        tree = rules.tree_clusters.output.tree,
        alignment = rules.align_clusters.output.alignment,
        metadata = "data/metadata_{lineage}_ha.tsv"
    output:
        tree = "results/clusters/tree_{lineage}_genome_{resolution}/{cluster}.nwk",
        node_data = "results/clusters/branch-lengths_{lineage}_genome_{resolution}/{cluster}.json"
    params:
        root = reference_strain
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root}
        """

def aggregate_trees(wildcards):
    checkpoint_output = checkpoints.clusters_fasta.get(**wildcards).output[0]
    return expand("results/clusters/tree_{lineage}_genome_{resolution}/{cluster}.nwk",
           lineage=wildcards.lineage,
           resolution=wildcards.resolution,
           cluster=glob_wildcards(os.path.join(checkpoint_output, "{cluster}.fasta")).cluster)

rule aggregate_cluster_trees:
    message:
        """
        aggregate_cluster_trees: merging individual cluster trees
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        trees = aggregate_trees
    output:
        tree = "results/aggregated/tree-raw_{lineage}_genome_{resolution}.nwk"
    params:
        outgroup = reference_strain
    shell:
        """
        python3 scripts/merge_trees.py \
            --trees {input.trees} \
            --outgroup {params.outgroup} \
            --output {output.tree}
        """

def aggregate_alignments(wildcards):
    checkpoint_output = checkpoints.clusters_fasta.get(**wildcards).output[0]
    return expand("results/clusters/aligned_{lineage}_genome_{resolution}/{cluster}.fasta",
           lineage=wildcards.lineage,
           resolution=wildcards.resolution,
           cluster=glob_wildcards(os.path.join(checkpoint_output, "{cluster}.fasta")).cluster)

rule align_aggregated:
    message:
        """
        align_aggregated: Concatenates clusters alignment into single fasta file
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        alignments = aggregate_alignments
    output:
        aggregated_alignment = "results/aggregated/aligned_{lineage}_genome_{resolution}.fasta"
    shell:
        "cat {input.alignments} > {output.aggregated_alignment}"

rule refine_aggregated:
    message:
        """
        refine_aggregated: Refining aggregate tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
          - Do not reroot
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        tree = rules.aggregate_cluster_trees.output.tree,
        alignment = rules.align_aggregated.output.aggregated_alignment,
        metadata = "data/metadata_{lineage}_ha.tsv"
    output:
        tree = "results/aggregated/tree_{lineage}_genome_{resolution}.nwk",
        node_data = "results/aggregated/branch-lengths_{lineage}_genome_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --keep-root
        """

rule ancestral_aggregated:
    message:
        """
        ancestral_aggregated: Reconstructing ancestral sequences and mutations for aggregated trees
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        tree = rules.refine_aggregated.output.tree,
        alignment = rules.align_aggregated.output.aggregated_alignment
    output:
        node_data = "results/aggregated/nt-muts_{lineage}_genome_{resolution}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate_aggregated:
    message:
        """
        translate_aggregated: Translating amino acid sequences for aggregated trees
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        tree = rules.refine_aggregated.output.tree,
        node_data = rules.ancestral_aggregated.output.node_data,
        reference = rules.reference_genome.output.ref_genome
    output:
        node_data = "results/aggregated/aa-muts_{lineage}_genome_{resolution}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits_aggregated:
    message:
        """
        traits_aggregated: Inferring ancestral traits for {params.columns!s}
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        tree = rules.refine_aggregated.output.tree,
        metadata = "data/metadata_{lineage}_ha.tsv"
    output:
        node_data = "results/aggregated/traits_{lineage}_genome_{resolution}.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def _get_node_data_for_export_aggregated(wildcards):
    """Return a list of node data files to include for aggregated build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine_aggregated.output.node_data,
        rules.ancestral_aggregated.output.node_data,
        rules.translate_aggregated.output.node_data,
        rules.traits_aggregated.output.node_data
    ]

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export_aggregated:
    message:
        """
        export_aggregated: Exporting Ausice JSONs
        {wildcards.lineage} {wildcards.resolution}
        """
    input:
        tree = rules.refine_aggregated.output.tree,
        metadata = "data/metadata_{lineage}_ha.tsv",
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export_aggregated
    output:
        auspice_tree = "auspice/seattle_flu_seasonal_{lineage}_genome_{resolution}_tree.json",
        auspice_meta = "auspice/seattle_flu_seasonal_{lineage}_genome_{resolution}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
