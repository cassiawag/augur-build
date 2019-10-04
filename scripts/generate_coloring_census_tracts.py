import argparse
import math
import colorsys
import pandas as pd
import matplotlib.colors as colors

def adjust_color_shading(r, g, b, factor):
    # converts rgb to hls in order to affect luminance by an inputed factor. returns new rgb
    h, l, s = colorsys.rgb_to_hls(r , g , b)
    l = max(min(l * factor, 1.0), 0.0)
    rgb = colorsys.hls_to_rgb(h, l, s)
    return rgb

def mapping(coor, label):
    # creates dictionary mapping of lat/long based on a label (i.e region, country, etc.)
    cfile = pd.read_csv(coor, sep = '\t')
    location_dict = {}
    for row in cfile.index:
        if cfile.iloc[row,0] == label:
            location_dict[cfile.iloc[row,1]] = {'latitude': cfile.iloc[row,2], 'longitude': cfile.iloc[row,3]}

    location_df = pd.DataFrame.from_dict(location_dict, orient = 'index')
    location_df.index.name = label
    location_df = location_df[(location_df['latitude'] >= 47.126646) & (location_df['latitude'] <= 48.094264)]
    location_df = location_df[(location_df['longitude'] >= -122.511497) & (location_df['longitude'] <= -121.918611)]
    location_df = location_df.reset_index()
    location_df['lat_rank'] = location_df['latitude'].rank()
    location_df['long_rank'] = location_df['longitude'].rank()
    location_lat_sorted = location_df.sort_values(by ='lat_rank')
    location_long_sorted = location_df.sort_values(by = 'long_rank')
    return location_lat_sorted, location_long_sorted

def color_assignment(lat_df, long_df, label, output_fname):

    # Evenly assigns a color on the Nextstrain 36 color scale
    color_scale = ["#511EA8", "#4928B4", "#4334BF", "#4041C7", "#3F50CC", "#3F5ED0", "#416CCE", "#4379CD", "#4784C7", "#4B8FC1", "#5098B9", "#56A0AF", "#5CA7A4", "#63AC99", "#6BB18E", "#73B583", "#7CB878", "#86BB6E", "#90BC65", "#9ABD5C", "#A4BE56", "#AFBD4F", "#B9BC4A", "#C2BA46", "#CCB742", "#D3B240", "#DAAC3D", "#DFA43B", "#E39B39", "#E68F36", "#E68234", "#E67431", "#E4632E", "#E1512A", "#DF4027", "#DC2F24"]
    counter = 0
    color_boundary = lat_df['lat_rank'].max()/36.0
    for index, row in lat_df.iterrows():
        lat_df.loc[index,'color'] = color_scale[counter]
        if math.floor(row['lat_rank'] % color_boundary) == 0:
            counter += 1

    # the luminance will vary based on a factor ranging from 0.25 to 1.25 depending on the longitude ranking
    # lat_df = lat_df.set_index('location')
    # long_df = long_df.set_index('location')
    for index, row in long_df.iterrows():
        factor = (0.5 * row['long_rank'] / long_df['long_rank'].max()) + 0.75
        print(factor)
        lat_df.loc[index,'shade_factor'] = factor

    # outputs new hex color after appropriate shading
    for index, row in lat_df.iterrows():
        r, g, b = colors.hex2color(lat_df.loc[index,'color'])
        new_rgb = adjust_color_shading(r, g, b, lat_df.loc[index,'shade_factor'])
        new_hex = colors.rgb2hex(new_rgb)
        lat_df.loc[index, 'new_color'] = new_hex.upper()

    # subsets dataframe to only include location and color and prints to tsv
    final_colors_df = lat_df[['location', 'new_color']]
    final_colors_df.insert(0, 'label', label)
    final_colors_df.to_csv(output_fname, sep='\t', index = False)

    return lat_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assigns coloring for census tract",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--coordinates', type=str, default="config/lat_longs.tsv", help="tsv file containing lat/long coordinates and cooresponding region")
    parser.add_argument('--label', type=str, default="location", help="geographic entity of interest like region, country or location")
    parser.add_argument('--output', type=str, default="new-colors.tsv", help="name of output file for colors.tsv")
    args = parser.parse_args()

# makes Dataframe mapping of lat_longs file
lat_df, long_df = mapping(args.coordinates, args.label)

# equally distributes colors among the census tracks and outputs a color TSV file
color_assignment(lat_df, long_df, args.label, args.output)
