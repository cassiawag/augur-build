import argparse
import pandas as pd

cast_to_int = ["residence_census_tract"]

region_mapping = {
    "Connecticut": "Northeast USA",
    "Maine": "Northeast USA",
    "Massachusetts": "Northeast USA",
    "New Hampshire": "Northeast USA",
    "Rhode Island": "Northeast USA",
    "Vermont": "Northeast USA",
    "New Jersey": "Northeast USA",
    "New York": "Northeast USA",
    "Pennsylvania": "Northeast USA",
    "Illinois": "Midwest USA",
    "Indiana": "Midwest USA",
    "Michigan": "Midwest USA",
    "Ohio": "Midwest USA",
    "Wisconsin": "Midwest USA",
    "Iowa": "Midwest USA",
    "Kansas": "Midwest USA",
    "Minnesota": "Midwest USA",
    "Missouri": "Midwest USA",
    "Nebraska": "Midwest USA",
    "North Dakota": "Midwest USA",
    "South Dakota": "Midwest USA",
    "Delaware": "South USA",
    "Florida": "South USA",
    "Georgia": "South USA",
    "Maryland": "South USA",
    "North Carolina": "South USA",
    "South Carolina": "South USA",
    "Virginia": "South USA",
    "District Of Columbia": "South USA",
    "West Virginia": "South USA",
    "Alabama": "South USA",
    "Kentucky": "South USA",
    "Mississippi": "South USA",
    "Tennessee": "South USA",
    "Arkansas": "South USA",
    "Louisiana": "South USA",
    "Oklahoma": "South USA",
    "Texas": "South USA",
    "Arizona": "West USA",
    "Colorado": "West USA",
    "Idaho": "West USA",
    "Montana": "West USA",
    "Nevada": "West USA",
    "New Mexico": "West USA",
    "Utah": "West USA",
    "Wyoming": "West USA",
    "Alaska": "West USA",
    "California": "West USA",
    "Hawaii": "West USA",
    "Oregon": "West USA",
    "Washington": "West USA",
    "Canada": "Canada",
    "Anguilla": "Central America",
    "Bahamas": "Central America",
    "Barbados": "Central America",
    "Belize": "Central America",
    "Bermuda": "Central America",
    "Cayman Islands": "Central America",
    "Costa Rica": "Central America",
    "Cuba": "Central America",
    "Dominica": "Central America",
    "Dominican Republic": "Central America",
    "El Salvador": "Central America",
    "Grenada": "Central America",
    "Guadeloupe": "Central America",
    "Guatemala": "Central America",
    "Haiti": "Central America",
    "Honduras": "Central America",
    "Jamaica": "Central America",
    "Martinique": "Central America",
    "Mexico": "Central America",
    "Nicaragua": "Central America",
    "Panama": "Central America",
    "Puerto Rico": "Central America",
    "Saint Kitts": "Central America",
    "Saint Lucia": "Central America",
    "Saint Vincent And The Grenadines": "Central America",
    "Trinidad And Tobago": "Central America",
    "Turks And Caicos": "Central America",
    "Usa": "?"
}

def get_region(row):
    region = "?"
    if "country" in row and "division" in row:
        if row["division"] in region_mapping:
            region = region_mapping[row["division"]]
        elif row["country"] in region_mapping:
            region = region_mapping[row["country"]]
        else:
            region = row["region"]
    else:
        region = row["region"]
    return region

def get_location(row):
    location = "other"
    if "residence_census_tract" in row:
        location = str(int(row["residence_census_tract"]))
    return location

def select(file, mergeby, fields):
    entries = pd.read_csv(file, sep='\t')
    mapping = {}
    for index, row in entries.iterrows():
        if row[mergeby]:
            values = []
            for field in fields:
                if field in row:
                    if str(row[field]) != 'nan':
                        if field == "region":
                            values.append(get_region(row))
                        elif field == "location":
                            values.append(get_location(row))
                        else:
                            if field in cast_to_int:
                                values.append(str(int(row[field])))
                            else:
                                values.append(str(row[field]))
                    else:
                        values.append("?")
                else:
                    values.append("?")
            # Remove quotes to avoid read_metadata error later
            values = list(map(lambda x: x.replace('"', ''), values))
            mapping[row[mergeby]] = values
    for key, values in mapping.items():
        print(str(key) + "\t" + "\t".join(values))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Concatenate multiple tsvs, merging specified columns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--files', nargs='+', type=str, required=True, help="list of tsv files")
    parser.add_argument('--mergeby', type=str, help="column name to merge on")
    parser.add_argument('--fields', nargs='+', type=str, help="column names to include")
    args = parser.parse_args()

    print(args.mergeby + "\t" + "\t".join(args.fields))
    for file in args.files:
        select(file, args.mergeby, args.fields)
