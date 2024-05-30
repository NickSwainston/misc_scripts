import json
import yaml
import argparse
import math

from alphashape import alphashape
import numpy as np

import matplotlib.pyplot as plt


def mesh_to_earth_approx(ref, xy):
    lng, lat = ref
    east, north = xy
    earth_radius = 6_378_000

    return [
        lng + math.degrees(east / earth_radius) / math.cos(math.radians(lat)),
        lat + math.degrees(north / earth_radius),
    ]

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Generate a GeoJSON file from a YAML file")
    parser.add_argument("yaml_file", help="The YAML file to process")
    args = parser.parse_args()

    # Open the YAML file and load its contents
    with open(args.yaml_file, 'r') as file:
        yaml_data = yaml.safe_load(file)
    station_name = list(yaml_data["platform"]["array"]["stations"].keys())[0]

    # Grab the position data
    filtered_data = {}
    for ant in yaml_data["platform"]["array"]["stations"][station_name]["antennas"]:
        antenna = yaml_data["platform"]["array"]["stations"][station_name]["antennas"][ant]
        smartbox = antenna["smartbox"]
        if smartbox in filtered_data:
            filtered_data[smartbox][ant] = {
                "east": antenna["location_offset"]["east"],
                "north": antenna["location_offset"]["north"],
                "up": antenna["location_offset"]["up"],

            }
        else:
            filtered_data[smartbox] = {
                ant: {
                    "east": antenna["location_offset"]["east"],
                    "north": antenna["location_offset"]["north"],
                    "up": antenna["location_offset"]["up"],
                }
            }

    # Convert to GeoJSON format and plot
    output_geojson = {
        "type": "FeatureCollection",
        "features": [],
    }
    for smartbox in filtered_data:
        x = []
        y = []
        box_points = []
        for ant in filtered_data[smartbox]:
            x.append(filtered_data[smartbox][ant]["east"])
            y.append(filtered_data[smartbox][ant]["north"])
            box_points.append(
                [
                    filtered_data[smartbox][ant]["east"],
                    filtered_data[smartbox][ant]["north"],
                ]
            )
        alpha = 0.33
        alpha_shape = alphashape(np.array(box_points), alpha)
        alpha_shape_x = []
        alpha_shape_y = []
        alpha_shape_lat_long = []
        for coord in alpha_shape.exterior.coords:
            alpha_shape_x.append(coord[0])
            alpha_shape_y.append(coord[1])
            alpha_shape_lat_long.append(mesh_to_earth_approx(
                (
                    yaml_data["platform"]["array"]["stations"][station_name]["reference"]["longitude"],
                    yaml_data["platform"]["array"]["stations"][station_name]["reference"]["latitude"],
                ),
                (coord[0], coord[1])
            ))

        output_geojson["features"].append(
            {
                "type": "Feature",
                "properties": {
                    "name": smartbox,
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [alpha_shape_lat_long],
                },
            }
        )

        plt.scatter(x, y, label=smartbox)
        plt.plot(alpha_shape_x, alpha_shape_y)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1.15))
    plt.savefig(f"{station_name}.png")

    with open(f"{station_name}_geo.json", 'w') as file:
        json.dump(output_geojson, file, indent=4)

    # import geopandas as gpd
    # gdf = gpd.read_file(f"{station_name}_geo.json")

    # # Plot the GeoDataFrame
    # gdf.plot()
    # plt.show()


if __name__ == "__main__":
    main()