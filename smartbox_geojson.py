import json
import yaml

from alphashape import alphashape
import numpy as np

import matplotlib.pyplot as plt


# Open the YAML file and load its contents
with open("data/s8-6.yaml", 'r') as file:
    yaml_data = yaml.safe_load(file)


filtered_data = {}
for ant in yaml_data["platform"]["array"]["stations"]["s8-6"]["antennas"]:
    antenna = yaml_data["platform"]["array"]["stations"]["s8-6"]["antennas"][ant]
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
    print(box_points)
    alpha_shape_x = []
    alpha_shape_y = []
    for coord in alpha_shape.exterior.coords:
        alpha_shape_x.append(coord[0])
        alpha_shape_y.append(coord[1])
    print(list(alpha_shape.exterior.coords))

    output_geojson["features"].append(
        {
            "type": "Feature",
            "properties": {
                "name": smartbox,
            },
            "geometry": {
                "type": "Polygon",
                "coordinates": [list(alpha_shape.exterior.coords)],
            },
        }
    )

    plt.scatter(x, y, label=smartbox)
    plt.plot(alpha_shape_x, alpha_shape_y)

plt.gca().set_aspect('equal', adjustable='box')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1.15))
plt.savefig("s8-6.png")

with open("s8-6_geo.json", 'w') as file:
    json.dump(output_geojson, file, indent=4)

import geopandas as gpd
gdf = gpd.read_file("s8-6_geo.json")

# Plot the GeoDataFrame
gdf.plot()
plt.show()