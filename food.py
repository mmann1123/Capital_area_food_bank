# %% os_prog_v2 environment

# %% [markdown]
# The Capital Area Food Bank procures and distributes food to a partner network made up of hundreds of different organizations across the DMC (soup kitchens, community centers, faith-based organizations, etc). These partner organizations are unique from each other in how they serve clients. They operate under different models, have different levels of capacity, and make food available during different days, different times of day, and different frequencies in a given month. Capital Area Food Bank is trying to understand the level of access that clients have to food assistance network-wide, and where there are gaps in access or saturations/high concentrations of partner sites. This analysis of client access should ideally take into account the time it takes for a client to travel to any partner location, and the method of transportation (walking or public transportation in urban areas, driving in more suburban areas). For this exercise, we have provided an isochrone file that shows areas that are within a 15-minute drive of a CAFB partner location, since we believe clients should not need to travel farther than 15 minutes by car to access the help they need. In more urban areas, clients should be within a 15-minute walk or 20-minutes via public transportation. For simplicity, we have just provided the driving-time layer. We have also included the operating hours of our partners. If a client is within a 15-minute drive of a partner organization, but that organization is only open once a month at a time of day the client cannot attend, we do not equate that to adequate client access.
# Objective Statement
#

# %% [markdown]
#
# We want to understand where there are access gaps in our service area or an oversaturation of resources based on both the geographic access/proximity (15-minute drive) to partner organizations and their hours of operation (considering day of week, time of day - morning, afternoon, evening, and frequency of operations - weekly, biweekly, once a month). With the attached isochrone layer, we can visualize basic client access gaps by seeing where in our service area is not within a 15-minute drive of a CAFB resource, but have no way to assess/incorporate operating hours, or understand the level of access in terms of how many organizations operate in an area. It is our goal to have a visual way to assess client access, incorporating hours of operation and the concentration or lack of partner organizations within a 15-minute drive.

# %%
import geopandas as gpd
from helpers import *

drivepoly = gpd.read_file("./data/2023_CAFB_15_Minute_Drive_Areas.geojson")

print(drivepoly.shape)
drivepoly.head()

# %%
# drivepoly.explore(column="FacilityOID")

# %%
import pandas as pd

availability = pd.read_csv("./data/Agency Data.csv")
print(availability.shape)
availability.head()

availability["AgencyRef1"].unique()
len(availability["AgencyRef1"].unique())


drivepoly["AgencyRef1"].unique()
len(drivepoly["AgencyRef1"].unique())

# %% compare unique values in both dataframes
all(drivepoly.groupby("AgencyRef1").size() == 1)

all(availability.groupby("AgencyRef1").size() == 1)

# %%
drive_hours = drivepoly.merge(
    availability, on="AgencyRef1", how="left", validate="one_to_one"
)
drive_hours.head()


# %% count missing with weekday in column
drive_hours["Weekend_AM_y"].isna().sum()
# %%
drive_hours["Weekend_PM_y"].isna().sum()


# %%  Clean, project data
# select only rows with weekend_am_y == 1
gdf = drive_hours[drive_hours["Weekend_AM_y"] == 1]
# keep only polygon and multipolygon geometries

gdf = gdf[gdf["geometry"].geom_type.isin(["Polygon", "MultiPolygon"])]
gdf.to_crs("epsg:32618", inplace=True)
gdf["geometry"] = gdf.geometry.make_valid()
gdf.head()
# gdf.explore()

# %% union of polygons counting overlapping areas

for index, row in gdf.iterrows():
    print("INDEX:", index)
    # Assuming gdf.geometry.name is correctly pointing to the geometry column.
    row_df = pd.DataFrame(row).transpose()  # Transpose to get the correct shape
    row_gdf = gpd.GeoDataFrame(
        row_df[["AgencyRef1", "geometry"]], geometry=row_df.geometry.name, crs=gdf.crs
    )
    # keep only polygons drop lines and points
    row_gdf["geometry"] = row_gdf["geometry"].apply(extract_polygons)

    if index == gdf.index[0]:
        unioned = row_gdf.copy()
    else:
        # keep only polygons drop lines and points
        unioned["geometry"] = unioned["geometry"].apply(extract_polygons)
        # check validity
        unioned = unioned[unioned.is_valid]
        try:
            unioned = my_union(
                unioned, row_gdf.reset_index(drop=True)
            )  # Union the GeoDataFrames
        except Exception as e:
            print(f"Caught a topology exception {e} ")
            pass

# Create a new column 'AgencyRefCount' to count the comma-delimited strings in 'AgencyRef1'
unioned["AgencyRefCount"] = unioned["AgencyRef1"].str.count(",") + 1
unioned = unioned[unioned.is_valid & ~unioned.is_empty]

dissolved_unioned = unioned.dissolve(by="AgencyRefCount", as_index=False)
dissolved_unioned["geometry"] = dissolved_unioned["geometry"].apply(extract_polygons)

# %% write out file to shapefile
dissolved_unioned.to_file("./data/weekend_am_Y.shp")


# %%
amap = dissolved_unioned.explore(
    column="AgencyRefCount",
    # scheme="EqualInterval",
    # k=1,
    legend_kwds={"caption": "Number of Partner Organizations Weekend_AM_y"},
    style_kwds={
        "fillOpacity": 0.75
    },  # Adjust opacity here; 0 is fully transparent, 1 is fully opaque
    tooltip="AgencyRefCount",
)

# %%
amap
# %%
amap.save("./data/Weekend_AM_y.html")
