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

drivepoly = gpd.read_file("./data/2023_CAFB_15_Minute_Drive_Areas.geojson")

print(drivepoly.shape)
drivepoly.head()

# %%
drivepoly.explore(column="FacilityOID")

# %%
import pandas as pd

availability = pd.read_csv("./data/Agency Data.csv")
print(availability.shape)
availability.head()

# %%
availability["AgencyRef1"].unique()
# %%
len(availability["AgencyRef1"].unique())

# %%
drivepoly["AgencyRef1"].unique()
# %%
len(drivepoly["AgencyRef1"].unique())

# %%
all(drivepoly.groupby("AgencyRef1").size() == 1)

# %%
all(availability.groupby("AgencyRef1").size() == 1)

# %%
drive_hours = drivepoly.merge(
    availability, on="AgencyRef1", how="left", validate="one_to_one"
)
drive_hours.head()

# %%
drive_hours.shape

# %% count missing with weekday in column
drive_hours["Weekend_AM_y"].isna().sum()
# %%
drive_hours["Weekend_PM_y"].isna().sum()


# %%
for i in drive_hours.columns:
    print(i)

# %%


# %%
# union of polygons counting overlapping areas

gdf = drive_hours[drive_hours["Weekend_AM_y"] == 1]
# keep only polygon and multipolygon geometries
gdf = gdf.explode()
gdf = gdf[gdf["geometry"].geom_type.isin(["Polygon", "MultiPolygon"])]
gdf.to_crs("epsg:32618", inplace=True)
gdf.head()
gdf.explore()


# %%

from shapely.errors import TopologicalError
from shapely import errors as se
from shapely.errors import ShapelyError  # This includes most Shapely-related errors


def overlay_union(df1, df2):
    """
    Overlay Union operation used in overlay function
    """

    dfinter = gpd.overlay(
        df1,
        df2,
        how="intersection",
        make_valid=True,
    )

    # Create a list of all columns that start with "AgencyRef1_"
    agency_ref_columns_dfinter = [
        col for col in dfinter.columns if col.startswith("AgencyRef1_")
    ]

    # Step 1: Combine AgencyRef1_1 and AgencyRef1_2 into a single column AgencyRef1
    # Handle NaN values to avoid 'nan' in the result. Only concatenate non-NaN values.
    dfinter["AgencyRef1"] = dfinter.apply(
        lambda row: ",".join(
            filter(pd.notnull, [row[col] for col in agency_ref_columns_dfinter])
        ),
        axis=1,
    )
    dfinter = dfinter.drop(columns=agency_ref_columns_dfinter)
    # try:
    dfsym1 = gpd.overlay(
        df1,
        df2,
        how="difference",
        make_valid=True,
    )
    dfsym2 = gpd.overlay(
        df2,
        df1,
        how="difference",
        make_valid=True,
    )
    dfunion = pd.concat([dfinter, dfsym1, dfsym2], ignore_index=True, sort=False)

    # dfsym = gpd.overlay(
    #     df1,
    #     df2,
    #     how="symmetric_difference",
    #     make_valid=True,
    # )  # _overlay_symmetric_diff(df1, df2)
    # # same for dfsym
    # agency_ref_columns_dfsym = [
    #     col for col in dfsym.columns if col.startswith("AgencyRef1_")
    # ]
    # dfsym["AgencyRef1"] = dfsym.apply(
    #     lambda row: ",".join(
    #         filter(pd.notnull, [row[col] for col in agency_ref_columns_dfsym])
    #     ),
    #     axis=1,
    # )
    # dfsym = dfsym.drop(columns=agency_ref_columns_dfsym)
    # # merge parts
    # dfunion = pd.concat([dfinter, dfsym], ignore_index=True, sort=False)

    # except Exception as e:  # This will catch any exception, including GEOSException

    # print(f"Caught a topology exception {e} ")
    # dfsym1 = gpd.overlay(
    #     df1,
    #     df2,
    #     how="difference",
    #     make_valid=True,
    # )
    # dfsym2 = gpd.overlay(
    #     df2,
    #     df1,
    #     how="difference",
    #     make_valid=True,
    # )
    # dfunion = pd.concat([dfinter, dfsym1, dfsym2], ignore_index=True, sort=False)

    # keep geometry column last
    columns = list(dfunion.columns)
    columns.remove("geometry")
    columns.append("geometry")
    return dfunion.reindex(columns=columns)


# %%

gdf.reset_index(drop=True, inplace=True)
for index, row in gdf.iterrows():
    print("INDEX:", index)
    # Assuming gdf.geometry.name is correctly pointing to the geometry column.
    row_df = pd.DataFrame(row).transpose()  # Transpose to get the correct shape
    row_gdf = gpd.GeoDataFrame(
        row_df[["AgencyRef1", "geometry"]], geometry=row_df.geometry.name, crs=gdf.crs
    )
    if index == gdf.index[0]:
        unioned = row_gdf.copy()
    else:
        unioned = overlay_union(unioned, row_gdf)  # Union the GeoDataFrames

# %%
unioned.explore()

# %%
# Create a new column 'AgencyRefCount' to count the comma-delimited strings in 'AgencyRef1'
unioned["AgencyRefCount"] = unioned["AgencyRef1"].str.count(",") + 1

dissolved_unioned = unioned.dissolve(by="AgencyRefCount", as_index=False)

# %%

amap = dissolved_unioned.explore(
    column="AgencyRefCount",
    scheme="EqualInterval",
    k=7,
    stroke=False,
    legend_kwds={"caption": "Number of Partner Organizations Weekend_AM_y"},
    style_kwds={
        "fillOpacity": 0.75
    },  # Adjust opacity here; 0 is fully transparent, 1 is fully opaque
    tooltip="AgencyRefCount",
)
# amap
# %%
amap.save("./data/Weekend_AM_y.html")
# %%https://github.com/geopandas/geopandas/issues/2792
import shapely


merged = []

tree = shapely.STRtree(df2.geometry.values)
left, right = tree.query(df1.geometry.values, predicate="intersects")

if len(left):
    pairs = (
        df1.take(left)
        .join(pd.DataFrame({"index_right": right}, index=df1.index.values.take(left)))
        .join(
            df2.rename(columns={"geometry": "geom_right"}),
            on="index_right",
            rsuffix="_2",
        )
    )

    intersections = pairs.copy()
    intersections["geometry"] = shapely.intersection(
        intersections.geometry.values, intersections.geom_right.values
    )
    intersections = intersections.drop(columns=["index_right", "geom_right"])
    merged.append(intersections)

    # aggregate areas in right by unique values of left, then use those to clip
    # areas out of left
    clip_left = gp.GeoDataFrame(
        pairs.groupby(level=0).agg(
            {
                "geom_right": lambda g: shapely.union_all(g) if len(g) > 1 else g,
                **{
                    c: "first"
                    for c in df1.columns
                    if not c in ["index_right", "geom_right"]
                },
            }
        ),
        geometry="geometry",
        crs=df1.crs,
    )
    clip_left["geometry"] = shapely.difference(
        clip_left.geometry.values, clip_left.geom_right.values
    )
    clip_left = clip_left.drop(columns=["geom_right"])
    merged.append(clip_left)

    clip_right = (
        gp.GeoDataFrame(
            pairs.rename(columns={"geometry": "geom_left", "geom_right": "geometry"})
            .groupby(by="index_right")
            .agg(
                {
                    "geom_left": lambda g: shapely.union_all(g) if len(g) > 1 else g,
                    "geometry": "first",
                }
            ),
            geometry="geometry",
            crs=df2.crs,
        )
        .join(df2.drop(columns=["geometry"]))
        .rename(
            columns={
                c: f"{c}_2" if c in df1.columns and c != "geometry" else c
                for c in df2.columns
            }
        )
    )
    clip_right["geometry"] = shapely.difference(
        clip_right.geometry.values, clip_right.geom_left.values
    )
    clip_right = clip_right.drop(columns=["geom_left"])
    merged.append(clip_right)


# add any from left or right data frames that did not intersect
diff_left = df1.take(np.setdiff1d(np.arange(len(df1)), left))
merged.append(diff_left)

diff_right = df2.take(np.setdiff1d(np.arange(len(df2)), right)).rename(
    columns={
        c: f"{c}_2" if c in df1.columns and c != "geometry" else c for c in df2.columns
    }
)
merged.append(diff_right)

# merge all data frames
merged = pd.concat(merged, ignore_index=True)

# push geometry column to the end
merged = merged.reindex(
    columns=[c for c in merged.columns if c != "geometry"] + ["geometry"]
)
