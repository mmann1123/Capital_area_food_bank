import geopandas as gpd
import shapely
import pandas as pd
import numpy as np


def clean_names(df):
    # Create a list of all columns that start with "AgencyRef1_"
    agency_ref_columns_dfinter = [
        col for col in df.columns if col.startswith("AgencyRef1")
    ]

    # Step 1: Combine AgencyRef1_1 and AgencyRef1_2 into a single column AgencyRef1
    # Handle NaN values to avoid 'nan' in the result. Only concatenate non-NaN values.
    df["AgencyRef1"] = df.apply(
        lambda row: ",".join(
            filter(
                pd.notnull,
                [row[col] for col in agency_ref_columns_dfinter],
            )
        ),
        axis=1,
    )
    return df.drop(
        columns=[item for item in agency_ref_columns_dfinter if item != "AgencyRef1"],
        errors="ignore",
    )


def extract_polygons(geom):
    # Check if the geometry is 'Polygon' or 'MultiPolygon' and return it directly
    if geom.geom_type in ["Polygon", "MultiPolygon"]:
        return geom
    # If the geometry is a 'GeometryCollection', extract Polygons and MultiPolygons
    elif geom.geom_type == "GeometryCollection":
        # Use the 'geoms' property to iterate over individual geometries
        polygons = [
            part for part in geom.geoms if part.geom_type in ["Polygon", "MultiPolygon"]
        ]
        # If there are any Polygons or MultiPolygons, return their union. Otherwise, return None.
        if polygons:
            return gpd.GeoSeries(polygons).unary_union
    # Return None for non-polygonal geometries
    return None


def my_union(df1, df2):
    merged = []
    # An STRtree spatial index is created for df2 to efficiently find geometries in df2 that intersect with geometries in df1.
    # The .query() method of the STRtree is used with the intersects predicate to find pairs of intersecting geometries between df1 and df2. The indices of intersecting geometries are stored in left and right.

    tree = shapely.STRtree(df2.geometry.values)
    left, right = tree.query(df1.geometry.values, predicate="intersects")

    if len(left):

        # For each pair of intersecting geometries, a DataFrame is constructed that combines information from both df1 and df2, including their geometries.
        pairs = pd.concat(
            [
                df1.take(left),
                (
                    pd.DataFrame(
                        {"index_right": right}, index=df1.index.values.take(left)
                    )
                ),
            ],
            axis=1,
        ).join(
            df2.rename(columns={"geometry": "geom_right"}).reset_index(drop=True),
            on="index_right",
            rsuffix="_2",
        )
        # For each pair of intersecting geometries, the actual intersection geometry is computed and stored in a new column, replacing the original geometry from df1. Other related columns are prepared for the merged output.
        intersections = pairs.copy()
        intersections["geometry"] = shapely.intersection(
            intersections.geometry.values, intersections.geom_right.values
        )
        intersections = intersections.drop(columns=["index_right", "geom_right"])
        merged.append(clean_names(intersections))

        # The geometries from df1 that intersect with any geometry from df2 are "clipped" by removing the intersecting parts, essentially performing a geometric difference operation. This results in geometries from df1 excluding the parts that overlap with df2.
        clip_left = gpd.GeoDataFrame(
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
        merged.append(clean_names(clip_left))
        # Similarly, geometries from df2 that intersect with df1 are clipped by removing the intersecting parts.
        clip_right = (
            gpd.GeoDataFrame(
                pairs.rename(
                    columns={"geometry": "geom_left", "geom_right": "geometry"}
                )
                .groupby(by="index_right")
                .agg(
                    {
                        "geom_left": lambda g: (
                            shapely.union_all(g) if len(g) > 1 else g
                        ),
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

        merged.append(clean_names(clip_right))

    # Any geometries from df1 and df2 that did not intersect with any geometry from the other DataFrame are identified and prepared to be included in the merged output without modification.
    diff_left = df1.take(np.setdiff1d(np.arange(len(df1)), left))
    merged.append(clean_names(diff_left))

    diff_right = df2.take(np.setdiff1d(np.arange(len(df2)), right)).rename(
        columns={
            c: f"{c}_2" if c in df1.columns and c != "geometry" else c
            for c in df2.columns
        }
    )
    merged.append(clean_names(diff_right))

    # merge all data frames
    merged = pd.concat(merged, ignore_index=True)

    # push geometry column to the end
    merged = merged.reindex(
        columns=[c for c in merged.columns if c != "geometry"] + ["geometry"]
    )

    # reset geometry column
    # merged = merged.set_geometry("geometry")
    # drop all empty geometries
    merged = merged[merged.is_valid]
    return merged
