"""
Import masks of choanoflagellate colonies from LaVision Davis software sets to record
their physical measurements and export as csv.
"""

import argparse
import lvreader
import numpy as np
import pandas as pd
import skimage
from tqdm import tqdm

from process import pixels_to_um


# table properties for main dataframe
list_properties = [
    "timestamp_s",
    "area_px2",
    "bbox_y_min_px",
    "bbox_y_max_px",
    "bbox_x_min_px",
    "bbox_x_max_px",
    "major_axis_length_px",
    "minor_axis_length_px",
    "centroid_y_px",
    "centroid_x_px",
    "scale_um_px",
    "eccentricity",
    "orientation_rad",
    "area_um2",
    "major_axis_length_um",
    "minor_axis_length_um",
    "centroid_y_um",
    "centroid_x_um",
    "velocity_y_um_s",
    "velocity_x_um_s",
    "velocity_mag_um_s",
    "velocity_angle_rad",
]

# dict for renaming blob series names
dict_property_renames = {
    "area": "area_px2",
    "bbox-0": "bbox_y_min_px",
    "bbox-1": "bbox_y_max_px",
    "bbox-2": "bbox_x_min_px",
    "bbox-3": "bbox_x_max_px",
    "major_axis_length": "major_axis_length_px",
    "minor_axis_length": "minor_axis_length_px",
    "centroid-0": "centroid_y_px",
    "centroid-1": "centroid_x_px",
    "orientation": "orientation_rad",
}


def import_set(path_lv: str) -> pd.DataFrame:
    """Import lavision mask set. Measure colony as blob for each frame and write to a
    pandas dataframe.

    Args:
        path_lv (str): import path to .set file

    Returns:
        pd.DataFrame: colony data for all frames
    """
    blob_properties = [
        "area",
        "bbox",
        "major_axis_length",
        "minor_axis_length",
        "centroid",
        "eccentricity",
        "orientation",
    ]
    lv_masks = lvreader.read_set(path_lv)
    df_main = pd.DataFrame(columns=list_properties)
    df_main.index.name = "frame"

    for frame_count, buffer in tqdm(enumerate(lv_masks), total=len(lv_masks)):
        timestamp = np.float64(buffer[0].attributes["AcqTimeSeries"][0:-3]) / 1000000
        scale = buffer[0].scales.x.slope * 1000
        mask = buffer[0].as_masked_array().mask
        blobs = skimage.measure.label(mask)
        df = pd.DataFrame(
            skimage.measure.regionprops_table(blobs, properties=blob_properties)
        )
        largest_blob = df.sort_values("area", ascending=False).iloc[0]
        largest_blob = largest_blob.rename(dict_property_renames)
        largest_blob = pd.concat(
            [
                largest_blob,
                pd.Series([timestamp, scale], index=["timestamp_s", "scale_um_px"]),
            ]
        )

        # largest_blob.index.name = "frame"
        largest_blob.name = frame_count
        # df_main = df_main.append(largest_blob)
        df_main = pd.concat([df_main, largest_blob.to_frame().transpose()])
    return df_main


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LaVision mask set to blob tracking")
    parser.add_argument(
        "--input",
        "-i",
        help="input path to .set file",
        type=str,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="./output.csv",
        help="output path",
        required=False,
    )
    args = parser.parse_args()

    df = import_set(args.input)
    df = pixels_to_um(df)
    df.to_csv(args.output)
