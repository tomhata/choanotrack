"""
Import masks of choanoflagellate colonies from LaVision Davis software sets to record
their physical measurements and export as csv.
"""

import argparse
from collections import OrderedDict
from datetime import datetime
import json
import lvreader
import numpy as np
import pandas as pd
import pathlib
from scipy import signal
import skimage
from tqdm import tqdm

# properties to record using blob analysis
blob_properties = [
    "area",
    "bbox",
    "major_axis_length",
    "minor_axis_length",
    "centroid",
    "eccentricity",
    "orientation",
]

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
    "rotation_rad_s",
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
    "bbox-1": "bbox_x_min_px",
    "bbox-2": "bbox_y_max_px",
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
    lv_masks = lvreader.read_set(path_lv)
    df_main = pd.DataFrame(columns=list_properties)

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
        largest_blob.name = frame_count
        df_main = pd.concat([df_main, largest_blob.to_frame().transpose()])
    return df_main


def import_dataframe(
    path_input: str, start_frame: int = 0, end_frame: int = -1
) -> pd.DataFrame:
    """Import csv of colony data. Optionally, limit range with start_frame and end_frame
    parameters.

    Args:
        path_input (str): Path to csv
        start_frame (int, optional): First frame to include. Defaults to 0.
        end_frame (int, optional): Last frame to include, or -1 to end. Defaults to -1.

    Returns:
        pd.DataFrame: colony data
    """
    df = pd.read_csv(path_input, index_col=0)
    if end_frame == -1:
        end_frame = len(df)
    df = df.loc[start_frame:end_frame]
    return df


def filter_positions(
    df: pd.DataFrame,
    filt_order: int = 4,
    wn: float = 0.12,
) -> pd.DataFrame:
    """Filter position data using Butterworth lowpass filter and filtfilt.

    Args:
        df (pd.DataFrame): colony position data imported by |import_dataframe|.
        filt_order (int, optional): order of filter. Defaults to 4.
        wn (float, optional): critical frequency. Defaults to 0.12.

    Returns:
        pd.DataFrame: dataframe with updated positions.
    """
    b, a = signal.butter(filt_order, wn)
    df.centroid_y_px = signal.filtfilt(b, a, df.centroid_y_px)
    df.centroid_x_px = signal.filtfilt(b, a, df.centroid_x_px)
    df.major_axis_length_px = signal.filtfilt(b, a, df.major_axis_length_px)
    df.minor_axis_length_px = signal.filtfilt(b, a, df.minor_axis_length_px)

    df = pixels_to_um(df)
    df.loc[df.index[1]::, "rotation_rad_s"] = signal.filtfilt(
        b, a, df.loc[df.index[1]::, "rotation_rad_s"]
    )
    return df


def pixels_to_um(df: pd.DataFrame) -> pd.DataFrame:
    """Convert pixel measurements to um, calculate velocities, and add to dataframe.

    Args:
        blob_series (pd.DataFrame): dataframe of blob data in pixels

    Returns:
        pd.DataFrame: dataframe of blob data with additional calculated data
    """
    last_timestamp = -1
    last_y_um = 0
    last_x_um = 0
    last_orientation = 0
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        df.at[index, "area_um2"] = row.area_px2 * (row.scale_um_px**2)
        df.at[index, "major_axis_length_um"] = (
            row.major_axis_length_px * row.scale_um_px
        )
        df.at[index, "minor_axis_length_um"] = (
            row.minor_axis_length_px * row.scale_um_px
        )
        df.at[index, "centroid_y_um"] = row.centroid_y_px * row.scale_um_px
        df.at[index, "centroid_x_um"] = row.centroid_x_px * row.scale_um_px

        dt = row.timestamp_s - last_timestamp
        dy = df.centroid_y_um[index] - last_y_um
        dx = df.centroid_x_um[index] - last_x_um
        dorientation = df.orientation_rad[index] - last_orientation
        if abs(dorientation) > np.pi / 4:
            if dorientation > 0:
                dorientation = dorientation - np.pi
            else:
                dorientation = np.pi - dorientation

        if last_timestamp == -1:
            df.at[index, "rotation_rad_s"] = 0
            df.at[index, "velocity_x_um_s"] = 0
            df.at[index, "velocity_y_um_s"] = 0
            df.at[index, "velocity_mag_um_s"] = 0
            df.at[index, "velocity_angle_rad"] = 0
        else:
            df.at[index, "rotation_rad_s"] = np.float64(dorientation / dt)
            df.at[index, "velocity_x_um_s"] = np.float64(dx / dt)
            df.at[index, "velocity_y_um_s"] = np.float64(dy / dt)
            df.at[index, "velocity_mag_um_s"] = (
                np.float64(
                    df.velocity_x_um_s[index] ** 2 + df.velocity_y_um_s[index] ** 2
                )
                ** 0.5
            )
            df.at[index, "velocity_angle_rad"] = np.arctan2(
                np.float64(df.velocity_y_um_s[index] / df.velocity_mag_um_s[index]),
                np.float64(df.velocity_x_um_s[index] / df.velocity_mag_um_s[index]),
            )

        last_timestamp = row.timestamp_s
        last_y_um = row.centroid_y_um
        last_x_um = row.centroid_x_um
        last_orientation = row.orientation_rad
    return df


def write_log(dict_params: OrderedDict, path_out: pathlib.PurePath):
    """Write log as json.

    Args:
        dict_params (OrderedDict): parameters to log
        path_out (pathlib.PurePath): path to log file
    """
    if pathlib.Path(path_out).exists():
        with open(path_out) as f:
            data = json.load(f)
        data.update(dict_params)
    else:
        data = dict_params
    with open(path_out, "w+") as fp:
        json.dump(data, fp, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Track colonies using LaVision mask files"
    )
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
        help="output path. Defaults to ./output.csv",
        required=False,
    )
    parser.add_argument(
        "--start",
        "-s",
        type=int,
        default=0,
        help="starting frame. Defaults to 0",
        required=False,
    )
    parser.add_argument(
        "--end",
        "-e",
        type=int,
        default=-1,
        help="end frame. defaults to -1",
        required=False,
    )
    parser.add_argument(
        "--filter",
        "-f",
        type=int,
        default=4,
        help="Butterworth filter order. Defaults to 4.",
        required=False,
    )
    parser.add_argument(
        "--wn",
        "-w",
        type=float,
        default=0.12,
        help="Butterworth filter critical frequency. Defaults to 0.12",
        required=False,
    )
    args = parser.parse_args()

    now = datetime.now()
    dt_str = now.strftime("%Y%m%d%H%M")
    path_out = pathlib.PurePath(args.output)
    dir_out = path_out.parent
    fname = path_out.stem
    fname_suffix = path_out.suffix

    path_out_raw = pathlib.PurePath(dir_out, fname + "_raw" + fname_suffix)
    path_out_filt = pathlib.PurePath(
        dir_out,
        fname
        + "_filt_start_"
        + str(args.start)
        + "_end_"
        + str(args.end)
        + fname_suffix,
    )
    path_log = pathlib.PurePath(dir_out, "choanotrack_log.json")

    df = import_set(args.input)
    df = pixels_to_um(df)
    df.to_csv(path_out_raw, index_label="frame")

    if args.end == -1:
        end_frame = max(df.index)
    else:
        end_frame = args.end

    df = df.loc[args.start: end_frame]
    df = filter_positions(df, args.filter, args.wn)
    df.to_csv(path_out_filt, index_label="frame")

    dict_param_logs = {
        dt_str: OrderedDict(
            {
                "operation": "track",
                "input": args.input,
                "output_raw": str(path_out_raw),
                "output_filtered": str(path_out_filt),
                "start_frame": args.start,
                "end_frame": args.end,
                "butter_filter_order": args.filter,
                "butter_wn": args.wn,
            }
        )
    }
    write_log(dict_param_logs, path_log)
