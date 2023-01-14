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


def measure_all_frames(
    lv_masks: lvreader.client.set.Set,
    min_size: int = 100,
    x: int = -1,
    y: int = -1,
) -> pd.DataFrame:
    """Measure and record colony data for each frame into a pandas dataframe.

    Args:
        lv_masks (lvreader.client.set.Set): set of lavision frame masks.
        min_size (int, optional): minimum blob size in pixels. Defaults to 100.
        x (int, optional): initial x position (pixels). -1=center. Defaults to -1.
        y (int, optional): initial y position (pixels). -1=center. Defaults to -1.

    Returns:
        pd.DataFrame: colony data for all frames
    """
    df = pd.DataFrame(columns=list_properties)

    (y_center, x_center) = lv_masks[0].as_masked_array().mask.shape
    if x == -1:
        curr_x = round(x_center / 2)
    else:
        curr_x = x
    if y == -1:
        curr_y = round(y_center / 2)
    else:
        curr_y = y

    for frame_count, buffer in tqdm(enumerate(lv_masks), total=len(lv_masks)):
        colony_entry = measure_blob(buffer, curr_x, curr_y, min_size)
        colony_entry.name = frame_count
        df = pd.concat([df, colony_entry.to_frame().transpose()])
    return df


def measure_blob(
    buffer: lvreader.buffer.Buffer,
    curr_x: int,
    curr_y: int,
    min_size: int = 0,
) -> pd.DataFrame:
    """Measure colony data for a single frame buffer and output to a new dataframe.

    Args:
        buffer (lvreader.buffer.Buffer): lavision frame buffer from set.
        curr_x (int):  current x position of colony centroid
        curr_y (int):  current y position of colony centroid
        min_size (int, optional): minimum blob size in pixels. Defaults to 0.

    Returns:
        pd.DataFrame: _description_
    """
    timestamp = np.float64(buffer[0].attributes["AcqTimeSeries"][0:-3]) / 1000000
    scale = buffer[0].scales.x.slope * 1000
    mask = buffer[0].as_masked_array().mask
    blobs = skimage.measure.label(mask)
    df_temp = pd.DataFrame(
        skimage.measure.regionprops_table(blobs, properties=blob_properties)
    )
    if min_size > 0:
        df_temp = df_temp[df_temp["area"] > min_size]
    df_temp.loc[:, "dx_px"] = abs(curr_x - df_temp["centroid-1"])
    df_temp.loc[:, "dy_px"] = abs(curr_y - df_temp["centroid-0"])
    df_temp.loc[:, "dist_px"] = (
        abs(df_temp["dx_px"] ** 2 + df_temp["dy_px"] ** 2) ** 0.5
    )
    colony_entry = df_temp.sort_values("dist_px", ascending=True).iloc[0]
    curr_x = colony_entry["centroid-1"]
    curr_y = colony_entry["centroid-0"]
    colony_entry = colony_entry.rename(dict_property_renames)
    colony_entry = pd.concat(
        [
            colony_entry,
            pd.Series([timestamp, scale], index=["timestamp_s", "scale_um_px"]),
        ]
    )
    return colony_entry


def import_dataframe(
    path_input: str, start_frame: int = 0, end_frame: int = -1
) -> pd.DataFrame:
    """Import csv of colony data. Optionally, limit range with start_frame and end_frame
    parameters.

    Args:
        path_input (str): path to csv
        start_frame (int, optional): first frame to include. Defaults to 0.
        end_frame (int, optional): last frame to include, or -1 to end. Defaults to -1.

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
    df.loc[df.index[1] : :, "rotation_rad_s"] = signal.filtfilt(
        b, a, df.loc[df.index[1] : :, "rotation_rad_s"]
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
        centroid_y_um = row.centroid_y_px * row.scale_um_px
        centroid_x_um = row.centroid_x_px * row.scale_um_px
        df.at[index, "centroid_y_um"] = centroid_y_um
        df.at[index, "centroid_x_um"] = centroid_x_um

        dt = row.timestamp_s - last_timestamp
        dy = centroid_y_um - last_y_um
        dx = centroid_x_um - last_x_um
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
        last_y_um = centroid_y_um
        last_x_um = centroid_x_um
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
        "--x",
        "-x",
        type=int,
        default=-1,
        help="initial colony x position (pixels). defaults to -1 (center)",
        required=False,
    )
    parser.add_argument(
        "--y",
        "-y",
        type=int,
        default=-1,
        help="initial colony y position (pixels). defaults to -1 (center)",
        required=False,
    )
    parser.add_argument(
        "--min_size",
        "-m",
        type=int,
        default=100,
        help="minimum colony size (pixels). defaults to 100",
        required=False,
    )
    parser.add_argument(
        "--filter",
        "-f",
        type=int,
        default=4,
        help="Butterworth filter order. 0 to skip filtering. Defaults to 4.",
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

    lv_masks = lvreader.read_set(args.input)
    df = measure_all_frames(lv_masks, min_size=args.min_size, x=args.x, y=args.y)
    df = pixels_to_um(df)
    df.to_csv(path_out_raw, index_label="frame")
    print(f"Output raw file created: {path_out_filt}")

    if args.end == -1:
        end_frame = max(df.index)
    else:
        end_frame = args.end

    df = df.loc[args.start : end_frame]
    if args.filter > 0:
        df = filter_positions(df, args.filter, args.wn)
        df.to_csv(path_out_filt, index_label="frame")
        print(f"Output filt file created: {path_out_filt}")

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
