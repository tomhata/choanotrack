"""
Filter raw position data of the colony and recalculate parameters.
"""
import argparse
import numpy as np
import pandas as pd
from scipy import signal
from tqdm import tqdm


def import_dataframe(path_input: str, start_frame=0, end_frame=-1) -> pd.DataFrame:
    """Import csv of colony data. Optionally, limit range with start_frame and end_frame
    parameters.

    Args:
        path_input (str): [description]
        start_frame (int, optional): [description]. Defaults to 0.
        end_frame (int, optional): [description]. Defaults to -1.

    Returns:
        pd.DataFrame: [description]
    """
    df = pd.read_csv(path_input, index_col=0)
    if end_frame == -1:
        end_frame = len(df)
    df = df.loc[start_frame:end_frame]
    return df


def filter_positions(df: pd.DataFrame, filt_order: int, wn: float) -> pd.DataFrame:
    """Filter position data using Butterworth lowpass filter and filtfilt.

    Args:
        df (pd.DataFrame): colony position data imported by |import_dataframe|.
        filt_order (int): order of filter
        wn (float): critical frequency

    Returns:
        pd.DataFrame: dataframe with updated positions.
    """
    b,a = signal.butter(filt_order, wn)
    df.centroid_y_px = signal.filtfilt(b, a, df.centroid_y_px)
    df.centroid_x_px = signal.filtfilt(b, a, df.centroid_x_px)
    print(df.head())
    df = pixels_to_um(df)
    return df


def pixels_to_um(df: pd.DataFrame) -> pd.DataFrame:
    """Convert pixel measurements to um, calculate velocities, and add to dataframe.

    Args:
        blob_series (pd.DataFrame): dataframe of blob data in pixels

    Returns:
        pd.DataFrame: dataframe of blob data with additional calculated data
    """
    last_y_um = -1
    last_x_um = -1
    last_timestamp = -1
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        df.at[index, "area_um2"] = row.area_px2 * (row.scale_um_px ** 2)
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

        if index == 0:
            df.at[index, "velocity_x_um_s"] = 0
            df.at[index, "velocity_y_um_s"] = 0
            df.at[index, "velocity_mag_um_s"] = 0
            df.at[index, "velocity_angle_rad"] = 0
        else:
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
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and recalculate colony data")
    parser.add_argument(
        "--input",
        "-i",
        help="input path to .csv",
        type=str,
    )
    parser.add_argument(
        "--filter",
        "-f",
        help="Butterworth filter order",
        type=int,
    )
    parser.add_argument(
        "--wn",
        "-w",
        help="Butterworth filter critical frequency",
        type=float
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="./output_filtered.csv",
        help="output path",
        required=False,
    )
    parser.add_argument(
        "--start",
        "-s",
        type=int,
        default=0,
        help="starting frame",
        required=False,
    )
    parser.add_argument(
        "--end",
        "-e",
        type=int,
        default=-1,
        help="end frame",
        required=False,
    )
    args = parser.parse_args()

    df = import_dataframe(args.input, start_frame=args.start, end_frame=args.end)
    df = filter_positions(df, args.filter, args.wn)
    df.to_csv(args.output, index_label="frame")
