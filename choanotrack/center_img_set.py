"""Recenter a single image pair using lavision set data for analysis in a singular,
instantaneous PIV measurement. This code co-opts/duplicates existing code.
Candidate for refactoring if robustness becomes necessary.
"""

import argparse
import imageio as iio
import lvreader
import numpy as np
import pandas as pd
import pathlib

from centerize import recenter_img
import track


def center_set(
    path_lv: str,
    curr_x: int,
    curr_y: int,
    dir_out: str = "./output/",
    min_size: int = 0,
    fill: int = 255,
):
    """_summary_

    Args:
        path_lv (str): _description_
        curr_x (int): _description_
        curr_y (int): _description_
        dir_out (str, optional): _description_. Defaults to "./output/".
        min_size (int, optional): _description_. Defaults to 0.
        fill (int, optional): _description_. Defaults to 255.
    """
    p = pathlib.Path(path_lv)
    base_name_out = p.with_suffix("").parts[-1]
    dir_out = dir_out + base_name_out
    path_out_csv = pathlib.Path(dir_out, base_name_out + "_colony.csv")
    pathlib.Path(path_out_csv).parents[0].mkdir(parents=True, exist_ok=True)

    lv_set = lvreader.read_set(path_lv)
    center_x_frame = int(lv_set[0].as_masked_array().mask.shape[1] / 2)
    center_y_frame = int(lv_set[0].as_masked_array().mask.shape[0] / 2)
    tilt = 0
    if curr_x == -1:
        curr_x = center_x_frame
    if curr_y == -1:
        curr_y = center_y_frame

    df = pd.DataFrame(columns=track.list_properties)
    df = track.measure_all_frames(lv_set, min_size=min_size, x=curr_x, y=curr_y)
    df = track.pixels_to_um(df)
    df.to_csv(path_out_csv, index_label="frame")

    dt = df.loc[df.index[1], "timestamp_s"] - df.loc[df.index[0], "timestamp_s"]

    for idx, buffer in enumerate(lv_set):
        tilt -= df.loc[idx, "rotation_rad_s"] * dt * 180 / np.pi
        img = buffer[0].as_masked_array().data
        centroid_x = int(round(df.loc[idx]["centroid_x_px"]))
        centroid_y = int(round(df.loc[idx]["centroid_y_px"]))
        x_diff = centroid_x - center_x_frame
        y_diff = centroid_y - center_y_frame
        img_centered = recenter_img(img, x_diff, y_diff, fill, tilt)
        path_img_out = pathlib.PurePath(dir_out, f"{base_name_out}_{idx:04}.tif")
        iio.imwrite(str(path_img_out), img_centered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create new image stack centered on colony from LaVision set file"
    )
    parser.add_argument(
        "--input",
        "-i",
        help="input path lavision set file",
        type=str,
    )
    parser.add_argument(
        "--output",
        "-o",
        default="./output/",
        help="Path to output directory. Defaults to ./output/",
        type=str,
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
        "--size",
        "-s",
        type=int,
        default=0,
        help="minimum colony size in pixels squared",
        required=False,
    )
    parser.add_argument(
        "--fill",
        "-f",
        default=255,
        help="grayscale value to fill background of centered images. Defaults to 255",
        type=int,
        required=False,
    )

    args = parser.parse_args()
    center_set(
        path_lv=args.input,
        curr_x=args.x,
        curr_y=args.y,
        dir_out=args.output,
        min_size=args.size,
        fill=args.fill,
    )
