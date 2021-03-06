"""Re-render videos as image stacks centered around colony centroid. Currently does not
adjust for rotation of the colony (only translation).
"""
import argparse
import imageio as iio
import numpy as np
import pandas as pd
import pathlib
import skimage
from tqdm import tqdm


def center_colonies(
    path_video: str,
    path_csv: str,
    path_output: str = "./output/",
    bottom_pad: int = 20,
    fill: int = 255,
    rotate: bool = False,
):
    """Create an image stack, with each image centered around the colony centroid.

    Args:
        path_video (str): Path to input video
        path_csv (str): Path to colony data csv
        path_output (str): Path to output directory for image stack. Defaults to "./output/".
        bottom_pad (int, optional): Pad pixels in video not part of image. Defaults to 20.
        fill (int, optional): Color value to fill empty background. Defaults to 255.
        rotate (bool, optional): Rotate image based on colony orientation. Deafults to False.
    """
    df = pd.read_csv(path_csv, index_col=0)
    reader = iio.get_reader(path_video)
    (frame_w, frame_l) = reader.get_meta_data()["source_size"]
    frame_l = frame_l - bottom_pad
    img_blank = np.uint8(fill) * np.ones([frame_l, frame_w, 3], dtype=np.uint8)
    center_x_frame = int(round(frame_w) / 2)
    center_y_frame = int(round(frame_l) / 2)
    pathlib.Path(path_output).mkdir(exist_ok=True)

    dt = df.loc[df.index[1], "timestamp_s"] - df.loc[df.index[0], "timestamp_s"]
    tilt = -df.loc[df.index[0], "orientation_rad"] * 180 / np.pi

    for idx, img in tqdm(enumerate(reader), total=df.shape[0]):
        if idx in df.index:
            img_cropped = img[:(-bottom_pad), :, :]
            centroid_x = int(round(df.loc[idx]["centroid_x_px"]))
            centroid_y = int(round(df.loc[idx]["centroid_y_px"]))
            x_diff = centroid_x - center_x_frame
            y_diff = centroid_y - center_y_frame

            if x_diff > 0:
                x_min_img = x_diff
                x_max_img = frame_w
                x_min_blank = 0
                x_max_blank = frame_w - x_diff
            else:
                x_min_img = 0
                x_max_img = frame_w + x_diff
                x_min_blank = -x_diff
                x_max_blank = frame_w
            if y_diff > 0:
                y_min_img = y_diff
                y_max_img = frame_l
                y_min_blank = 0
                y_max_blank = frame_l - y_diff
            else:
                y_min_img = 0
                y_max_img = frame_l + y_diff
                y_min_blank = -y_diff
                y_max_blank = frame_l

            img_centered = img_blank.copy()
            img_centered[
                y_min_blank:y_max_blank, x_min_blank:x_max_blank
            ] = img_cropped[y_min_img:y_max_img, x_min_img:x_max_img]

            if rotate:
                tilt -= df.loc[idx, "rotation_rad_s"] * dt * 180 / np.pi
                img_centered = np.uint8(
                    skimage.transform.rotate(
                        img_centered,
                        tilt,
                        mode="constant",
                        cval=fill,
                        preserve_range=True,
                    )
                )
            path_img_out = pathlib.PurePath(path_output, f"{idx:04}.tif")
            iio.imwrite(str(path_img_out), img_centered)
        if idx >= max(df.index):
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create new image stack centered on colony"
    )
    parser.add_argument(
        "--video_in",
        "-v",
        help="input path to video",
        type=str,
    )
    parser.add_argument(
        "--csv_in",
        "-c",
        help="Input path to csv with colony data",
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
        "--bottom_pad",
        "-b",
        default=20,
        help="number of pad pixels on bottom of video to remove. Defaults to 20",
        type=int,
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
    parser.add_argument(
        "--rotate",
        "-r",
        action="store_true",
        help="rotate images based on changes in colony orientation if flag is present.",
        required=False,
    )

    args = parser.parse_args()
    center_colonies(
        args.video_in,
        args.csv_in,
        path_output=args.output,
        bottom_pad=args.bottom_pad,
        fill=args.fill,
        rotate=args.rotate,
    )
