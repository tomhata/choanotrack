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


def center_frames(
    path_video: str,
    path_csv: str,
    path_output: str = "./output/",
    bottom_pad: int = 0,
    fill: int = 255,
    rotate: bool = False,
    vertical: bool = False,
):
    """Create an image stack, with each image centered around the colony centroid.

    Args:
        path_video (str): Path to input video
        path_csv (str): Path to colony data csv
        path_output (str): Path to output dir for image stack. Defaults to "./output/".
        bottom_pad (int, optional): Remove bottom pixels. Defaults to 0.
        fill (int, optional): Color value to fill empty background. Defaults to 255.
        rotate (bool, optional): Rotate image based on orientation. Defaults to False.
        vertical (bool, optional): Rotate colony to be vertical. Defaults to False.
    """
    df = pd.read_csv(path_csv, index_col=0)
    reader = iio.get_reader(path_video)
    (frame_w, frame_l) = reader.get_meta_data()["source_size"]
    frame_l = frame_l - bottom_pad
    # img_blank = np.uint8(fill) * np.ones([frame_l, frame_w, 3], dtype=np.uint8)
    center_x_frame = int(round(frame_w) / 2)
    center_y_frame = int(round(frame_l) / 2)
    pathlib.Path(path_output).mkdir(exist_ok=True)

    dt = df.loc[df.index[1], "timestamp_s"] - df.loc[df.index[0], "timestamp_s"]
    if vertical:
        tilt = -df.loc[df.index[0], "orientation_rad"] * 180 / np.pi  # vertical colony
    else:
        tilt = 0.0  # don't rotate first frame.

    for idx, img in tqdm(enumerate(reader), total=df.shape[0]):
        if idx in df.index:
            if rotate:
                tilt -= df.loc[idx, "rotation_rad_s"] * dt * 180 / np.pi
            img_cropped = img[:(-bottom_pad), :, :]
            centroid_x = int(round(df.loc[idx]["centroid_x_px"]))
            centroid_y = int(round(df.loc[idx]["centroid_y_px"]))
            x_diff = centroid_x - center_x_frame
            y_diff = centroid_y - center_y_frame
            img_centered = recenter_img(img_cropped, x_diff, y_diff, fill, tilt)
            path_img_out = pathlib.PurePath(path_output, f"{idx:04}.tif")
            iio.imwrite(str(path_img_out), img_centered)
        if idx >= max(df.index):
            break


def recenter_img(
    img: np.ndarray,
    x_diff: int,
    y_diff: int,
    fill: int = 255,
    tilt: float = 0.0,
) -> np.ndarray:
    """Recenter and optionally rotate a single image.

    Args:
        img (np.ndarray): input image.
        x_diff (int): x coordinate colony centroid - frame center in pixels.
        y_diff (int): y coordinate colony centroid - frame center in pixels.
        fill (int, optional): Color value to fill empty background. Defaults to 255.
        tilt (float, optional): Tilt angle to rotate image if not 0.. Defaults to 0..

    Returns:
        np.ndarray: recentered image.
    """
    if len(img.shape) > 2:
        img = img[:, :, 0]
    frame_l = img.shape[0]
    frame_w = img.shape[1]
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

    img_centered = np.uint8(fill) * np.ones([frame_l, frame_w], dtype=np.uint8)
    img_centered[y_min_blank:y_max_blank, x_min_blank:x_max_blank] = img[
        y_min_img:y_max_img, x_min_img:x_max_img
    ]

    if tilt:
        img_centered = np.uint8(
            skimage.transform.rotate(
                img_centered,
                tilt,
                mode="constant",
                cval=fill,
                preserve_range=True,
            )
        )
    return img_centered


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create new image stack centered on colony from original video"
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
        default=0,
        help="number of pad pixels on bottom of video to remove. Defaults to 0",
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
    parser.add_argument(
        "--vertical",
        "-vt",
        action="store_true",
        help="rotate images based on changes in colony orientation if flag is present.",
        required=False,
    )

    args = parser.parse_args()
    center_frames(
        args.video_in,
        args.csv_in,
        path_output=args.output,
        bottom_pad=args.bottom_pad,
        fill=args.fill,
        rotate=args.rotate,
        vertical=args.vertical,
    )
