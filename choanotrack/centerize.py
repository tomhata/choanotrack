"""Re-render videos as image stacks centered around colony centroid. Currently does not
adjust for rotation of the colony (only translation).
"""
import imageio as iio
import numpy as np
import os
import pandas as pd


def center_colonies(
    path_video: str,
    path_csv: str,
    path_output: str,
    bottom_pad: int = 20,
    fill: int = 0,
):
    """Create an image stack, with each image centered around the colony centroid.

    Args:
        path_video (str): Path to input video
        path_csv (str): Path to colony data csv
        path_output (str): Path to output directory for image stack
        bottom_pad (int, optional): Pad pixels in video not part of image. Defaults to 20.
        fill (int, optional): Color value to fill empty background. Defaults to 0.
    """
    df = pd.read_csv(path_csv, index_col=0)
    if not os.path.exists(path_output):
        os.makedir(path_output)
    reader = iio.get_reader(path_video)
    for idx, img in enumerate(reader):
        if idx in df.index:
            img_blank = np.uint8(fill) * np.ones(
                [img.shape[0] - bottom_pad, img.shape[0], 3],
                dtype=np.uint8
                )
        




img = np.zeros([100, 100, 3], dtype=np.uint8)
