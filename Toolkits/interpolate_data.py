import cv2
from multiprocessing import Pool
import numpy as np
from tifffile import TiffWriter
import os
import argparse
from tqdm import tqdm
from scipy.interpolate import interp1d
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("folder", help="The folder to process.")
    parser.add_argument(
        "--prefix", "-p", default="", help="A prefix to add to processed files."
    )
    parser.add_argument(
        "--sample_rate", "-s", type=int, default=5, help="The sample rate to use."
    )
    parser.add_argument(
        "--behaviour",
        "-b",
        action="store_false",
        help="Whether to enable the behavior feature.",
    )
    parser.set_defaults(behaviour=True)
    return parser.parse_args()


def read_video(fname):
    cap = cv2.VideoCapture(fname)
    frames = []
    while True:
        ret, frame = cap.read()
        if not ret:
            break
        frames.append(frame[..., 0])
    cap.release()
    return frames


def interpolate_video(frames, ts, fs):
    interp_f = interp1d(ts, frames, axis=0)
    ts_new = np.arange(ts[0], ts[-1], 1 / fs)
    frames_new = interp_f(ts_new)
    return frames_new, ts_new[0]


def process_video(vv, tt, folder, fs):
    frames = read_video(os.path.join(folder, vv))
    with open(os.path.join(folder, tt), 'rb') as f:
        ts = pd.read_pickle(f)
        ts = np.array(ts["timestamp"])
    frames_new, ts_start = interpolate_video(frames, ts, fs)
    return frames_new, ts_start, len(frames_new)


def concat_cal_vid(folder, save_folder, cal_list, cal_ts_list, prefix, fs):
    print("Fs:", fs)
    cal_ts = []
    ts_start = []
    save_name = os.path.join(save_folder, prefix + "_calcium_concat.tiff")
    n_frames = 0
    with TiffWriter(save_name, bigtiff=True) as tif, Pool() as pool:
        results = list(
            tqdm(pool.starmap(process_video, [(vv, tt, folder, fs) for vv, tt in zip(cal_list, cal_ts_list)]),
                 desc="combining calcium video"))
        for frames_new, start_ts, n_frames_new in results:
            for ff in range(frames_new.shape[0]):
                tif.write(frames_new[ff, :, :])
                n_frames += 1
            cal_ts.extend(list(np.linspace(start_ts, start_ts + n_frames_new / fs, n_frames_new)))
            ts_start.append(start_ts)
    return n_frames, cal_ts, ts_start


def main():
    num = 0
    output_dir = "output"
    args = get_args()
    save_folder = os.path.join(args.folder, output_dir)

    # Check if the save directory already exists
    while os.path.exists(save_folder):
        num += 1
        save_folder = os.path.join(args.folder, f"{output_dir}{num}")

    os.mkdir(save_folder)

    if num == 0:
        print(f"A directory '{output_dir}' has been created in {os.getcwd()}")
    else:
        print(f"A directory '{output_dir}{num}' has been created in {os.getcwd()}")

    file_list = os.listdir(args.folder)

    # Separate different file types
    cal_list, cal_ts_list = [], []
    for ff in file_list:
        if ff.endswith(".avi"):
            cal_list.append(ff)
        elif ff.endswith(".pkl"):
            cal_ts_list.append(ff)

    # Concatenate calcium video and save timestamps
    n_frames, cal_ts, ts_start = concat_cal_vid(args.folder, save_folder, cal_list, cal_ts_list, args.prefix,
                                                args.sample_rate)
    if np.all(np.sort(cal_ts) == cal_ts) and n_frames == len(cal_ts):
        print("Calcium timestamps check out")
    np.save(os.path.join(save_folder, args.prefix + "_cal_ts.npy"), cal_ts)
    np.save(os.path.join(save_folder, args.prefix + "_start_ts.npy"), ts_start)

    # Concatenate behavior video and save timestamps
    print("Finished concatenating files")


if __name__ == "__main__":
    main()