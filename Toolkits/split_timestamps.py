import os
import cv2
import pandas as pd
import re

path = os.getcwd()


def natural_keys(text):
    return [(lambda t: int(t) if t.isdigit() else t)(c) for c in re.split(r'(\d+)', text)]


def get_num_frames(path):
    video = cv2.VideoCapture(path)
    nframes = int(video.get(cv2.CAP_PROP_FRAME_COUNT))
    video.release()
    return nframes


video_files = sorted([f for f in os.listdir(path) if f.endswith('.avi')], key=natural_keys)

num_frames_dict = {}
fs = 0
for video_file in video_files:
    video_path = os.path.join(path, video_file)
    num_frames = get_num_frames(video_path)
    num_frames_dict[video_file] = num_frames

timestamp_df = pd.read_csv('timestamps.csv')
timestamp_df.rename(columns={'Time Stamp (ms)': 'timestamp'}, inplace=True)
video_dfs = {}

for video_file, num_frames in num_frames_dict.items():
    video_dfs[video_file] = timestamp_df.head(num_frames)
    timestamp_df = timestamp_df.tail(len(timestamp_df) - num_frames)
    timestamp_df["Frame Number"] -= timestamp_df["Frame Number"].min()
    timestamp_df["timestamp"] /= 1000
    timestamp_df.set_index("Frame Number")

for video_file, video_df in video_dfs.items():
    output_file = os.path.splitext(video_file)[0] + '.pkl'
    video_df.to_pickle(output_file)