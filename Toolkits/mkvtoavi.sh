#!/bin/bash
i=0
for file in *.mkv; do
    ffmpeg -i "$file" "${i}.avi"
    dim = $(ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of csv=s=x:p=0 "${i}.avi")
    if "$dim" != "600x600"; then
        ffmpeg -i "${i}.avi" -s 600x600 "${i}.avi"
        echo "Video resolution adjusted"
    i=$(i+1)
done