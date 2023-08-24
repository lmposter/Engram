#!/bin/bash
rm -f settings.mat
for i in {1..5}; do
    if [ -d "$i" ]; then
        rm -f "$i"/gray*.avi
        rm -rf "$i"/result
        rm -f "$i"/ms.mat
    fi
done
