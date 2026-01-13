#!/bin/bash

for ((i=0; i<36; i++)); do
    ../build/lwproj "gpt${i}" "power" "10" "0" "0"
done