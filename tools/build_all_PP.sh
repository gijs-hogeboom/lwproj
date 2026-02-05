#!/bin/bash

for gpt in {13..35}; do
	../build/lwproj "r3D${gpt}" "power" 20 0 0
done
