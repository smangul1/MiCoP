#!/bin/bash
wget https://zenodo.org/record/6399965/files/indices.zip?download=1
mv indices.zip\?download\=1 data/indices.zip
cd data
unzip indices.zip
rm indices.zip
cd ..
