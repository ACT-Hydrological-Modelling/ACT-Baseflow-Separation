#!/bin/bash
# Sync Jupytext paired notebooks
# Run this after editing .py files to update .ipynb files, or vice versa

cd "$(dirname "$0")"

# Activate conda environment and sync
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate ACTHYDRO

echo "Syncing Jupytext paired files..."

# Sync from notebooks (this works whether .py files exist or not)
jupytext --sync Baseflow_Separation_LBG.ipynb Baseflow_Separation_SWIOID.ipynb

echo "Done!"

