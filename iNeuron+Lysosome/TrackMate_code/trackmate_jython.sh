#!/bin/bash

#SBATCH --job-name=trackmate_jython
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --export=ALL
#SBATCH --time=40:00:00
#SBATCH --mem=60G
#SBATCH --output=./%x-%j.log

# Move into your working directory:
cd /scratch/welchjd_root/welchjd2/shared_data/Xenium_0421_Lyso

# Define the root folder where all FieldXXXX image‐folders live:
ROOT_IMAGE_DIR="path/to/lysosome/movement/images"

# Define where to dump the stacked TIFFs and where to store TrackMate results:
STACKED_DIR="$PWD/Stacked_tifs"
RESULTS_DIR="$PWD/TrackMateResults"

# Path to Fiji and to the parameterized Jython script:
FIJI="$HOME/Fiji.app/ImageJ-linux64"
JYTHON_SCRIPT="$PWD/trackmate_jython_all.py"

# Make sure the “Stacked_tifs” and “TrackMateResults” parent directories exist:
mkdir -p "$STACKED_DIR"
mkdir -p "$RESULTS_DIR"

# Loop over every “FieldXXXX” directory under ROOT_IMAGE_DIR:
#for FIELD_DIR in "$ROOT_IMAGE_DIR"/Field*; do
for FIELD_DIR in "$ROOT_IMAGE_DIR"/Field{0001..0480}; do
    # Extract just “FieldXXXX” (e.g. “Field0001”) from the full path:
    FIELD_NAME=$(basename "$FIELD_DIR")

    echo "=== Processing $FIELD_NAME ==="

    # 1) Build the stack‐filename:
    STACK_PATH="$STACKED_DIR/${FIELD_NAME}_stack.tif"

    # 2) Stack all TIFFs from that folder into one multi‐page TIFF:
    #    (same as in your original TrackMate code for Field0048 :contentReference[oaicite:0]{index=0})
    tiffcp "$FIELD_DIR"/*.tif "$STACK_PATH"
    #    Set the ImageJ axes (TYX):
    tiffset -s 270 "axes=TYX" "$STACK_PATH"
    tiffset -s 270 "ImageJ=1.54g
    images=50
    frames=50
    finterval=1
    loop=false" "$STACK_PATH"

    # 3) Create a subfolder for this field’s TrackMate outputs:
    FIELD_OUTPUT_DIR="$RESULTS_DIR/$FIELD_NAME"
    mkdir -p "$FIELD_OUTPUT_DIR"

    # 4) Launch Fiji headless + Jython, passing in:
    #      argv[1] = path to the stacked TIFF
    #      argv[2] = path to this field’s output folder
    "$FIJI" --ij2 --headless --console --jython "$JYTHON_SCRIPT" \
        "$STACK_PATH" "$FIELD_OUTPUT_DIR"

    echo "=== Finished $FIELD_NAME ==="
done

echo "All fields processed."
