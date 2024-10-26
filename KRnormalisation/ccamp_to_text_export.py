import gcMapExplorer.lib as gmlib
import argparse
import os
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Process ccmap and save the output vector.')
    parser.add_argument('--ccmap_input', type=str, help='The ccmap input data.')
    parser.add_argument('--output_file_path', type=str, help='Path to where the text file of ccmap will be saved.')

    # Parse the arguments
    args = parser.parse_args()

    # Create variable for output path
    ccmap_output = args.ccmap_input+".exported.txt"

    # Load ccmap
    ccmap = gmlib.ccmap.load_ccmap(args.ccmap_input)  # Load ccmap based on your actual data structure

    gmlib.ccmap.export_cmap(ccmap, ccmap_output, doNotWriteZeros=True)


if __name__ == '__main__':   
    main()
