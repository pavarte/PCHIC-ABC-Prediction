import gcMapExplorer.lib as gmlib
import numpy as np
import os
import argparse
from scipy import sparse
import scipy.io
#from oct2py import octave
from norm_KR_gpu_cupy import *

def get_norm_KR_vector_matlab(ccmap):
    # Discard rows and columns with missing data
    ccmap.make_readable()  # Assuming ccmap has this method to clean data
    bNonZeros = gmlib.ccmapHelpers.get_nonzeros_index(ccmap.matrix)  # Get indices of non-zero rows/cols
    A = (ccmap.matrix[bNonZeros, :])[:, bNonZeros]  # Subset matrix to only non-zero indices
    print("Successfully extracted A matrix for matlab")
    # Calculate normalization vector using original MATLAB code
    if isinstance(A, list):
        A = np.array([float(x) if isinstance(x, (int, float)) else np.nan for x in A])
    print("No weird instances found")
    A = A[:].astype(np.float32)
    print(A.shape)
    B = np.array([
    [0.0, 0.0, 0.17398819, 0.0, 0.0, 0.0, 0.0, 0.0, 0.27768434, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.21460403],
    [0.17398819, 0.0, 0.0, 0.0, 0.81572187, 0.70048771, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.81572187, 0.0, 0.0, 0.0, 0.44153182, 0.0, 0.0, 0.18628054],
    [0.0, 0.0, 0.70048771, 0.0, 0.0, 0.0, 0.0, 0.9058306, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.44153182, 0.0, 0.0, 0.07811406, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.9058306, 0.07811406, 0.0, 0.0, 0.0],
    [0.27768434, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.21460403, 0.0, 0.0, 0.18628054, 0.0, 0.0, 0.0, 0.0, 0.0]
])
    norm_KR(B)
    print("Norm KR test successful")
    vector = norm_KR(A)
    print(vector)
    print("Managed to get KR normalisation vector from python")
    return vector

def save_vector(vector, output_dir, output_file_name):
    # Save the vector to a file in the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create the output directory if it doesn't exist

    output_path = os.path.join(output_dir, output_file_name)  # Save it as a CSV file
    np.savetxt(output_path, vector, delimiter=',')  # Save the numpy array as a CSV

    print(f"Vector saved to {output_path}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process ccmap and save the output vector.')
    parser.add_argument('--ccmap_input', type=str, help='The ccmap input data.')
    parser.add_argument('--output_file_path', type=str, help='Path to where the vector will be saved.')

    # Parse the arguments
    args = parser.parse_args()

    # Create variable for output path
    output_file_name = args.ccmap_input+".norm_vector.txt"
    output_dir = args.output_file_path

    # Assuming ccmap is loaded from a file or passed directly; adjust as needed
    ccmap = gmlib.ccmap.load_ccmap(args.ccmap_input)  # Load ccmap based on your actual data structure
    output_dir = args.output_file_path
    
    # Call the function to compute vector
    vector = get_norm_KR_vector_matlab(ccmap)# Call the function to compute vector
    
    # Save the vector
    save_vector(vector, output_dir, output_file_name)
    # Export the ccmap as txt
    ccmap_output = os.path.join(output_dir,args.ccmap_input+".exported.txt")
    gmlib.ccmap.export_cmap(ccmap, ccmap_output, doNotWriteZeros=True)


if __name__ == '__main__':
    main()
