# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:58:14 2024

@author: jason.kong
"""

import os
import logging
import argparse
from PIL import Image

def split_tif(input_path, output_dir):
    """Split a multi-frame TIF file into individual TIF images.
    
    Args:
        input_path: Path to the input multi-frame TIF file.
        output_dir: Directory where the split TIF images will be saved.
    """
    try:
        input_fn, _ = os.path.splitext(os.path.basename(input_path))
        with Image.open(input_path) as tif:
            for i in range(tif.n_frames):
                tif.seek(i)
                output_path = os.path.join(output_dir, f'{input_fn}_{i}.tif')
                tif.save(output_path)
        logging.info(f'Successfully split and saved images to {output_dir}')
    except Exception as e:
        logging.error(f'Error splitting TIF file: {e}')

def main():
    # Configure log
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s'
        )
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Split a multi-frame TIF file into individual TIF images.'
        )
    parser.add_argument(
        '-i', '--input', metavar='Input Path', required=True, 
        help='Path to the input multi-frame TIF file.'
        )
    parser.add_argument(
        '-o', '--output', metavar='Output Directory', required=True, 
        help='Directory where the split TIF images will be saved.'
        )
    
    args = parser.parse_args()
    
    # Normalize to absolute paths
    input_path = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    
    # Validate input path
    if not os.path.isfile(input_path):
        logging.error('The input file does not exist or is not a file.')
        return
    
    if not input_path.lower().endswith(('.tif', '.tiff')):
        logging.error('The input file does not have a .tif or .tiff extension.')
        return
    
    # Validate and create output directory if necessary
    if not os.path.isdir(output_dir):
        logging.info(f'The output directory does not exist. Creating new directory at {output_dir}.')
        os.makedirs(output_dir, exist_ok=True)

    # Split file
    logging.info('Splitting file...')
    split_tif(input_path, output_dir)
    logging.info('Complete!\n')

if __name__ == "__main__":
    main()