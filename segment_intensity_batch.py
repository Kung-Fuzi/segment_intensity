# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:42:28 2024

@author: jason.kong

Cell Segmentation and Edge Intensity - 
Based off EpiTools for Napari (https://github.com/epitools/epitools)

Dependencies:
    - epitools == 0.0.12
    - napari == 0.4.19
    - numpy == 1.26.4
    - scikit-image == 0.22.0
"""

import os
import csv
import logging
import argparse
import numpy as np
from skimage.io import imread
from skimage.filters import gaussian
from skimage.measure import label, regionprops
from skimage.morphology import local_minima
from skimage.segmentation import watershed, relabel_sequential

def cell_segmentation(image, spot_sigma=6, outline_sigma=3, threshold=0.01):
    """Cell segmentation by thresholded local minima seeded watershed. Code
    derived from EpiTools for Napari (https://github.com/epitools/epitools).
    
    Args:
        image: Input IHC TIF image.
        spot_sigma: Gaussian sigma for controlling spot detection.
        outline_sigma: Gaussian sigma for controlling outline detection.
        threshold: Threshold for background cutoff.
        
    Returns:
        filtered_labels: Cell labels for perimeter and intensity calculations.
    """
    try:
        # Format image as NumPy array
        image = np.asarray(
            image, 
            dtype=np.float64
            )
        
        # Gaussian blur passes for spot and outline detection
        spot_blurred = gaussian(
            image, 
            sigma=spot_sigma
            )
        outline_blurred = gaussian(
            image, 
            sigma=outline_sigma if outline_sigma != spot_sigma else spot_sigma
            )
        
        # Spot and outline detection by local minima and watershed segmentation
        spots = label(local_minima(spot_blurred))
        outlines = watershed(
            outline_blurred, 
            spots, 
            mask=(image > 0)
            )
        
        # Initialize new labels for filtered outlines
        labels = np.zeros_like(outlines, dtype=np.int32)
        
        # Loop through each region to apply the intensity threshold filter
        for region in regionprops(outlines, intensity_image=image):
            if region.mean_intensity > threshold:
                labels[outlines == region.label] = region.label
        
        # Return filtered labels
        filtered_labels, _, _ = relabel_sequential(labels)
        logging.info('Successfully segmented image!')
        
    except Exception as e:
        logging.error(f'Error segmenting image: {e}')
    
    return filtered_labels

def quantify_edge_intensity(image, labels):
    """Quantify cell edge intensity by average intensity per pixel on cell
    perimeters.
    
    Args:
        image: Input IHC TIF image.
        labels: Cell labels for perimeter and intensity calculations.
    Returns:
        edge_intensity: Average intensity per pixel on cell edges
    """
    try:
        # Calculate the perimeter of each labeled region
        perimeters = np.zeros_like(image, dtype=np.float64)
        for region in regionprops(labels):
            perimeters[region.coords[:, 0], region.coords[:, 1]] = 1
        
        # Calculate the brightness at the edges
        edge_intensity = np.mean(image[perimeters > 0])
        logging.info('Successfully analyzed image!')
    
    except Exception as e:
        logging.error(f'Error analyzing image: {e}')
    
    return edge_intensity

def main():
    # Configure log
    logging.basicConfig(
        level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s'
        )
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Segment and quantify cell edge intensities from IHC data.'
        )
    parser.add_argument(
        '-i', '--input', metavar='Input Directory', required=True, 
        help='Directory containing input IHC TIF files.'
        )
    parser.add_argument(
        '-s', '--spot', metavar='Spot Sigma', type=float, default=6, 
        help='Gaussian sigma for controlling spot detection.'
        )
    parser.add_argument(
        '-o', '--outline', metavar='Outline Sigma', type=float, default=3, 
        help='Gaussian sigma for controlling outline detection.'
        )
    parser.add_argument(
        '-t', '--threshold', metavar='Threshold', type=float, default=0.1, 
        help='Threshold for background cutoff.'
        )   
    
    args = parser.parse_args()
    
    # Normalize to absolute paths
    input_dir = os.path.abspath(args.input)
    
    # Validate input directory
    if not os.path.isdir(input_dir):
        logging.error('The input directory does not exist.')
        return
    input_files = [fn for fn in os.listdir(input_dir)]
    
    # Initiate CSV writer
    results_path = os.path.join(input_dir, 'results.csv')
    with open(results_path, 'w', newline='') as results:
        writer = csv.writer(results)
        writer.writerow(['Image Name', 'Average Edge Intensity'])
    
        # Loop over files in input directory
        for input_fn in os.listdir(input_dir):
            
            logging.info(f'Processing image: {input_fn}')
            
            # Validate input path
            if not input_fn.lower().endswith(('.tif', '.tiff')):
                logging.error('The input file does not have a .tif or .tiff extension. Skipping.')
                continue
            
            input_path = os.path.join(input_dir, input_fn)
            
            # Load image from input path
            image = imread(input_path)
            
            # Segment image
            logging.info('Segmenting image...')
            filtered_labels = cell_segmentation(
                image, 
                spot_sigma=args.spot, 
                outline_sigma=args.outline, 
                threshold=args.threshold
                )
            
            # Quantify edge intensity
            logging.info('Analyzing image...')
            edge_intensity = quantify_edge_intensity(
                image, 
                filtered_labels
                )
            
            # Write results to CSV
            writer.writerow([input_fn, edge_intensity])
            
            logging.info(f'Average edge intensity: {edge_intensity}')
            logging.info('Complete!\n')

if __name__ == '__main__':
    main()
