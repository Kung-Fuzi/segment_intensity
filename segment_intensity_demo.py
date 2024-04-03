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

Roadmap:
    - Assess and compare different cell segmentation techniques, including 
      thresholding (e.g. binary, adaptive, Otsu, U-Net), edge detection (e.g. 
      Sobel, Canny), clustering (e.g. K-means, fuzzy C-means), and ML (DeepCell, 
      Cellpose). See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10208095/ and 
      https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.21079 for 
      additional info.
    - Improved image pre-processing (e.g. CLAHE contrast adjustment, convolution 
      v deconvolution, edge detection, background subtraction)
    - Optimize threshold cutoffs (e.g. global, local, Otsu, Voronoi-Otsu)
    - Add post-processing (e.g. region joining, fill holes)
"""

import os
import logging
import argparse
import napari
import numpy as np
from skimage.io import imread
from skimage.filters import gaussian
from skimage.measure import label, regionprops
from skimage.morphology import local_minima
from skimage.segmentation import watershed, relabel_sequential

def cell_segmentation(image, spot_sigma=6, outline_sigma=3, threshold=0.02):
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
        # Formats the image as NumPy array of type 'float64' for increased 
        # precision. Napari and other image viewers interpret TIF files as 
        # NumPy-like arrays with size N >= 2 dimensions. For 2D greyscale 
        # images, files are represented as a 2D matrix of brightness values 
        # across X (height) and Y (width) dimensions. Additional dimensions 
        # can be added, such as a C (channel) dimension for RGB images, a 
        # Z (depth) dimension for volumetric images, and a T (time) dimension 
        # for time-lapse images.
        original_image = np.asarray(
            image, 
            dtype=np.float64
            ).copy()
        logging.info(f'\tImage shape: {original_image.shape}')
        logging.info(f'\tMinimum original image intensity: {original_image.min()}')
        logging.info(f'\tMaximum original image intensity: {original_image.max()}')
        logging.info(f'\tAverage original image intensity: {original_image.mean()}')
        
        # Normalizes the image to a standard range.
        normalized_image = (original_image - original_image.min()) / (original_image.max() - original_image.min())
        logging.info(f'\tMinimum normalized image intensity: {normalized_image.min()}')
        logging.info(f'\tMaximum normalized image intensity: {normalized_image.max()}')
        logging.info(f'\tAverage normalized image intensity: {normalized_image.mean()}')

        # Applies two gaussian blurs using the 'gaussian' filter from the 
        # library 'scikit-image'. This filter smooths images by averaging pixel 
        # values with their neighbors weighted using a Guassian function, where 
        # the standard deviation (sigma) controls the degree of smoothing from 
        # the kernel - higher simgas blur over a larger neighborhood, while 
        # smaller sigmas retain more sharpness. Two passes are applied to the 
        # original image in order to optimize for spot seed detection and edge 
        # outline detection, respectively.
        spot_blurred = gaussian(
            normalized_image, 
            sigma=spot_sigma
            )
        outline_blurred = gaussian(
            normalized_image, 
            sigma=outline_sigma if outline_sigma != spot_sigma else spot_sigma
            )
        
        # Processes the blurred images for cell segmentation. First, the 
        # 'local_minima' function from the library 'scikit-image' is used to 
        # identify local minima in the spot_blurred image. These may be 
        # potential cell centers or other important spots in the context of 
        # cell segmentation by surface staining. These spots are subsequently 
        # used in the 'watershed' segmentation algorithm from the library 
        # 'scikit-image' to identify cell boundaries in the 'outline_blurred' 
        # image. The algorithm treats pixel intensity values as a landscape 
        # where high values form peaks and low values form valleys and assigns 
        # every pixel to the nearest seed from 'spots' based on this topology.
        spots = label(local_minima(spot_blurred))
        unique_spots = np.unique(spots)
        logging.info(f'\tNumber of unique spots: {len(unique_spots)}')
        outlines = watershed(
            outline_blurred, 
            spots, 
            mask=(normalized_image > 0)
            )
        unique_outlines = np.unique(outlines)
        logging.info(f'\tNumber of unique outlines: {len(unique_outlines)}')
        
        # Initializes new labels for filtered outlines. This prepares a blank 
        # NumPy array of type 'int32' with the same shape as 'outlines', but 
        # filled with zeros. This will be used to store new labels for segments 
        # that pass the subsequent intensity threshold filter.
        labels = np.zeros_like(outlines, dtype=np.int32)
        
        # Relabels segments sequentially. This loops through each segment 
        # 'region' in 'outlines' and calculates whether the mean intensity is 
        # greater than the threshold. Labels that pass the threshold intensity 
        # filter are saved into 'labels' and relabeled sequentially.
        intensities = []
        for region in regionprops(outlines, intensity_image=normalized_image):
            region_intensity = region.mean_intensity
            intensities.append(region_intensity)
            if region_intensity > threshold:
                labels[outlines == region.label] = region.label
        filtered_labels, _, _ = relabel_sequential(labels)
        unique_filtered = np.unique(filtered_labels)
        logging.info(f'\tNumber of unique labels: {len(unique_filtered)}')
        logging.info(f'\tMinimum label intensity: {min(intensities)}')
        logging.info(f'\tMaximum label intensity: {max(intensities)}')
        logging.info(f'\tAverage label intensity: {sum(intensities) / len(intensities)}')
        logging.info('Successfully segmented image!')
        
    except Exception as e:
        logging.error(f'Error segmenting image: {e}')
    
    return filtered_labels, spot_blurred, outline_blurred, spots, outlines

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
        # Initializes new labels for perimeters. This prepares a blank 
        # NumPy array of type 'float64' with the same shape as 'image', but 
        # filled with zeros. This will be used as a amask to identify the 
        # perimeters of segments.
        perimeters = np.zeros_like(image, dtype=np.float64)
        
        # Calculates the perimeter of each labeled region. This loops through 
        # each segment 'region' in 'labels' and sets the pixels corresponding 
        # to the region's coordinates in the 'perimeters' array to 1, marking 
        # the perimeter of each segment in the image.
        for region in regionprops(labels):
            perimeters[region.coords[:, 0], region.coords[:, 1]] = 1
        
        # Calculates the overall brightness. This calculates the mean intensity 
        # of all pixels in 'image'.
        overall_intensity = np.mean(image)
        logging.info(f'\tAverage overall intensity: {overall_intensity}')
        
        # Calculates the brightness at the edges. This filters 'image' using 
        # the 'perimeter' array mask and calculates the mean intensity of all 
        # perimeter pixels.
        edge_intensity = np.mean(image[perimeters > 0])
        logging.info(f'\tAverage edge intensity: {edge_intensity}')
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
        '-i', '--input', metavar='Input Path', required=True, 
        help='Path to input IHC TIF file.'
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
        '-t', '--threshold', metavar='Threshold', type=float, default=0.02, 
        help='Threshold for background cutoff.'
        )   
    
    args = parser.parse_args()
    
    # Normalize to absolute paths
    input_path = os.path.abspath(args.input)
    
    # Validate input directory
    if not os.path.exists(input_path):
        logging.error('The input file does not exist or is not a file.')
        return

    if not input_path.lower().endswith(('.tif', '.tiff')):
        logging.error('The input file does not have a .tif or .tiff extension.')
    
    input_fn, _ = os.path.splitext(os.path.basename(input_path))
    
    # Loads the image from input path using the 'imread' function from the 
    # library 'scikit-image'. This reads the TIF file from 'input_path' into 
    # memory - in this case, it is stored as the variable 'image'.
    input_image = imread(input_path)
    
    # Segments the image using the 'cell_segmentation' function defined above.
    logging.info('Segmenting image...')
    filtered_labels, spot_blurred, outline_blurred, spots, outlines = cell_segmentation(
        input_image, 
        spot_sigma=args.spot, 
        outline_sigma=args.outline, 
        threshold=args.threshold
        )
        
    # Quantifies the average edge intensity using the 'quantify_edge_intensity' 
    # function defined above.
    logging.info('Analyzing image...')
    edge_intensity = quantify_edge_intensity(
        input_image, 
        filtered_labels
        )

    logging.info('Opening in Napari')

    # Open results in Napari
    viewer = napari.Viewer()
    viewer.add_image(input_image, name='Original Image')
    viewer.add_image(spot_blurred, name='Spot Blurred')
    viewer.add_image(outline_blurred, name='Outline Blurred')
    viewer.add_labels(filtered_labels, name='Cells')
    viewer.add_labels(spots, name='Spots')
    viewer.add_labels(outlines, name='Outlines')
    napari.run()
    
    logging.info('Complete!\n')

if __name__ == '__main__':
    main()