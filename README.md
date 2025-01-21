# Digital Image Processing - Lab 1

Introduction to Color Images

Color images are composed of three primary color channels: Red (R), Green (G), and Blue (B). Each channel represents intensity values for its respective color, combining to form a full-color image. Understanding how to process these channels individually is fundamental in digital image processing.

This lab introduces you to:

Reading color images.

Extracting individual color components (R, G, B).


Converting color images to grayscale.




# Digital Image Processing - Lab 2



This lab focuses on implementing Histogram Equalization for grayscale images to enhance their contrast. The key objective is to improve image quality by redistributing the intensity levels more uniformly across the entire dynamic range.





# Steps Performed in Lab-2





-->Loading and Grayscale Conversion:


A color image is loaded using imread().
The image is converted to grayscale using rgb2gray() for easier processing.



-->Histogram Computation:


The intensity levels of the grayscale image are analyzed, and the histogram is computed manually.
Each intensity level's frequency is calculated using nested loops.



-->Cumulative Distribution Function (CDF):


The cumulative distribution function (CDF) of the normalized histogram is computed.
The minimum CDF value (cdf_min) is used to map intensity values to their new equalized levels.



-->Histogram Equalization:


The grayscale image is transformed using the normalized CDF values to produce an equalized image.
This ensures a more uniform distribution of pixel intensities.



-->Visualization:


The original and equalized histograms are plotted.
The original and equalized grayscale images are displayed side-by-side for comparison.
The CDF is overlaid on the histogram for both the original and equalized images to illustrate the redistribution of intensity levels.
