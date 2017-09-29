# SST Collation

Code to generate level 3 sea surface temperature data from L2C granules.

At a very basic level the algorithm has 2 steps:

1. Cloud Masking
2. Collation

### Cloud Masking

The cloud masking step performs a series of masking algorithms to remove
cloudy pixels from the SST data. the cloud masks include:

* ACSPO Cloud mask
* A selection of masks from the polar SPT algorithm
  - Median Mask
  - Standard Deviation Mask
* Nearest Neighbor Mask
  - Removes consecutive pixels whose difference is > T_NN
* Diagonal Mask
  - Removes pixels that are off of a three pixel diagonal
* Brightness Temperature Ratio Masks
* Least Squares Mask
  - Computes a  least square curve through each time seires and removes pixels below the ocmputed curve
* Eigen Mask
  - Multispectral mask which removes pixels with large second eigen values


### Collation Step

The collation step converts the SAMPLING_RATE SST values to an hourly SST product using the following steps:

1. Moving average of cloud masked SST values using a large 3 dimensional window
2. Reinstate misclassified cloud pixels back to ocean pixels by comparing to smooth curve
3. Project the smooth curve onto the masked cloud SST data producing an approximate SST curve
4. Weighted averaging of combined clear-sky and approximated data to produce hourly product