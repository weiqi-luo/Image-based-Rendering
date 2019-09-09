# Project Overview
The task is to generate a virtual image from an __arbitrary__ view angle based on stereo images.
![Alt text](img/L2.jpg?raw=true "Title")

## 1. Epipolar Rectification 

The images are rectified, such that the epipolar lines are parallel to each other, which constraints the search of corresponding points on a line instead of the whole image.

### 1.1 Harris Corner Detector
Detect the corner in images according to the change of image brightness with Harris detector.
![Alt text](img/harris.png?raw=true "Title")

### 1.2 Find Correspondences 

![Alt text](img/ncc.png?raw=true "Title")

### 1.3  Fundamental Matrix Estimation based on RanSaC Algorithm

For a robust estimate of the fundamental matrix $F$, the RanSaC algorithm was applied.

![Alt text](img/RANSAC.png?raw=true "Title")

### 1.4 Epipolar rectification

With the fundamental matrix we are able to apply epipolar rectification.

![Alt text](img/rectification.png?raw=true "Title")

## 2. Disparity Map Generation

Semi-global block matching (SGBM) [1] is used to obtain a dense stereo matching.

![Alt text](img/disp.png?raw=true "Title")


## 3. Virtual Image Generation
Resulting virtual image based on the Depth Based Image Rendering (DIBR) [2] algorithm.

![Alt text](img/dibr.png?raw=true "Title")

# References
[1] Hirschmüller H. Stereo Processing by Semiglobal Matching and Mutual Information[J]. IEEE Transactions on Pattern Analysis and Machine Intelligence, 2008, 30(2): 328-341.  
[2] Zinger S S, Do Q L. Free-viewpoint depth image based rendering[J]. Journal of Visual Communication and Image Representation, 2010, 21(5): 533-541.