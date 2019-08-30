# Project Overview
Given a pair of stereo images of a desk with groceries, generate virtual image from arbitrary view. 
![Alt text](IMG/L2.jpg?raw=true "Title")
The task can be seperated into three sub-tasks:



## 1. Epipolar Rectification 

The images are rectified, such that the epipolar lines are parallel to each other, which constraints the search of corresponding points on a line instead of the whole image.

### 1.1 Harris Corner Detector
Detect the corner in images according to the change of image brightness with Harris detector.
![Alt text](img/harris.png?raw=true "Title")

### 1.2 Find Correspondences 

![Alt text](img/ncc.png?raw=true "Title")

### 1.3  Fundamental Matrix Estimation based on RanSaC Algorithm

RanSaC algorithm is applied in order to precisely estimate fundamental matrix.  

![Alt text](img/RANSAC.png?raw=true "Title")

### 1.4 Epipolar rectification

With fundamental matrix we are able to apply epipolar rectification.

![Alt text](img/rectification.png?raw=true "Title")

## 2. Disparity Map Generation

A semi-global block matching (SGBM) [1] to accurate dense stereo matching

![Alt text](img/disp.png?raw=true "Title")


## 3. Virtual Image Generation
where depth based image rendering (DIBR) algorithm [2] helps a lot.

![Alt text](img/dibr.png?raw=true "Title")

# References
[1] Hirschmuller H. Stereo Processing by Semiglobal Matching and Mutual Information[J]. IEEE Transactions on Pattern Analysis and Machine Intelligence, 2008, 30(2): 328-341.  
[2] Zinger S S, Do Q L. Free-viewpoint depth image based rendering[J]. Journal of Visual Communication and Image Representation, 2010, 21(5): 533-541.