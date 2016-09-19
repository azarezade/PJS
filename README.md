# PJS: Patchwise Joint Sparse Tracking Matlab Framework
![Build Status](https://img.shields.io/teamcity/codebetter/bt428.svg)
![License](https://img.shields.io/badge/license-BSD-blue.svg)

PJS is a Matlab framework for visual tracking.
For more information about PJS tracker check out the [project site](https://azarezade.github.io/PJS).
This framework is designed for research purpose only. Please cite [1,2] if you use this framework.

## Features
* It includes implementations of state-of-the-art visual trackers including APG [3], IVT [4], MIL [5], MTT [6] and our proposed tracker PJSM, PJSS [1,2].
* A collection of standard evaluation measures of visual trackers in bin directory.
* A set of m-files for plotting trackers performance.

## Compile
Install SPAMS [7] and vlfeat [8] libraries and add them to MATLAB path. Then run compile.m.
If you need specific compile options, edit compile_mex_code.m from the `toolbox` folder.
Set parameters in "set_param.m" file. Then Install required dataset e.g. from [here](http://visual-tracking.net/). Change the folder naming and ground-truth as the sample dataset (trellis) provided in `dataset` folder:
* imgs: this folder contains all images of video sequence.
* *\<dataset_name\>\_gt.txt*: each row of this file is the [upper_left_corner_x, upper_left_corner_y, width, height] of the bounding box defining the target.
* *\<dataset_name\>\_gtInterv.txt*: this file contains a number defining the interval in which the ground truth are available. For example use 1 if ground-truth is available for every frame and 5 if ground-truth is available for every five frame.

Main file for running different trackers are named as *main_\<tracker_name\>.m*.
Results are saved in `results` folder.

Implementations of APG, IVT, MIL and MTT are taken from author's site, and only edited to run in our framework.

## References
1. Zarezade A., Rabiee H. R., Soltani-Farani A., and Khajenezhad A., “Patchwise Joint Sparse Tracking with Occlusion Detection”, IEEE Transactions on Image Processing (TIP), 2014. [download](http://ieeexplore.ieee.org/document/6873285/)
2. Soltani-Farani, Ali, Hamid R. Rabiee, and Ali Zarezade. "Collaborating frames: Temporally weighted sparse representation for visual tracking.", IEEE International Conference on Image Processing (ICIP), 2014. [download](http://ieeexplore.ieee.org/document/7025091/)
3. Bao, Chenglong, et al. "Real time robust l1 tracker using accelerated proximal gradient approach." Computer Vision and Pattern Recognition (CVPR), 2012 IEEE Conference on. IEEE, 2012.
4. Ross, David A., et al. "Incremental learning for robust visual tracking." International Journal of Computer Vision 77.1-3 (2008): 125-141.
5. Babenko, Boris, Ming-Hsuan Yang, and Serge Belongie. "Robust object tracking with online multiple instance learning." IEEE Transactions on Pattern Analysis and Machine Intelligence 33.8 (2011): 1619-1632.
6. Zhang, Tianzhu, et al. "Robust visual tracking via multi-task sparse learning." Computer Vision and Pattern Recognition (CVPR), 2012 IEEE Conference on. IEEE, 2012.
7. http://spams-devel.gforge.inria.fr/downloads.html
8. http://www.vlfeat.org/install-matlab.html
