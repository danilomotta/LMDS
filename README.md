# MDS and LMDS
A multidimensional scaling (MDS) and landmark multidimensional scaling (LMDS) implementation.

This is an implementation of the technique described in:
An implementation of the paper [Sparse multidimensional scaling using landmark points.](http://graphics.stanford.edu/courses/cs468-05-winter/Papers/Landmarks/Silva_landmarks5.pdf) by De Silva, Vin, and Joshua B. Tenenbaum.

TL;DR: MDS can be used (among other things) to project data to lower dimensions but it needs a distance matrix of NxN, where N is the number of samples. LMDS can project using a matrix of LxN, where L is the number of landmark points and L<<N.

## Files
* mds.py: MDS and LMDS implementation;
* MDS_LMDS_example: notebook with usage examples.

Please contact "ddanilomotta@gmail.com" if you have any questions.