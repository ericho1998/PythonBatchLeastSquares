# BLSQ Project
## Eric H, Started June 5, 2021

Attempt with Python to recreate a simple batch least squares estimation from a ENGO 585 lab.

Eventually, more versions of least squares, namely sequential and summation of normals, will be added to this project, to use as a free alternative to MATLAB.

### **What does it do?**
---
This code is meant to perform batch least squares on a set of 150 epochs, with each epoch containing 4 measurements of range to targets situated in a local coordinate system. The targets are at [0,0], [100,0], [100,100] and [0,100].

The initial starting point is [50,50], and the position is static for the first 50 measurements, moving to the right and then up.
The resulting least squares should be capable of modelling the path correctly, and output it similarly to a MATLAB script using the matplotlib library.

### **Bugs**
---
The stationary estimations seem good, but as soon as the positions start changing then the estimation goes wild and it doesn't follow the known path (right 50 m and then up 50 m)
will have to look into documentation further
