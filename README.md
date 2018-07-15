## Point Cloud segmentation using Octree ##

# 1. Prerequisites
I have tested the library in **Ubuntu 16.04**, but it should be easy to compile in other platforms.

## C++11 or C++0x Compiler
We use the new thread and chrono functionalities of C++11.

## Pangolin
I use [Pangolin](https://github.com/stevenlovegrove/Pangolin) for visualization and user interface. Dowload and install instructions can be found at: https://github.com/stevenlovegrove/Pangolin.

## Eigen3
Download and install instructions can be found at: http://eigen.tuxfamily.org. **Required at least 3.1.0**.
  
## Building ##

Find the latest version on [Github](https://github.com/georgebola/Octree):
```
git clone https://github.com/georgebola/Octree
```
mkdir build\
cd build\
cmake ..\
make \

## Testing ##

./sample/sample
