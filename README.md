This code is used to check if integrating object observation in the bundle adjustment of SLAM improves the SLAM results. Objects are represented as wireframe models wherein I use the concept of principle component analysis. Each object is represented as a set of ten eigen values characterising the shape of the object which deforms a mean wireframe.

I have added this eigen value(in other words 'shape' of the object) optimization in the bundle adjustment of the pount cloud, wherein the camera poses and the point cloud 3D locations are optimized.
