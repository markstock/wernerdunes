# wernerdunes
Application of Werner's algorithm to modeling of wind-blown dunes

![dunes dem example](example2.png?raw=true "Example of output, digital elevation model of virtual sand dunes")

## Build and run
With GCC, Eigen3, and libpng installed on a Linux machine, this should be as simple as

	make
	./werner.bin -o out
	./werner.bin -h

## Videos
Here's a high-resolution run of the code: https://www.youtube.com/watch?v=-bo7SWY1b-o

## Background
This code is based off the method of Werner BT. 1995. Eolian dunes: computer simulation and attractor interpretation. Geology 23: 1107–1110.

## Thanks
Eigen, CLI11, libpng

Other implementations of this or similar methods are: https://github.com/ebgoldstein/wDune95 https://github.com/weigert/SimpleWindErosion
