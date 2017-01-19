# Atom Mapping
Atom Mapping is a ROS package that implements a novel form of 3D mapping, intended for use in real-time robot navigation and planning. In essence, an Atom Map is a generalized occupancy grid that exists in continuous space, and is not restricted to a grid-like structure. Moreover, the general concept is amenable to metrics other than traditional metrics of occupancy, for example it can store the signed distance to an obstacle.

Atom Mapping is developed by **Erik Nelson** and **David Fridovich-Keil**. Erik is currently with [nuro.ai](https://nuro.ai), and David is with the [Hybrid Systems Lab](http://hybrid.eecs.berkeley.edu) and the [Berkeley Artificial Intelligence Research lab](http://bair.berkeley.edu).

## Status
AtomMap is set up as a single ROS package that can be downloaded and integrated into existing ROS projects with ease. We also provide a sample client as an example of how AtomMap could be used.

Although the code is well-documented, we encourage you to read our ICRA 2017 paper, linked [here](http://people.eecs.berkeley.edu/~dfk/atom_map_final.pdf).

## Build
To build AtomMap you will need to have ROS installed on your machine. Clone the repository and build as follows:

```
git clone https://github.com/ucberkeley-vip/atom_mapping
cd atom_mapping
cd internal
catkin_make -DCMAKE_BUILD_TYPE=Release
```

## Usage
The last command above will build both the `atom_map` package and the `atom_map_example` package. In order to run the example, you will also need to have installed the [BLAM](https://github.com/erik-nelson/blam) package.

You may need to adjust the `test.launch` file in `atom_map_example` for your system, and also the launch file for BLAM. Once you have all paths fixed for your system, though, running the example is very simple (using multiple terminal windows):

```
roscore
rviz
roslaunch blam_example test_offline.launch
roslaunch atom_map_example test.launch
```

## Wiki
Please visit our [wiki](https://github.com/ucberkeley-vip/atom_mapping/wiki) for more details.

## Issues
If you run into any issues building or using this package, please feel free to post an [issue](https://github.com/ucberkeley-vip/atom_mapping/issues) and one of us will get back to you as quickly as possible.

## Citation
If you find AtomMap useful in your own work, we request that you cite our [paper](http://people.eecs.berkeley.edu/~dfk/atom_map_final.pdf).
