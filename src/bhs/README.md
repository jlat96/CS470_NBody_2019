This is a parallel implementation of the Barnes-Hut algorithm in C++ with MPI support based on code found [here](https://github.com/barkm/n-body). Very few changes have been made to @barkm's code other than to the standardize the input and output with the rest of our research. The Barnes-Hut algorithm is an nlogn time algorithm for simulating N-body systems. You can read more about the algorithm here [_Barnes-Hut algorithm_](https://www.nature.com/nature/journal/v324/n6096/abs/324446a0.html). inspired by [_Salmon, John K.  Parallel hierarchical N-body methods. Diss. California Institute of Technology, 1991_](http://thesis.library.caltech.edu/6291/). The program can run on 2<sup>d</sup> process using MPI. 

A [simulation](https://www.youtube.com/watch?v=yFQX-5nmXYc) of two spinning disks and a [visualization](https://www.youtube.com/watch?v=KrtevnjgtgM) of the Barnes-Hut tree. 

## Building/Installation

Install a MPI library such as [OpenMPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/).

To build run

``` 
$ mkdir build
$ cd build
$ cmake ..
$ make
```
 
To check if the program runs (on 2 processes):

```
$ mpirun -np 2 ./n_body
```

## Running

To display the command line options run

```
$ mpirun -np 1 ./n_body -?
```

The following command runs a simulation of 1000 bodies for 10000 time steps. The sampled positions of the bodies are stored in positions.txt and the program is run with the verbose flag (-v).

```
$ mpirun -np 2 ./n_body -N 1000 -s 10000 -v -o positions.txt
```