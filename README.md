# PHSX815 Spring 2021 Project 2

## Monte Carlo Sampling and Integration

This repository contains several types of programs:

- `ExpHypoSim.x` : generates and exports exponentially distributed data samples
with either fixed or gamma-distributed rate parameters [C++]
- `ExpHypoTest.x` : performs analysis on two data sets to determine the measurements required to achieve certain levels of significance in distinguishing them [C++]

### Requirements

In order to compile (by typing `make`) and run the C++ examples, you
need the ROOT package installed (for visualization):
- [ROOT](https://root.cern/) (C++)

### Usage

All of the executables can be called from the
command line with the `-h` or `--help` flag, which will print the options
