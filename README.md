# Input-specific inhibitory plasticity improves decision accuracy under noise

__Authors: Bo Shen, Soomin Song, Robert Machold, Bernardo Rudy, Paul W. Glimcher, Kenway Louie, and Robert C. Froemke__

__New York University, Grossman School of Medicine__

---

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [License](./LICENSE)

# Overview

This is a repository for a manuscript under review, titled "Input-specific inhibitory plasticity improves decision accuracy under noise
". This is an unpublished work. The contents may change during the revision.

# Repo Contents

- [Simulation code](./Main.m): Matlab code for the modeling part of the paper. To replicate the simulation, please follow the [main file](./Main.m) for replication of singleton results in [Fig3](./Singleton.m), and binary-input results [Figs4-5](./Binary.m).


# System Requirements

## Hardware and Software Requirements

The `Matlab` code was developped with GPU parallel computation. For optimal performance, we recommend a computer with a graphic card. The code was developed on a Windows-11 laptop with the following specs:

RAM: 16.0 GB  
CPU: Intel(R) Core(TM) i7-10870H CPU @ 2.20GHz   
GPU: NVIDIA GeForce RTX 3070 Laptop GPU  

Matlab code was developped with Matlab 2023a (The Mathworks Inc., 2023).

### OS Requirements

The package development version is tested on *Mac OSX* and *Windows* operating systems. However, the 'MATLAB' GPU parrallel simupation part would be only compatible with Windows system where NVIDA CUDA is properly set.



