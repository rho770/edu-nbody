# A Simple Many Body Code

## Authors

- Thorsten Hater <t.hater@fz-juelich.de>
- Teodor Nikolov

## Introduction 

A very simple code for simulating the interaction of many particles through gravitational (or electrostatic) forces.
The total force acting on a particle is computed by summing up the individual contributions for all other particles.
This leads to $N^2$ complexity in the number of particles $N$.

Your mission is to make this code run as fast as possible!

The reference implementation is in C++11. You can pick and choose which C++11 features you want to use in your implementations. 


## Project organisation

- `src` : contains the vanilla source code
- `handout` : the presentation and handout