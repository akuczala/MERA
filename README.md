# MERA
Matlab code implementing the multiscale entanglement renormalization ansatz (MERA) in (1+1) dimensions. See [here](https://arxiv.org/abs/1109.5334) for algorithm details. As written, the algorithm finds the ground state energy per site of the [tranverse field Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model). It also computes the entanglement entropy and magnetization. 

## Usage
The function `MERA(g,chi,totSteps,stepSize,SIMERA,start)` runs the algorithm. The parameters are as follows:
- `g` : transverse field strength
- `chi` : list of bond dimensions in each layer. The last bond dimension sets either the top or scale invariant bond dimension
- `totSteps` : number of iteration steps
- `stepSize` : number of gradient descent steps per iteration
- `SIMERA` : if `true`, uses the scale-invariant MERA algorithm, otherwise uses a finite number of layers
- `start` : Optionally include initial values for the network tensors

The networks used in this algorithm were constructed with the graphical tensor network tool I wrote, found [here](https://github.com/akuczala/processing-sketches/tree/master/drawNetwork).

This algorithm uses code written by Pfiefer et al to efficently (at the time, in 2014) contract tensors. This code and the corresponding publications can be found [here](https://arxiv.org/abs/1402.0939) and [here](https://arxiv.org/abs/1304.6112).

