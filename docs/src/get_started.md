# Get Started

## Basics

We first have to make JuAFEM.jl available in our code (along with other packages you might use):
```julia
using JuAFEM
```
Following that, the Nodes that contain the coordinates and the element connectivity matrix (in text files) need to be loaded.
```julia

Elements,Nodes=Load("Elements.txt","Nodes.txt")

```
Note that you have to provide the appropriate location of the txt file. Here we have in the current folder

To view the mesh you can refer to Plots.jl documentation (It is a handy package used for ploting)

Next step is to get the traditional Finite Element stiffness matrix 
```julia

 kf=fem(Elements,Nodes,Nu,E)

```
Here Nu and E represent poissions ratio and elastic modulus respectively
Calculate the Node-based Smoothing Finite Element Method (NS-FEM) stiffness matrix
```julia

kn=nsfem(Elements,Nodes,0.33,2000000)

```
At last use them to get the stiffness matrix for AFEM
```julia

K=afem(kf,kn,alpha)

```
Alpha here is the optimum value for the algorithm at which it will produce nearly exact results (experimentation indicates that alpha lies between 0.5 and 0.6).

Then just solve the system of linear equations using the defined boundary conditions .

## Input Files

Similar to FEM, all calculations in JuAFEM.jl require an input file to describe the nodal coordinates , element connectivity and mechanical properties .

>Nodes Matrix => A two dimensional matrix containing the nodal coordinates with respect to index as nodal number 

> Elements Matrix => A two dimensional matrix containing the node numbers of the isoparametric element. 