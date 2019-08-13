# Curve Fitting and PDE Solving with Gradient Descent

To build all programs concurrently, run 'make' without any arguments
'make clean' will remove compiled executables and created data files

## Problem 1 Curve Fitting with Gradient Descent

To build: run the makefile from the command line
$ make curve_fit

To run: make sure there is a "points.txt" file in the same folder as the executable
either run the program with no arguments or with a learning rate and convergence criterion
$ ./curve_fit
OR
$ ./curve_fit 0.0004 5.0E-7

(when no parameters are entered, the default choices are 0.0004 and 5.0E-7)
(caution: the program skips the first header row of the points file)

## Problem 2 Conjugate Gradient Heat Equation Solver

To build conjugate gradient: run the makefile from the command line
$ make conj_grad

To build conjuage gradient with jacobi preconditioner: run the following
$ make precond_conj_grad

For both versions of conjugate gradient:
(to run with default parameters N=100, NT=200, L=4, T=1, a=1)
$ conj_grad
OR
$ precond_conj_grad

(to run with increasing grid sizes ranging from 10 to 200)
$ conj_grad demo
OR
$precond_conj_grad demo

(to run with your own parameters)
$ conj_grad <n> <NT> <L> <T> <a>

## Plotting and Analysis

Plotting was performed in Matlab with various scripts and can be provided upon request
PDF of the analysis was created in Pages

## Authors

* Alexi Chryssanthou