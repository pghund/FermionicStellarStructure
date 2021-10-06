# StellarStructure
This repository contains code to compute the structure of certain types of stars for my thesis

The file Williams is trying to recreate the algorithm in P.S. Williams' 1988 paper "analytical solutions for the rotating polytrope n=1". There is some issue with this
because he uses Newton-Raphson to find his coefficients, but he never gives his starting point. I can get close to his numerical results with my code, but not exactly. 

The Hachisu code is an attempt to generalize to multiple fermionic species types the single species code in Izumi Hachisu's 1986 paper 
"A versatile method for obtaining strucures of rapidly rotating stars". This is an iterative method and the code will output the computed potentials and 
densities for each iteration, with an initial starting guess the solution of the nonrotating system. The txt file SurfacePlotter will plot the surface of the densities using gnuplot. Unfortunately, I have not been successful in getting this iteration to converge, I suspect the issue is that the electric potential can change relative value very quickly. The pdf Hachisu Numerics gives a mathematical explantion of my generalization of this code.

Since my generalization of the Hachisu code has not worked, I have tried to solve for the densities using a simple gradient descent method. There is a file in my BosonicStars repository in which I explain how I have applied graident descent to that system; the application to this system is essentially the same. Unfortunately I cannot get this to converge either.  

