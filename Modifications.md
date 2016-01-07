# Modifying the Seismic.jl package

## Contents

*	Introduction to making modifications to the package
*	Naming conventions 
*	Formatting conventions
*	Documentation
*	Reproducible examples for all 
*	Use of internal functions (the package has 3 ricker wavelets, we should eliminate those outside of the wavelets definitions and make sure all codes use this one)
*	Examples

## Introduction to making modifications to the package

*	We should show how to fork the repository, make pull requests, and tips for making a pull request that is sure to be merged without problems. This guide can go into the documentation so future developers will know where to begin.
Naming conventions
*	Programs
*	SeisProgramName (eg SeisRadon)
*	Linear operators
*	OperationNameOp (eg RadonOp)
*	Optimization routines
*	SolverName (eg ConjugateGradients)

## Formatting conventions

*	General notes: if you want a function to do fancy disk operations, use trace headers, processing on gathers etc, first define a simple function that works on data in memory, then make other function definitions at the end of the file (multiple dispatch) where this function will accept an asciistring input and will do input  from disk ,call your function, then write to disk.  The idea is to always have a very simple matlab like code at least for the first function definition in the file. Multiple dispatch was designed for these situations.
*	Programs
function foo(in;parameter1=default1,parameter2-default2,parameterN=defaultN)
     return out
end
*	Linear operators
function foo(m,d,adj=true;parameter1=default1,parameter2=default2,parameterN=defaultN)
end

*	Optimization routines
function foo(m,d,op,params;parameter1=default1,parameter2=default2,parameterN=defaultN)
	return cost
end

## Documentation

There are two components to the documentation
*	Self documentation in the function. This is text that is placed at the top of the function in markdown format. 
*	The website: written in markdown (source files here) then converted to html using MkDocs. Then html is just placed in the website directory.
Examples
*	Whenever a program is added, a very small, very fast test should be added to the test directory and runtests.jl so that we can be sure the code will work for other people when they download the package. These tests are not so much about ensuring the geophysics or math is correct in a program, but making sure that a code you write or modify wont make another program in the package fail (for example if you were to modify conjugate gradients incorrectly it might cause MWNI to fail). This is tested automatically when a pull request is generated and the owner of the repository will get a notice of whether all tests have passed or not before deciding whether to merge the pull request to the main repository.
*	Need examples for 5d interpolation, migration, etc. Should be small, fast, and show the basic use of the program. If possible have the example download input data (small) from seismic.physics.ualberta.ca/data/. These tests can be written as IJulia notebooks so people can see the results inside Github and see how to use the package and what it can do.

