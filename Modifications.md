# Modifying the Seismic.jl package

## Contents

* Introduction to making modifications to the package
* Files organization
* Naming conventions
* Formatting conventions
* Documentation
* Style guidelines
* Tests
* Reproducible examples

## Introduction to making modifications to the package

* We show [here](http://seismic.physics.ualberta.ca/docs/develop_SeismicJulia.pdf) the basics of how to fork the main repository, edit it, commit the changes, and make pull requests.

## Files organization

New files must be placed in the proper location. Source files are organized in the following directories into the /src directory:
* API/c
* Imaging
* Modelling
* Operators
* Plotting
* Processing
* ReadWrite
* Solvers
* Tools
* Utils
* Wavelets
* Windows

Tests (see below) must me placed into /test directory and a line should be added in runtests.jl, while examples must be placed into /examples directory and IJulia examples must be placed in examples/IJulia.

## Naming conventions

* Programs: SeisProgramName (eg SeisRadon)
* Linear operators
* OperationNameOp (eg RadonOp)
* Optimization routines
* SolverName (eg ConjugateGradients)

## Formatting conventions

* General notes: if you want a function to do fancy disk operations, use trace headers, processing on gathers etc, first define a simple function that works on data in memory, then make other function definitions at the end of the file (multiple dispatch) where this function will accept an ASCII-string input and will do input from disk, call your function, then write to disk.  The idea is to always have a very simple matlab like code at least for the first function definition in the file. Multiple dispatch was designed for these situations.

* Programs:
```julia
function foo(in1, in2; parameter1=default1, parameter2=default2, parameterN=defaultN)
  return out
end
```

* Linear operators:
```julia
function foo(m, d, adj=true; parameter1=default1, parameter2=default2, parameterN=defaultN)
  return m or return d
end
```

* Optimization routines:
```julia
function foo(m, d, op, params; parameter1=default1, parameter2=default2, parameterN=defaultN)
  return cost
end
```

## Documentation

There are two components to the documentation:
* Self documentation in the function. This is text (docstrings) that is placed at the top of the function in markdown format. Guidelines for docstrings given [here](http://docs.julialang.org/en/release-0.4/manual/documentation/).
* The website: written in Markdown and converted to HTML using MkDocs. Then HTML is just placed in the website directory.

## Style guidelines

The codes should be written following the following guides:
* [Julia Style.jl guidelines](https://github.com/johnmyleswhite/Style.jl)
* [Style Guide for Python Codes](https://www.python.org/dev/peps/pep-0008/#whitespace-in-expressions-and-statements)

## Tests

* Whenever a program is added, a very small, minimalistic, very fast simple test should be added to the /test directory and to /test/runtests.jl so that one can be sure the code will work for other people when they download the package. These tests are not so much about ensuring the geophysics or math is correct in a program, but making sure that a code you write or modify wont make another program in the package fail (for example if you were to modify conjugate gradients incorrectly it might cause MWNI to fail). This is tested automatically by [Travis CI](https://travis-ci.org/) when a pull request is generated and the owners of the repository will get a notice of whether all tests have passed or not before deciding whether to merge the pull request to the main repository.

## Reproducible examples

* Examples showing the functionalities of the package should be added to the directory /examples. Should be small, fast, and show the basic use of the program. If possible have the example download input data (small) from seismic.physics.ualberta.ca/data/. These tests can be written as IJulia notebooks so people can see the results inside Github and see how to use the package and what it can do. 
* Need examples for 5d interpolation, migration, etc.
