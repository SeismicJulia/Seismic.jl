Solvers
=======

Solvers are used to invert for a model given data. Their parameters are given by a Dictionary (param = Dict()).

ConjugateGradients
^^^^^^^^^^^^^^^
Non-quadratic regularization with CG-LS. The inner CG routine is taken from Algorithm 2 from Scales, 1987. Make sure linear operator passes the dot product. This program can has two methods: it can act on data and a model stored in memory (for small-scale problems such as multidimensional interpolation of a small blocks of data), or on data and model stored on disk (for large-scale problems such as Least squares migration). For preconditioning, a vector of functions (linear operators) can be passed to the code provided each program follows the structure: 

    function mylinearoperator(m,param) 

      (param["adj"] == false)

      return d

    end

and
    
    function mylinearoperator(d,param) 

      (param["adj"] == true)

      return m

    end

(the above functionality can be written into one code).

or 

    function mylinearoperator(m::ASCIIString,d::ASCIIString,param) 

    ... do disk operations for param["adj"] == false or (param["adj"] == true.

    end

where param=Dict() uses the field "adj"=>true or "adj"=>false, and all other parameters needed by the operator are written into the param dictionary.

DotTest
^^^^^^^^^^^^^^^
Dot product test for a linear operator you wish to use in ConjugateGradients. This will ensure that your code using param["adj"]=>true and param["adj"]=false is behaving as a true forward adjoint pair which is a necessary condition for ConjugateGradients to converge.

DotTest3C
^^^^^^^^^^^^^^^
Dot product test for 3 component data.

CGStep
^^^^^^^^^^^^^^^
Simple arithmetic used in updating array in ConjugateGradients.


