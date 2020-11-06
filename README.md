# Free energy landscape of ...

In this repository you can find the code for the paper CITA.

## C++ 

In the subdirectory ```c++```, you can find an implementation in C++ of the gradient descent used in the paper (```gradient.cpp```).
You can also find a launch script for the [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html) (```SLURM_LAUNCH```).

The code uses [EIGEN](http://eigen.tuxfamily.org/index.php?title=Main_Page) to deal with linear algebra.

A small number of helper function implementing the model of the paper are also available.

## Julia

In the subdirectory ```julia```, you can find:
- an implementation in Julia of the gradient descent used in the paper, completely analogous to the C++ implementation;
- an interface to the non-positive-definite QP solver of [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) through [JuMP](https://jump.dev/), used in the paper to perform global optimization of the model;
- an interface to the interior point solver of [IPOPT](https://coin-or.github.io/Ipopt/) through [JuMP](https://jump.dev/), not used in the paper, but possibly useful.

A number of helper function implementing the model of the paper are also available.

