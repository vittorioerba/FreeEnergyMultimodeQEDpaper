# Free energy landscape of ...

In this repository you can find the code for the paper CITA.

## C++ 

In the subdirectory ```c++```, you can find an implementation in C++ of the gradient descent used in the paper (```gradient.cpp```).
You can also find a launch script for the [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html) (```SLURM_LAUNCH```).

The code uses [EIGEN](http://eigen.tuxfamily.org/index.php?title=Main_Page) to deal with linear algebra.

A small number of helper function implementing the model of the paper are also available.

## Julia

Most of the code used in the paper is written in Julia, and packaged in the package [FreeEnergyMultimodeQED.jl](https://github.com/vittorioerba/FreeEnergyMultimodeQED.jl).
See the README there for install instructions.

## Notebooks

In this subdirectory you will find a Jupyter (Julia) notebook with the code for the data processing and for the data plots.
You will also find a Mathematica notebook used to make the qualitative representations presented in the paper.

Finally, the file ```notebooks/recapNumberSimulations.csv``` contains, for given values of q and p, the total number of gradient descent runs considered.

## Data

The raw and processed data and all the plots produced by the notebooks can be found [here](https://unimi2013-my.sharepoint.com/:f:/g/personal/vittorio_erba_unimi_it/EgDqAdK9WBRFuuI3qoUuJawBhMPM0Ac71lA-scwLNtMTPg).
