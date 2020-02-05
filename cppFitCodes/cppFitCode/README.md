# C++ code to fit data

The code implements a parallel Nelder-Mead algorithm in C++ to fit data to functions.
The main functions are inside the _hFiles_ folder. There are two main codes,
one named *getData.cpp* and other one named *fitMain.cpp*. The first one allows
you to calculate the best estimation of the correlation function using
a bootstrap resample. The second one fits the data while also bootstrapping the
data provided.

All the input parameters are written inside the file *input.in*.

In order to compile the code you would need a MPI implemetation for C++, such
as _openMPI_ or _intel_mpi_.

```bash
mpic++ getData.cpp/fitMain.cpp -o nameProgram.out
```

And to run it you would use

```bash
mpiexec -n N_proc ./nameProgram.out inputFile.in
```

The script named *convergeFit.sh* allows you to generate an interative
process to fit the data. It uses a initial guess to calculate a first estimation
of the fitted parameters, then it uses those parameters as the initial parameters
for the next iterations. It continues until the process reaches a maximum
number of iterations or a tolerance is achieved.
