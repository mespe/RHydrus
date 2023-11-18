# HYDRUS-H1 Model

[Original source](https://www.pc-progress.com/en/Default.aspx?h1d-description)

Fortran source files were converted to C files
using [f2c](http://www.netlib.org/f2c/). *NOTE: Per instructions on
source page, lines depending on MSFLIB (Microsoft PowerStation) were commented out -
should not affect model.*

### Other packages

[hyrdusR](https://github.com/shoebodh/hydrusR)



Originally, we planned to embed the FORTRAN code directly within R and avoid
the stand-alone application. We would then pass the input parameters and data
directly from R objects to the FORTRAN code and run the HYDRUS1D simulations.
This would avoid modifying the input files with the new values of the parameter settings.

Unfortunately, the HYDRUS1D code doesn't simply read the input files and then later compute
derived values. Instead, these steps are understandbly combined. This makes it less feasible
to merely pass the inputs from R.

However, this package does make small but helpful changes to the HYDRUS1D code that is
loaded into R.
The R user can specify the full path to each input file.
This means that they don't have to be called, e.g. SELECTOR.IN
and all have to be in the same directory.
Instead, the names of the input files can have identifying characteristics.
Also, we can use input files from different directories related to do different scenarios
and or projects.

The R user can specify an R function (or a C routine) that is called at the end of each iteration
of the simulation.

There are many opportunities to provide much richer access to the HYDRUS1D code and its
functionality from R.  This will, at times, involve changing some of the HYDRUS1D code.

