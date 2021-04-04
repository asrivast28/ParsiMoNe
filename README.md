# ParsiMoNe - Parallel Construction of Module Networks
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

ParsiMoNe (**Par**allel Con**s**truct**i**on of **Mo**dule **Ne**tworks) supports learning of module networks in parallel.

## Requirements
* **gcc** (with C++14 support) is used for compiling the project.  
_This project has been tested only on Linux platform, using version [10.1.0](https://gcc.gnu.org/gcc-10/changes.html)._
* **[Boost](http://boost.org/)** libraries are used for parsing the command line options, logging, and a few other purposes.  
_Tested with version [1.74.0](https://www.boost.org/users/history/version_1_74_0.html)._
* **[TRNG](https://www.numbercrunch.de/trng/)** is used for generating pseudo random numbers sequentially and in parallel.  
_Tested with version [4.22](https://github.com/rabauke/trng4/releases/tag/v4.22)._
* **[Armadillo](http://arma.sourceforge.net/)** is used for executing linear algebra operations during consensus clustering.  
_Tested with version [9.800.3](http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz)._
* **[MPI](https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/mpi31-report.htm)** is used for execution in parallel.  
_Tested with [MVAPICH2 version 2.3.3](http://mvapich.cse.ohio-state.edu/static/media/mvapich/mvapich2-2.3.3-userguide.html)._
* **[SCons](http://scons.org/)** is required for building the project.  
_Tested with version [3.1.2](https://scons.org/doc/3.1.2/HTML/scons-user.html)._
* The following repositories are used as submodules:
  * **[BN Utils](https://github.com/asrivast28/bn-utils)** contains common utilities for learning in parallel and scripts for post-processing.  
  * **[mxx](https://gitlab.com/patflick/mxx)** is used as a C++ wrapper for MPI.  
  * **[C++ Utils](https://github.com/asrivast28/cpp-utils)** are used for logging and timing.  

## Building
After the dependencies have been installed, the project can be built as:  
<pre><code>scons
</code></pre>  
This will create an executable named `parsimone`, which can be used for constraint-based structure learning.  
By default, all the paths from the environment in `CPATH` and `LIBRARY_PATH` variables are used as include paths and library paths.  
Path to external includes and libraries at non-default locations can also be specified as:  
<pre><code>scons LOCALINCLUDES=&lt;comma-delimited list of paths&gt; LOCALLIBS=&lt;comma-delimited list of paths&gt;
</code></pre>

#### Debug
For building the debug version of the executable, the following can be executed:
<pre><code>scons DEBUG=1
</code></pre>  
Debug version of the executable is named `parsimone_debug`.

#### Logging
By default, logging is disabled in the release build and enabled in the debug build.
In order to change the default behavior, `LOGGING=[0,1]` argument can be passed to `scons`:  
<pre><code>scons LOGGING=1 # Enables logging in the release build
</code></pre>
Please be aware that enabling logging will affect the performance.

#### Timing
Timing of high-level operations can be enabled by passing `TIMER=1` argument to `scons`.

## Execution
Once the project has been built, please execute the following for more information on all the options that the executable accepts:
<pre><code>./parsimone --help
</code></pre>
For running in parallel, the following can be executed:
<pre><code> mpirun -np 8 ./parsimone ...
</code></pre>  

## Algorithms
Currently, the only supported algorithm for learning module networks is `lemontree` that corresponds to the algorithm by [Bonnet et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983) originally implemented in [_Lemon-Tree_](https://github.com/erbon7/lemon-tree).

## Publication
_Currently under double-blind review._

## Licensing
Our code is licensed under the Apache License 2.0 (see [`LICENSE`](LICENSE)).
