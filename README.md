# SignedStabilityBenchmark
Generation of random signed networks with a planted optimal partition and the evaluation of some partitioning methods with respect to the *Correlation Clustering (CC) Problem*

* Copyright 2020-21 Nejat Arınık

*SignedStabilityBenchmark* is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see the file `LICENCE`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/arinik9/SignedStabilityBenchmark
* Contact: Nejat Arınık <arinik9@gmail.com>


-----------------------------------------------------------------------

# Description
This set of `R` and `Julia` scripts was designed to generate random signed networks, where we know their optimal solutions by construction thanks to the definition of [stability range](https://doi.org/10.1145/1553374.1553473), and to apply some heuristic and exact partitioning methods onto these networks by solving the Correlation Clustering Problem. In this repository, we include several exact and heuristic methods that we list below.

## Exact methods

* CoNS(nbMaxEdit, pruningWithoutMVMO):  integrated in the [EnumCC](https://github.com/arinik9/EnumCC) method, which aims to explore the direct neighbor optimal solutions of a given optimal solution up to distance *nbMaxEdit*. See Sections 5.4.2 and 5.6.2 of *[Arınık'21]* for more details.

## Heuristic methods

* **(stochastic) Iterated Local Search:** ILS-CC_l1_a<ALPHA>_g<GAIN>_p3_t<TIME_LIMIT>_i<ITER_NUMBER>  ([available upon request from the authors](https://doi.org/10.1007/s13675-017-0082-6))
* **(stochastic) Greedy Randomized Adaptive Search Procedure:** GRASP-CC_l1_a<ALPHA>_g<GAIN>_t<TIME_LIMIT>_i<ITER_NUMBER>_n1  ([available upon request from the authors](https://doi.org/10.1007/s13675-017-0082-6))
* **(stochastic) Variable neighborhood search:** Brusco-VNS-CC_t<TIME_LIMIT> ([source code](https://github.com/arinik9/HeuristicsCC))
* **(stochastic) Tabu Search:** TS-CC_t<TIME_LIMIT>  ([source code](https://github.com/arinik9/HeuristicsCC))
* **(stochastic) Simulated Annealing:** SA-CC_t<TIME_LIMIT> ([source code](https://github.com/arinik9/HeuristicsCC))
* **(stochastic) Memetic:** MLMSB-CC_r1_t<TIME_LIMIT>  ([available upon request from the authors](https://doi.org/10.1016/j.knosys.2015.05.006))
* **(deterministic) Message Passing, followed by Greedy Additive Edge Contraction, followed by Kernighan-Lin with joins:** MP-GAEC-KLj-CC_t<TIME_LIMIT> ([source code](https://github.com/LPMP/))
* **(deterministic) Iterative Cycle Packing, followed by Greedy Additive Edge Contraction, followed by Kernighan-Lin with joins:** ICP-GAEC-KLj-CC_t<TIME_LIMIT> ([source code](https://github.com/LPMP/))
* **(deterministic) Greedy Additive Edge Contraction, followed by Kernighan-Lin with joins:** GAEC-KLj-CC_t<TIME_LIMIT> ([source code](https://github.com/LPMP/))




# Data
The details about the generator are explained in Section 2.4 of *[Arınık'21]* (see Dataset 2.3 therein). All our results, as well as our generated signed networks with their optimal solutions, are publicly available on [FigShare](https://doi.org/10.6084/m9.figshare.14551113.v3) (*Chapter 2 >> Dataset 2.3*).




# Organization
Here are the folders composing the project:
* Folder `src`: contains the source code (R scripts).
* Folder `in`: contains the generated signed networks. 
* Folder `lib`: contains executable files related to the used external partitioning methods.
  * Folder `ExCC`: Executable file of the method `ExCC` whose the name will be `cplex-partition.jar`. See the *Installation* section for more details.
* Folder `out`: contains the folders and files produced by our scripts. See the *Use* section for more details.


# Installation
1. Install the [`R` language](https://www.r-project.org)
2. Install the [`Julia` language](https://julialang.org)
3. Install the following R packages (R is tested with the version 4.1):
   * [`igraph`](http://igraph.org/r/) Tested with the version 1.2.6.
   * [`XML`](https://cran.r-project.org/web/packages/XML/index.html)
   * [`expm`](todo)
   * [`ggplot2`](todo)
   * [`scales`](todo)
   * [`ggallin`](todo)
4. Install the following Julia packages (Julia is tested with the version 1.6.2):
   * [`DelimitedFiles`](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/)
   * [`DataStructures`](https://github.com/JuliaCollections/DataStructures.jl)
   * [`SparseArrays`](https://docs.julialang.org/en/v1/stdlib/SparseArrays/)
   * [`JLD`](https://github.com/JuliaIO/JLD.jl)
   * [`JuMP`](https://jump.dev/JuMP.jl/stable/) Tested with the version 0.21.9
   * [`CPLEX`](https://github.com/jump-dev/CPLEX.jl) Tested with the version 0.7.7
5. Install [`IBM CPlex`](https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en). Tested with the versions 12.8 and 20.1. Set correctly the variable `CPLEX.BIN.PATH` in `define-algos.R` (e.g. `/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/`).
   * For ubuntu, type the following command:
     * `sudo ./cplex_studio<YOUR_VERSION>.linux-x86-64.bin` 
       * The default installation location for education version is: `/opt/ibm/ILOG/CPLEX_Studio<YOUR_VERSION`.
       * The default installation location for trial version is: `/opt/ibm/ILOG/CPLEX_Studio_Community<YOUR_VERSION/cplex/bin/x86-64_linux/`.
6. Download the project of `ExCC` on [github](https://github.com/arinik9/ExCC). First, configure and then compile it. To test it, you can run the file `run.sh`. If everything works (i.e. if a file `sol0.txt` created in the output folder), move the executable file `ExCC.jar`, which is in `exe`, into the `lib/ExCC` folder in this project. `ExCC` will be used only for retrieving from CPLEX the ILP model of a perfectly structurally balanced signed graph. The set of Julia scripts rely on this exported ILP model.
7. Download the project of `EnumCC` on [github](https://github.com/arinik9/EnumCC). Move the executable files `EnumCC.jar` and `RNSCC.jar` into the `lib/EnumCC` folder in this project. To execute `CoNS(nbMaxEdit, pruningWithoutMVMO)`, we just need `RNSCC.jar`.


# Use
1. Set correctly the variable `CPLEX.BIN.PATH` in the file `src/define-algos.R`.
2. Set correctly the MATLAB installation path in `get.MLMSB.command()`, `get.BDCC.command()` and `get.ZONOCC.command()` in the file `src/define-algos.R` 
3. Open the `R` console.
4. Set the current directory as the working directory, using `setwd("<my directory>")`.
5. Run the main script `src/main.R`.


The script will produce the following subfolders in the folder `out`:
* `benchmark-analysis/partitions`: Folder containing all obtained partitions.
* `benchmark-analysis/csv`: Folder containing all csv results, as well as their corresponding plots.



# References
* **[Arınık'21]** N. Arınık, *Multiplicity in the Partitioning of Signed Graphs*. PhD thesis in Avignon Université (2021).
