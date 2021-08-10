# Experiments
This file is intended to serve as a guide for reproducing the results presented in our paper.

## Contents
* [Artifact Description](#artifact-description)
* [Run-time Environment](#run-time-environment)
  * [Hardware](#hardware)
  * [Software](#software) 
* [Initial Setup](#initial-setup)
  * [Cloning and Building](#cloning-and-building)
  * [Data sets](#data-sets)
* [Validating Setup](#validating-setup)
  * [Generating Smaller Data sets](#generating-smaller-data-sets)
  * [Running _ParsiMoNe_](#running-parsimone)
  * [Running _Lemon-Tree_](#running-lemon-tree) 
* [Measuring Performance](#measuring-performance)
  * [Sequential Performance](#sequential-performance)
  * [Parallel Performance](#parallel-performance) 

## Artifact Description
We have developed a parallel algorithm for learning module networks in parallel, based on the sequential _Lemon-Tree_ algorithm by [Bonnet et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983) In the experiments described below, we compare the performance of the original [_Lemon-Tree_](https://github.com/erbon7/lemon-tree) implementation with that of our software and also measure the parallel performance of our algorithm. We developed two main artifacts for this purpose:
1. [_ParsiMoNe_](https://github.com/asrivast28/ParsiMoNe/releases/tag/v1.0.0)  
This software implements the original sequential algorithm and our parallel algorithm for generating module networks.
3. [Modified _Lemon-Tree_](https://github.com/asrivast28/lemon-tree/tree/MatchOutput)  
In order to compare the performance of _ParsiMoNe_ with that of _Lemon-Tree_, we modified _Lemon-Tree_ to use the same PRNG
as _ParsiMoNe_ and made some other optimizations in _Lemon-Tree_ that we had implemented in _ParsiMoNe_.
These changes were made for the two implementation to produce the same output for the same input data set and parameters.  


## Run-time Environment
We provide a summary of the hardware and the software required for running the experiments. The details of the run-time environment that we used for our experiments (collected using the file [`collect_environment.sh`](https://github.com/SC-Tech-Program/Author-Kit/blob/master/collect_environment.sh)) are available in [`phoenix_environment.log`](phoenix_environment.log).

### Hardware
We used the [Phoenix cluster at Georgia Tech](https://docs.pace.gatech.edu/phoenix_cluster/gettingstarted_phnx/) for our experiments. Each node in the cluster has a 2.7 GHz 24-core Intel Xeon Gold 6226 processor and main memory of 192 GB or more. The nodes are connected via HDR100 (100 Gbps) InfiniBand. More details on the cluster resources can be found in the [cluster documentation](https://docs.pace.gatech.edu/phoenix_cluster/resources_phnx/).

### Software
We conducted our experiments on nodes running Linux [**RHEL** _v7.6_](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/7.6_release_notes/index) operating system. We used the following versions of the compiler and other libraries for experimenting with _ParsiMoNe_. 
* [**gcc** _v10.1.0_](https://gcc.gnu.org/gcc-10/changes.html)
* [**MVAPICH** _v2.3.3_](http://mvapich.cse.ohio-state.edu/static/media/mvapich/mvapich2-2.3.3-userguide.html)
* [**Boost** _v1.74.0_](https://www.boost.org/users/history/version_1_74_0.html)
* [**TRNG** _v4.22_](https://github.com/rabauke/trng4/releases/tag/v4.22)
* [**Armadillo** _v9.800.3_](http://sourceforge.net/projects/arma/files/armadillo-9.800.3.tar.xz)
* [**SCons** _v3.1.2_](https://scons.org/doc/3.1.2/HTML/scons-user.html)

The purpose of all the libraries is explained in more detail in [`README.md`](README.md#requirements). We also used [**OpenJDK** _v1.8.0\_262_](https://openjdk.java.net/install/) and the corresponding server VM for executing the original _Lemon-Tree_ implementation.

## Initial Setup
### Cloning and Building
#### _ParsiMoNe_
_ParsiMoNe_ can be downloaded by cloning this Github repo, along with all the submodules, by executing the following:
<pre><code>git clone --recurse-submodules git@github.com:asrivast28/ParsiMoNe.git
</code></pre>
Once all the requirements have been installed, the executable for measuring performance can be built by executing the following:
<pre><code>scons TIMER=1
</code></pre>
More information on different build options can be found in [`README.md`](README.md#building)

#### Lemon-Tree
The modified _Lemon-Tree_ can be downloaded and built by executing the following:
<pre><code>git clone git@github.com:asrivast28/lemon-tree.git
cd lemon-tree/LemonTree
git checkout -b MatchOutput origin/MatchOutput
ant jar
</code></pre>

The modified _Lemon-Tree_ uses the random number generators from the TRNG library, made available in _Java_ using [_Java Native Interface_](https://docs.oracle.com/javase/8/docs/technotes/guides/jni/). The corresponding library needs to be built separately by executing the following:
<pre><code>cd lemon-tree/LemonTree/src/lemontree/utils
make
</code></pre>

### Data sets
We used the following two gene expression data sets from two model organisms for our experiments: _Saccharomyces cerevisiae_, or Baker's yeast, and the _Arabidopsis thaliana_ plant.
* A data set created by [Tchourine et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5987223/) from multiple RNA-seq studies of _S. cerevisiae_.
The data set contains 2,577 observations for 5,716 genes and it can be downloaded from [Zenodo](https://zenodo.org/record/3355524#.Xpx0t1NKhhE).
* A manually curated data set created from multiple microarray studies of _A. thaliana_ focusing only on the studies of the development process in the plant. 
The data set contains 5,102 observations for 18,373 genes and it can also be downloaded from [Zenodo](https://zenodo.org/record/4672797#.YG9TQhNKhQI).

## Validating Setup
The experimental setup can be validated using smaller data sets as described below.

### Generating Smaller Data sets
Smaller data sets from the complete yeast or _A. thaliana_ data set can be generated for validation purpose, using common Linux command line utilities. For example, a data set with the first 100 observations for the first 100 variables in the yeast data set can be obtained by executing the following:
<pre><code>n=100;m=100;head -$(($n+1)) yeast_microarray_expression.tsv | cut -d $'\t' -f 1-$(($m+1)) > yeast_n${n}_m${m}.tsv
</code></pre>
The values of `n` and `m` in the command can be varied to generate multiple smaller data sets.

### Running _ParsiMoNe_
Using the above generated smaller data set, module network can be generated sequentially using _ParsiMoNe_ by executing the following:
<pre><code>./parsimone -f yeast_n100_m100.tsv -n 100 -m 100 -c -v -i -s $'\t' -g experiment_configs_seed0.json -o sequential_parsimone
</code></pre>
This will create a directory called `sequential_parsimone` with all the files. Similarly, it can be executed in parallel for generating module networks as follows:
<pre><code>mpirun -np 4 ./parsimone -f yeast_n100_m100.tsv -n 100 -m 100 -c -v -i -s $'\t' -g experiment_configs_seed0.json -o parallel_parsimone
</code></pre>
The generated networks are expected to be the same, irrespective of the number of processors used. This can be verified using the script [`compare_lemontree.py`](https://github.com/asrivast28/bn-utils/blob/main/scripts/compare_lemontree.py) that we created for the purpose, by executing the following:
<pre><code>common/scripts/compare_lemontree.py sequential_parsimone parallel_parsimone
</code></pre>

### Running _Lemon-Tree_
The original _Lemon-Tree_, built as described above, can be used for generating module networks using the script [`parsimone_lemontree.py`](https://github.com/asrivast28/bn-utils/blob/main/scripts/parsimone_lemontree.py). This script accepts the same arguments as _ParsiMoNe_ and can be executed as:
<pre><code>common/scripts/parsimone_lemontree.py -f yeast_n100_m100.tsv -n 100 -m 100 -c -v -i -s $'\t' -g experiment_configs_seed0.json -o lemontree_parsimone
</code></pre>
Again, given the same input data set and parameters, we expect _Lemon-Tree_ to generate the same network as _ParsiMoNe_. This can be verified as:
<pre><code>common/scripts/compare_lemontree.py sequential_parsimone lemontree_parsimone
</code></pre>

## Measuring Performance
We provide a Python script, [`parsimone_experiments.py`](https://github.com/asrivast28/bn-utils/blob/main/scripts/parsimone_experiments.py), for easily
experimenting with _ParsiMoNe_ as well as _Lemon-Tree_ and measuring their performance. 
The commands using the script below expect that the script is executed from the `ParsiMoNe` directory cloned above and the two data sets are available at the following paths in the directory: `data/yeast/yeast_microarray_expression.tsv` and `data/athaliana/athaliana_development_exp.tsv`

### Sequential Performance
We compared the sequential performance of _ParsiMoNe_ with that of _Lemon-Tree_ for learning module networks.  
We obtained the run-times of our implementation for 15 different subsamples of the yeast data set using three different random seeds by executing the following:
<pre><code>for seed in {0,1,2}; do
  common/scripts/parsimone_experiments.py -r 1 -p 1 -d yeast -n 1000 2000 3000 -m 125 250 500 750 1000 -g "\-g experiment_configs_seed${seed}.json" --results ours_yeast_sequential_seed${seed}.csv -b . -s . --output-suffix _seed${seed}
done
</code></pre>
The above run is expected to take about six days and will generate three CSV files with the run-times of different components in our implementation: `ours_yeast_sequential_seed0.csv`, `ours_yeast_sequential_seed1.csv`, and `ours_yeast_sequential_seed2.csv`.

Using the script `parsimone_lemontree.py` described earlier, the sequential performance of _Lemon-Tree_ for the 15 subsampled data sets can be measured similar to that of our implementation by executing the following:
<pre><code>for seed in {0,1,2}; do
  common/scripts/parsimone_experiments.py -r 1 -p 1 -d yeast -n 1000 2000 3000 -m 125 250 500 750 1000 -g "\-g experiment_configs_seed${seed}.json" --results lemontree_yeast_sequential_seed${seed}.csv -b . -s . --output-suffix _seed${seed} --lemontree
done
</code></pre>
Again, this will generate three CSV files with the run-times of different _Lemon-Tree_ components: `lemontree_yeast_sequential_seed0.csv`, `lemontree_yeast_sequential_seed1.csv`, and `lemontree_yeast_sequential_seed2.csv`. This run is expected to take about 21 days.

Then, the outputs generated by _Lemon-Tree_ and our implementation in the above runs for different subsampled data sets and different PRNG seeds can be compared, using `compare_lemontree.py`, by executing the following:
<pre><code>for n in {1000,2000,3000}; do
  for m in {125,250,500,750,1000}; do
    for seed in {0,1,2}; do
      common/scripts/compare_lemontree.py yeast_n${n}_m${m}_lemontree_seed${seed} yeast_n${n}_m${m}_seed${seed}
    done
  done
done
</code></pre>

### Parallel Performance
#### Smaller Data sets
First, we measured the strong scaling parallel performance of _ParsiMoNe_ for learning the network for all the variables in the yeast data set using a subset of observations in the data set.
The following can be executed for the purpose:
<pre><code>for seed in {0,1,2}; do
  common/scripts/parsimone_experiments.py -r 1 --ppn 24 -p 1 2 4 8 16 32 64 128 256 512 1024 -d yeast -n 5716 -m 125 250 500 750 1000 -g "\-g experiment_configs_seed${seed}.json" --results parallel_yeast_small_seed${seed}.csv -b . -s . --output-suffix _seed${seed}
done
</code></pre>
Similar to the sequential performance, this will generate seed-specific CSV files (`parallel_yeast_small_seed0.csv`, etc.) with run-times of _ParsiMoNe_ when using different number of cores. All these runs are expected to take about 24 days.

This command automatically compares the output generated when using different number of processors with the first output generated for every combination of `n` and `m`.
Therefore, in this case, it compares against the output generated by `p=1` and errors out in case of any mismatches.

### Big Data sets
Then, we measured the parallel run-times for learning module networks from the two complete data sets.  
For the yeast data set, we conducted strong scaling exxperiments by executing the following:
<pre><code>for seed in {0,1,2}; do
  common/scripts/parsimone_experiments.py -r 1 --ppn 24 -p 4 8 16 32 64 128 256 512 1024 2048 4096 -d yeast -g "\-g experiment_configs_seed${seed}.json" --results parallel_yeast_complete_seed${seed}.csv -b . -s . --output-suffix _seed${seed}
done
</code></pre>
This run is expected to take approximately 26 days.

Since learning network from the _A. thaliana_ development data set requires a lot of time, we experimented with the data set only on larger number of cores as follows:
<pre><code>for seed in {0,1,2}; do
  common/scripts/parsimone_experiments.py -r 1 --ppn 24 -p 1024 2048 4096 -d development -g "\-g experiment_configs_seed${seed}.json" --results parallel_development_complete_seed${seed}.csv -b . -s . --output-suffix _seed${seed}
done
</code></pre>
This run is expected to take approximately 13 days.
