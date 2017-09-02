# Smatch (semantic match) tool

This is source code of [smatch](http://amr.isi.edu/evaluation.html), an evaluation tool for AMR (Abstract Meaning Representation). 

The code was mostly developed during 2012-2013, and has undergone many fixes and updates. It is now hosted on github for better collaboration.

More details and updates about AMR and smatch can be found in USC/ISI's AMR site: http://amr.isi.edu/index.html


---
## Smatch ILP
Smatch has an implementation using 0-1 Integer Linear Programming technique
The `smatch_ilp.py` provides this ILP based algorithm which uses [Gurobi](http://www.gurobi.com/) ILP solver


    # Install miniconda (if not already installed)
    # create python env for python2.7
    # source  activate python2.7
    conda config --add channels http://conda.anaconda.org/gurobi

    conda install gurobi

### Usage:

    $ python smatch_ilp.py -h
        usage: smatch_ilp.py [-h] [-v] [-vv] amrfile amrfile

        Smatch ILP

        positional arguments:
          amrfile     Path to file having AMR

        optional arguments:
          -h, --help  show this help message and exit
          -v          Verbose (log level = INFO)
          -vv         Verbose (log level = DEBUG)

#### Examples

    python smatch_ilp.py test_input{1,2}.txt                    # minimal output
    python smatch_ilp.py test_input1.txt test_input2.txt -v     # Output includes variable and triple matchings
    python smatch_ilp.py test_input1.txt test_input2.txt -vv    # Output includes ILP optimization logs