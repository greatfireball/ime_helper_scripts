# cluster\_build.pl

If you are using this script, please cite [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1012747.svg)](https://doi.org/10.5281/zenodo.1012747).

A small program to generate miRNA clusters for our database based on a
table of miRNA positions.

# SYNOPSIS

    # build a list of clusters with max distance 7500 nt from input.csv
    cluster_build.pl input.csv 7500

# PARAMETER

- Input file

    First parameter is the input csv file.

- cluster size

    Second parameter is the maximum distance between two miRNAs still
    belonging to the same cluster. If no parameter is provided, 6000 is
    assumend. If multiple clustersizes should be testes, just put them
    comma/space separated as second parameter.

# INPUT

A csv file with comma as field separator containing the following
information:

- position\_entry\_id \[mandatory, unique\]
- genome\_id
- start\_position \[mandatory\]
- stop\_position \[mandatory\]
- chromosome \[mandatory\]
- strand \[mandatory\]
- miRNA \[mandatory\]

Quotes will be stripped out.

# OUTPUT

Output will be printed to STDOUT and additional information are
printed to STDERR. The csv content of STDOUT contains the following
columns:

- cluster\_id
- position\_entry\_id
- maximal\_distance

The columns are comma separated.

# AUTHOR

This script was written by Frank Förster `frank.foerster@ime.fraunhofer.de`.

# CHANGELOG

- 2017-10-18 v0.1.2

    First working version.

- 2017-10-18 v0.1.3

    Updated documentation/allows multiple cluster sizes

# COPYRIGHT AND LICENCE

MIT License

Copyright (c) 2017 Frank Förster

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
