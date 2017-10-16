# fastq-shuffle.pl

If you are using this script, please cite [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1012747.svg)](https://doi.org/10.5281/zenodo.1012747).

A small program to convert and combine bismark coverage output files
into DSS input files.

# SYNOPSIS

    generate_dss_input bismark1.cov.gz bismark2.cov

# OUTPUT

The output is a csv file containing a header and the columns

- 1. chr
- 2. pos
- 3. N
- 4. X

as required by DSS. The counts for total coverage and methylated
coverage of all input files is combined and the sum will be printed to
stdout.

# SEE ALSO

- Bismark

    A mapper for bs-seq `https://www.bioinformatics.babraham.ac.uk/projects/bismark/`

- DSS

    A R package for differential methylation experiments `https://www.bioconductor.org/packages/release/bioc/html/DSS.html`

# AUTHOR

This script was written by Frank Förster `frank.foerster@ime.fraunhofer.de`.

# CHANGELOG

- 2017-10-16 v0.1

    First working version.

- 2017-10-16 v0.1.1

    New version using array instead of hashes to reduce memory
    consumption.

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
