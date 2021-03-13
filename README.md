# AraSim
ARA simulation and monte carlo

## Attribution

If you use AraSim in your work, please cite:

> P. Allison et al. for the ARA Collaboration. "First Constraints on the Ultra-High Energy Neutrino Flux from a Prototype Station of the Askaryan Radio Array." 
> Astroparticle Physics, Vol 70, 2015, Pg 62-80, https://doi.org/10.1016/j.astropartphys.2015.04.006.

And please add the following sentence to your acknowledgement section:

> We are grateful to the ARA Collaboration for making available the AraSim simulation program used in this work.

And provide the link to this repository (https://github.com/ara-software/AraSim) in your paper, e.g. as a footnote or reference.

## Prerequisites
* [ROOT](https://root.cern/install/) -- ROOT 6 is preferred
* [Boost](https://www.boost.org/users/download/) -- Greater than 1.71 required; necessary for ray tracing and interpolation
* [FFTW](http://www.fftw.org/download.html) -- required for FFTs
* libRootFftwWrapper -- a ROOT wrapper for FFTW 3 downloadable from [Ryan Nichol's GitHub](https://github.com/nichol77/libRootFftwWrapper), with documentation [here](http://www.hep.ucl.ac.uk/uhen/libRootFftwWrapper/)

## Installation
0. Install all dependencies, and ensure the libraries and headers are appropariately included in e.g. your `$LD_LIBRARY_PATH`.

1. Checkout the code from the ARA Git repository, eg.: `git clone https://github.com/ara-software/AraSim.git`

2. Run `make`. 
