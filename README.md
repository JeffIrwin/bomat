
![](https://github.com/JeffIrwin/bomat/workflows/CI/badge.svg)

# bomat

*Cultural learnings of Bohemia for make benefit glorious library of LAPACK*

A set of Fortran functions for visualizing Bohemian matrices 

![](https://raw.githubusercontent.com/JeffIrwin/bomat/main/examples/bomat-7c--2022-02-01.png)

See the original project:  https://github.com/BohemianMatrices/BHIME-Project

## Download

    git clone --recursive https://github.com/JeffIrwin/bomat
    cd bomat

## Build

Use the provided CMake wrapper script.  If using `gfortran` and reference LAPACK, build that first one time too:

    ./build-lapack.sh
    ./build.sh

If using `ifort`, there is no need to build LAPACK.  Just export the environment variable to tell bomat's CMake list to use the Intel Fortran compiler:

    export BOMAT_INTEL=true
    ./build.sh

## Run

The `bomat` program loads its input from a JSON file:

    ./build/bomat examples/bomat-7c--2022-02-01.json

where `examples/bomat-7c--2022-02-01.json` configures a few input options:

    {
    	# This is a comment (non-standard JSON extension)
        
    	# Complex numbers
    	"Population":
    	[
    		# Re       , Im
    		 0         ,  1,
    		-0.8660254 , -0.5,
    		 0.8660254 , -0.5,
    		 0.1       ,  0.1,
    		 0.1       , -0.1,
    		-0.1       ,  0.1,
    		-0.1       , -0.1
    	],
    
    	# Integers-only.  Zeros will remain 0, non-zeros will be sampled randomly
    	# from population.  This JSON array is rank-1, but it is reshaped in
    	# a row-major sense into a rank-2 matrix
    	"Template matrix":
    	[
    		1, 1, 0, 0, 0, 0, 0, 0,
    		1, 1, 1, 0, 0, 0, 0, 0,
    		0, 1, 1, 1, 0, 0, 0, 0,
    		0, 0, 1, 1, 1, 0, 0, 0,
    		0, 0, 0, 1, 1, 1, 0, 0,
    		0, 0, 0, 0, 1, 1, 1, 0,
    		0, 0, 0, 0, 0, 1, 1, 1,
    		0, 0, 0, 0, 0, 0, 1, 1
    	],
    
    	"Samples": 100000000,
    	"Image size": 3840,
    
    	"Colormap file": "submodules/colormapper/submodules/colormaps/ColorMaps5.6.0.json",
    	"Colormap name": "Magma (matplotlib)"
    }

The path to the `Colormap file`, if not absolute, must be relative to the runtime directory.

### Command-line arguments

Use `-h` for help:

     Starting bomat
    
    Usage: bomat [-h] [-p] [-e] FILE.JSON
    
    Calculate Bohemian matrix eigenvalues and export a plot to a PNG file
    
    Positional arguments:
    FILE.JSON   Configuration filename for setting inputs
    
    Optional arguments:
    -h, --help  Show this help message and exit
    -e          Calculate and export eigenvalues without plotting
    -p          Plot eigenvalues from previous job

    Sample FILE.JSON contents are like this:
    {
    [omitted from README]
    }

I recommend initially running with a low number of `Samples`, around 1 million,
while playing with the `Population` and `Template matrix` structure.  Once you
find something interesting, use a higher number of `Samples` for
a higher-quality image.

Finally, you can use `-p` to re-plot with a different colormap without
recalculating eigenvalues from scratch.
 
