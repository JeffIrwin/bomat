
# bomat

*Cultural learnings of Bohemia for make benefit glorious library of LAPACK*

A set of Fortran functions for visualizing Bohemian matrices 

![](https://raw.githubusercontent.com/JeffIrwin/bomat/main/examples/bomat-7c--2022-02-01.png)

See the original project:  https://github.com/BohemianMatrices/BHIME-Project

## Download

    git clone --recursive https://github.com/JeffIrwin/bomat
    cd bomat

Use CMake, or run the provided CMake wrapper script:

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

