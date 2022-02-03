#!/bin/bash

# Compile with ifort, with colormapper routines removed:
ifort main.f90 -qmkl -qopenmp

