# Lung Model #

## Installing

Build it with

`make`

The executable will be installed in this working directory.

## Running
Usage: lungmodel <dim_x> <dim_y> <dim_z> <offset_x> <offset_y> <offset_z>

To run the full lung model, execute:

`./lungmodel 25255 21031 43734 0 0 0`

It will create an 'airway.csv' whose full filepath is required in the simcov config file. For example:

```
; Directory containing files for lung model
  lung-model =                  /users/projects/simcov/lungmodel

```
