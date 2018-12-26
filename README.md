# AMMP-PDEVS
Automatic Modeling of Metabolic Pathways in Parallel DEVS

## Dependencies
 1. C++17 >
 2. Boost 1.57 >
 3. TinyXML-2
 4. Cadmium
 5. DESTimes
 6. DEVSDiagrammer

*Note:* Dependencies from 3 to 6 comes in this project as submodules: git submodule update -i --recursive

## Standard names
Some variables and concepts have standarized names that are used across the project code, the list of this names is:
 * cid: Compartment ID
 * sid: Specie ID
 * eid: Enzyme ID
 * rid: Reaction ID
 * rsn: Reaction Set Name

## Hot to install the model generator Python module
 1. cd model_generator
 2. python setup install

## How to generate models from SBML files (once installed the model generator)
 1. pmgbp_generate_model -f <sbml_file_path>
*Note:* For a complete list of the model generator parameters do: pmgbp_generate_model --help

## How to compile a generated model
 1. Having the model in the project root dir (where the main.cpp file is palced) run: make all

## Notes:
*The directory structure format:* The structure used in this project for the directory structure was taken from
https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/ 
