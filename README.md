# AMMP-PDEVS
Automatic Modeling of Metabolic Pathways in Parallel DEVS

## Dependencies
 1. C++14 >
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

## Notes:

*The directory structure format:* The structure used in this project for the directory structure was taken from
https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/ 
