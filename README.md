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


## INSTALL MongoDB C++ drivers

Source: http://mongocxx.org/mongocxx-v3/installation/

Step 1:

$sudo apt install cmake libssl-dev libsasl2-dev

Step 2 (install mongo-c-driver):

source: http://mongoc.org/libmongoc/current/installing.html

$ wget https://github.com/mongodb/mongo-c-driver/releases/download/1.17.4/mongo-c-driver-1.17.4.tar.gz
$ tar xzf mongo-c-driver-1.17.4.tar.gz
$ cd mongo-c-driver-1.17.4
$ mkdir cmake-build
$ cd cmake-build
$ cmake -DENABLE_AUTOMATIC_INIT_AND_CLEANUP=OFF ..
$ make
$ sudo make install

step 3 (install mongo-cxx drivers in /usr/local)

source: http://mongocxx.org/mongocxx-v3/installation/linux/

$curl -OL https://github.com/mongodb/mongo-cxx-driver/releases/download/r3.6.2/mongo-cxx-driver-r3.6.2.tar.gz
$tar -xzf mongo-cxx-driver-r3.6.2.tar.gz
$cd mongo-cxx-driver-r3.6.2/build

$ cmake ..                                \
$     -DCMAKE_BUILD_TYPE=Release          \
$     -DCMAKE_INSTALL_PREFIX=/usr/local

$ sudo cmake --build . --target EP_mnmlstc_core
$ cmake --build .
$ sudo cmake --build . --target install

Note: If when running the model with memore, there is an error like this one 
(error while loading shared libraries: libmongocxx.so.\_noabi: cannot open shared object library file: No such file or directory)
try to upgrade the system:
$ sudo apt upgrade (this updates g++ and c++ and that is the problem)