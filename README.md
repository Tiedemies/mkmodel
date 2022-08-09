# mkmodel
This package contains my simulator and statistic analyzer for trading information research.

Filestructure is complex and in dire need of documentation, but now some architectural sense has been made.
The build now uses cmake.

$ cmake .
$ make 

Binaries that are generated are:
$ RunEmpirical
This will try to fit the data into a rudimetary version of the hidden cascade. Currently it is defunct. The optimizer is retained only for reference.

$ Vaccinator
This is work in progress which uses the IndustryCascade - wrapper for the HiddenCascade class. It supports simulations where the number of simulations is proportional to the number of announcements by the company. 

$ test_io
This binary will generate the data matrices used in our seminal paper based on the data files that are (at this moment) hard coded in the file
util/defs.hpp 
If you wish to generate the matrices with different data, you must change defs.hpp for the appropriate constants. TODO-list includes creating an installation script that only requires the user to specify the root of the data directory. 



