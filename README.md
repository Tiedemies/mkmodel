# mkmodel
This package contains my simulator and statistic analyzer for trading information research.

Filestructure is complex and in dire need of documentation. 

$ make test_io 
will make the binary test_io, which will produce the regressor matrices for all transactions. It also does a check of data
integrity 

$ make RunEmpirical
will make the RunEmpirical binary which will try to fit the two-parameter model. This is currently not functional
due to overhaul in the I/O

$ make vaccinator
will make the single-run "vaccinator" program that is intended (eventually) to fit the multiparameter model and 
solve the so-called "vaccination problem", i.e., find out which nodes or connections need to be shut down in order to
minimize the spread of information as measured by the number of expected infections


