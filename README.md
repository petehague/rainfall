Rainfall
========

Rainfall is a platform for Bayesian analysis, currently only supporting MCMC.

Usage
-----

Look at testfunc.cpp - this is the likelihood evaluation function. Make a copy of this and edit it to match your own problem. Create an extra line in the makefile along the lines of:

rainfall-function: $(COREOBJS) bin/function.o
	$(CXX) $(CPPFLAGS) $(COREOBJS) bin/function.o -orainfall-function

and then type make rainfall-function. The resultant executable will perform a simple MCMC analysis of your function
