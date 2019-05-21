# CS 470 Research - NBody Problem
This problem was researched in Spring 2019 by Ben Bole, John Latino, Richard Bimmer, and Kevin Kelly under Dr. Mike Lam. We wished to compare parallel methods of various approaches to the N-Bodies Simulation Problem. In this repository we contain code for all of the algorithms which we compared as well as a directory for testing those algorithms.

A more detailed summary of the selected algorithms and approaches, our experiment methodology, results, and future reccomendations can be found [here](https://drive.google.com/open?id=1pRMRGz9M3wwcX9fkrnGJekzg1Cn7ktZB) (Please cite this as a resource, as well as our referenced works if you plan to continue research with these applications)

# Repository Contents

### Barnes-Hut Similation
This directory contains the Barnes-Hut simulation code, which is written in C++. The code was originally implemented by @barkm on Github and the original repository can be found [here](https://github.com/barkm/n-body). We manipulated the code so that we could time it. Further details on this repository can be found within the folder it is currently in.
### Fast-Multipole Method
This directory contains an implementation of the Fast Multipole Method, written in C. This code was originally implemented by @cfkane24 on Github and the original repository can be found [here](https://github.com/cfkane24/fast-multipole-method). Further details on this repository can be found within the folder it is currently in.
### Parker-Sochaki Method
This directory contains an implementation of the Parker-Sochaki Method, written in Fortran90. This code was originally implemented by Dr. C. David Pruett, Dr. William H. Ingham, and Dr. Ralph D. Herman and can be found written about [here](http://educ.jmu.edu/~sochacjs/PruettInghamHerman.pdf). Further details on this repository and tips on understanding Fortran code can be found within the folder it is currently in.
### Input Generation
This directory contains all scripts for generating input for each algorithm listed above. The folders "input_files_log_10", "input_files_pow", and "input_files_swarm" all contain the input for various test suites for each algorithm.
