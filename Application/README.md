Code for reproducing the application. Written by Darjus Hosszejni, 2021.

The computations were made on a cluster of computers, operated by the Sun Grid Engine (SGE) version 8.1.9. The following line of code is needed for the reproduction of the Markov chains:

qsub cluster-job-sge-application.qsub

This executes the R code in 'geweke-application.R' 42 times, each time providing the envinronment variable SGE_TASK_ID with a unique value between 1 and 42. The R code figures out which setting belongs to that task id and produces the Markov chain for the given input variable and parameterization.
This behavior can be reproduced without the SGE; however, this may take hours on a regular personal computer.
At the end, the results are collected in a folder called 'results-application/results-app5'.

Finally, by processing the output from the cluster computation, the R code in 'evaluate-results-application.R' reproduces the figures and tables in the paper.
