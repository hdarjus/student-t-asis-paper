Code for reproducing the simulation study. Written by Darjus Hosszejni, 2021.

The computations were made on a cluster of computers, operated by the Sun Grid Engine (SGE) version 8.1.9. The following line of code is needed for the reproduction of the Markov chains:

qsub cluster-job-sge-simulation.qsub

This executes the R code in 'simulate_student_asis.R' 7425 times, each time providing the envinronment variable SGE_TASK_ID with a unique value between 1 and 7425. The R code figures out which setting belongs to that task id and produces four Markov chains.
This behavior can be reproduced without the SGE; however, this may take days on a regular personal computer.
At the end, the results are collected in a folder called 'results-simulation/results3'.

Finally, by processing the output from the cluster computation, the R code in 'evaluate-results-simulation.R' reproduces the figures and tables in the paper.
