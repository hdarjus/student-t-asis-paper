Code for reproducing the Fisher information plot. Written by Darjus Hosszejni, 2021.

The computations were made on a cluster of computers, operated by the Sun Grid Engine (SGE) version 8.1.9. The following line of code is needed for the reproduction of the Markov chains:

`qsub cluster-job-sge-fisher.qsub`

This executes the R code in 'fisher-info.R' 10100 times, each time providing the envinronment variable SGE_TASK_ID with a unique value between 1 and 10100. The R code figures out which setting belongs to that task id and evaluates the difference of the Fisher informations $I_{\tau}$ and $I_u$ (using the notation of the paper) using Monte Carlo estimation.
This behavior can be reproduced without the SGE; however, this may take hours on a regular personal computer.
At the end, the results are collected in a folder called 'results-fisher/results-fisher2'.

Finally, by processing the output from the cluster computation, the R code in 'evaluate-results-fisher.R' reproduces the figure in 'efficiency-contour.pdf' and some further insightful figures.
