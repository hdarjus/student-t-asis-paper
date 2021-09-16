Description of the Supplementary Material for the attached paper. Written by Darjus Hosszejni,
2021.

1. LICENSE: GNU General Public Lincese version 3 or higher. You must have received a copy of the license as part of this material.
2. The implementations of algorithms AA, SA, and ASIS can be found in the source R package 'sample.student.asis', version 0.0.1, which is part of the supplementary material as the archive file 'sample.student.asis_0.0.1.tar.gz'. This file can be opened with free and open source software.
To install this package, first install R and the 'progress' package, then this package. In the R console, execute:
```r
install.packages("progress")
install.packages("sample.student.asis_0.0.1.tar.gz", method = "source", repos = NULL)
```
3. The contour plots of the conditional densities in the paper have been generated using the R code in 'generate-fig.R'.
4. The computations behind the Fisher information contour plot can be found in folder 'Fisher-information'. Please see the README file there for more information.
5. The computations behind the simulation study can be found in folder 'Simulation-study'. Please see the README file there for more information.
6. Finally, the computations behind the application can be found in folder 'Application'. Please see the README file there for more information.
