# Genetic algorithm for the search of cancer subtypes with clinical significance according to their gene expression patterns

**Motivation**: Clustering analysis has been long used to find underlying structures in different omics data such as gene expression profiles. This data typically presents high number of dimensions and has been used successfully to find co-expressed genes in samples that share similar molecular and clinical characteristics. Nevertheless, the clustering results are highly dependent of the features used and the number of clusters considered, while the partition obtained does not guarantee clinically relevant findings.

**Methods**: We propose a multi-objective optimization algorithm for disease subtype discovery based on a non-dominated sorting genetic algorithm. Our proposed framework combines the advantages of clustering algorithms for grouping heterogeneous omics data and the searching properties of genetic algorithms for feature selection and optimal number of clusters determination to find features that maximize the survival difference between subtypes while keeping cluster consistency high.

### To run the algorithm all the data must be in the *Data* folder

* In the *Functions* folder are all the files to run the algorithm.
* First, be sure you have all the needed libraries (listed in the 'libraries.R' file) installed and properly configured.
* To run the code, open the R terminal and set the *Genetic_algo* folder as the main path and then source the 'RUN.R' file.
* To edit the hyperparamethers, change the repective values of the parameters in the 'Hyperparameters.R' file.
