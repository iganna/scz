# sem_config: get configuration of a gene pathway

## Description

This repository contains the pipeline to select and estimate a SEM model based on the pathway with ambiguous nodes. Topologies of pathways are incorporated in `pathways.R`.  

To run the pipeline, one can use two options:
* Notebook `pipeline.ipynb`
* R script `pipeline.R`

## Data

The table with gene expression data is in `dExpr.txt`

## Requirements

* R package `lavaan`   
* R package `DiagrammeR`, **version 8.2**. The folder with this version is cloned to the current repository.

## References

Anna. A. Igolkina , Chris Armoskus, Jeremy. R. B. Newman, Oleg. V. Evgrafov, 
Lauren. M. McIntyre, Sergey. V. Nuzhdin, Maria. G. Samsonova, *Analysis of gene expression variance in schizophrenia using structural equation modeling* , 2018, Frontiers in Molecular Neuroscience.


## Authors

**Anna Igolkina** [e-mail](mailto:igolkinaanna11@gmail.com)

## License information

This repository contains open-sourced pipeline licensed under the [MIT license](https://opensource.org/licenses/MIT).
