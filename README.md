# ABEILLE
a novel method for ABerrant Expression Identification empLoying machine LEarning from sequencing data

## Main advantage of ABEILLE

ABEILLE is fast and new method to identify aberrant expression in count table using only machine learning. Since ABEILLE is build to take any kind of normalized count, every function is flexible to adapt to the dataset. As ABEILLE is based on machine learning, you can feed your dataset with healthy samples or other dataset, but it is important to normalize (batch correction) everything before using the model.

### Prerequisites

We highly recommand to read the vignette explaining how to use every function.

### Installing

Build on the whole architecture or load the tar.gz file as a package.

From Rstudio:

```
setwd("/my/path/to/ABEILLE")
devtools::build()
devtools::install()
library(ABEILLE)
```


