
## Using a Whole Genome Co-expression Network to Inform the Functional Characterisation of Predicted Genomic Elements from Mycobacterium tuberculosis Transcriptomic Data

Jennifer J. Stiens, Yen Yi Tan, Rosanna Joyce, Kristine Arnvig, Sharon Kendall, Irilenia Nobeli

The data can be explored using a R ShinyApp. Users can choose a module, transcript or set of genomic coordinates to explore the module associations, associated predicted non-coding RNAs and expression profiles of all of the transcripts analysed in the paper.

The ShinyApp can be run in R directly from the GitHub server.

The following packages must be installed in the R environment:

```
install.pakcages("here")
install.packages("shiny")
install.packages("shinyjs")
install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("DT")
install.packages("JBrowseR")

```
Run the following code:

```
library(shiny)

runGitHub("mtb_wgcna", "jenjane118")

```

Open in browser (Chrome, preferably) to get optimal functionality.

### Abstract

A whole genome co-expression network was created using Mycobacterium tuberculosis transcriptomic data from publicly available RNA-sequencing experiments covering a wide variety of experimental conditions. The network includes expressed regions with no formal annotation, including putative sRNAs and UTRs, along with the protein-coding genes. Non-coding RNA were among the most well-connected members of the module sub-networks, making up more than half of the ‘hub’ genes in modules that include protein-coding genes known to be part of regulatory systems involved in stress response and host adaptation. This dataset provides a valuable resource for investigating the role of non-coding RNA in transcriptomic remodelling. Based on their connections to genes with known functional groupings and correlations with replicated host conditions, predicted non-coding RNA elements can be screened as suitable candidates for further experimental validation.


### The summary and overview of scripts relevant for paper are found in Mtb_module_overview.Rmd

Including:

1) RNA-seq processing and mapping

2) Transcript prediction

3) Feature quantification

4) Normalisation, transformation and batch correction

5) Network creation with WGCNA

6) Determining categories and positions of ncRNA transcripts

7) Further analysis including functional enrichment
