# MiCoP

### Basic usage

First, use a mapping tool that produces sam output, such as BWA.

Then, run the compute-abundances.py tool as such:  
python compute-abundances.py alignments.sam [--virus OR --fungi] [options]

To see options, run "python compute-abundances.py -h"

### Simulated data

The data used for our simulations can be found in the simulated\_data directory.

The mock community datasets are taken from the following papers:

Viral data: Conceicao-Neto, N., Zeller, M., Lefrere, H., De Bruyn, P., Beller, L., Deboutte, W., ... & Matthijnssens, J (2015). Modular approach to customise sample preparation procedures for viral metagenomics: a reproducible protocol for virome analysis. *Scientific reports*, 5, e16532. doi:10.1038/srep16532

Fungi data: Tonge, D. P., Pashley, C. H., & Gant, T. W. (2014). Amplicon–Based Metagenomic Analysis of Mixed Fungal Samples Using Proton Release Amplicon Sequencing. *PloS one*, 9(4), e93849. https://doi.org/10.1371/journal.pone.0093849

Finally, the HMP data was downloaded using the instructions in the MetaPhlAn tutorial: https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial

