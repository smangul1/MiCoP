# MiCoP

### Setup

Simply run "./setup.sh".

This will download about 8.4 GB of data, which will be uncompressed into files totalling about 13 GB of disk space.

### Basic usage

First, use a mapping tool that produces sam output, such as BWA. For convenience, we include a copy of BWA kit in our repository (see License Info). The setup script will also download pre-indexed NCBI databases for viruses and fungi. For info on how to use BWA, see: http://bio-bwa.sourceforge.net/bwa.shtml#13

Then, run the compute-abundances.py tool as such:  
python compute-abundances.py alignments.sam [--virus OR --fungi] [options]

To see options, run "python compute-abundances.py -h"

### Simulated data

The data used for our simulations can be found in the simulated\_data directory.

The mock community datasets are taken from the following papers:

Viral data: Conceicao-Neto, N., Zeller, M., Lefrere, H., De Bruyn, P., Beller, L., Deboutte, W., ... & Matthijnssens, J (2015). Modular approach to customise sample preparation procedures for viral metagenomics: a reproducible protocol for virome analysis. *Scientific reports*, 5, e16532. doi:10.1038/srep16532

Fungi data: Tonge, D. P., Pashley, C. H., & Gant, T. W. (2014). Amplicon–Based Metagenomic Analysis of Mixed Fungal Samples Using Proton Release Amplicon Sequencing. *PloS one*, 9(4), e93849. https://doi.org/10.1371/journal.pone.0093849

Finally, the HMP data was downloaded using the instructions in the MetaPhlAn tutorial: https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial

### License Info

For convenience, we include a version of the BWA kit in our repository. This is allowed under BWA's GPLv3 license, under the requirement that we also adopt GPLv3. We also use a code snippet in the setup script that is adapted from a Stack Overflow answer, permitted for inclusion under the Creative Commons license.
#
