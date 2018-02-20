# MiCoP

MiCoP (Microbial Community Profiling) is a method for high-accuracy profiling of viral and fungal metagenomic communities. MiCoP uses a mapping-based approach that produces more accurate results relative to state of the art general purpose metagenomic profiling methods, which tend to focus mostly on bacteria. This repository is set up to use BWA to map reads to the full NCBI Virus and Fungi RefSeq databases, and then filter and profile these results using the compute-abundances.py script. Usage is simple and details are provided below. For more details on MiCoP, see the paper link below.

### Paper

The MiCoP paper is currently available on bioRxiv: https://www.biorxiv.org/content/early/2018/01/04/243188

If you use MiCoP, please cite the paper. For instance:
LaPierre, N., Mangul, S., Alser, M., Mandric, I., Wu, N. C., Koslicki, D., & Eskin, E. (2018). MiCoP: Microbial Community Profiling method for detecting viral and fungal organisms in metagenomic samples. *bioRxiv*, 243188.

### Setup

Simply run "./setup.sh".

This will download about 8.4 GB of data, which will be uncompressed into files totalling about 13 GB of disk space. These are the pre-indexed NCBI Virus and Fungi databases for BWA.

### Basic usage

First, use a mapping tool that produces sam output, such as BWA mem. For convenience, we include a copy of BWA kit in our repository (see License Info). The setup script will also download pre-indexed NCBI databases for viruses and fungi. It is strongly suggested that these index files be used for BWA mem, to ensure proper profiling results. For info on how to use BWA, see: http://bio-bwa.sourceforge.net/bwa.shtml#13

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
