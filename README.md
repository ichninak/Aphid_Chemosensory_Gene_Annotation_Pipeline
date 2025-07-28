# Aphids_Annotation_OR_GR
Workflow: a pipeline adapted from Ants annotation. 

in progress not done 

# Installation

## requirement

  - Java (17 or newer):
  - [Nextflow](https://www.nextflow.io/docs/latest/install.html)
  - [HAPpy-ABCENTH](https://github.com/biorover/HAPpy-ABCENTH)
  - [Genome tools](https://genometools.org/)
  - Blast
  - InterProScan
  - [GeMoMa v1.7.1](https://www.jstacs.de/index.php/GeMoMa) (Note that newer versions will not work in this pipeline, as the module CompareTranscripts has been modified)

## conda installation

Installation of essential dependance for the main dependance:

```bash
conda create -n aphids_nf -c bioconda -c etetoolkit python=3 numpy ete3 wise2 mafft=7 hmmer=3 intervaltree pandas perl
conda activate aphids_nf
```

installation of all main dependencie:

```bash
conda install -c bioconda -c conda-forge openjdk nextflow genometools-genometools blast interproscan gemoma=1.7.1
```

installation of HAPpy-ABCENTH:
```bash
pip install HAPpy-ABCENTH
```
