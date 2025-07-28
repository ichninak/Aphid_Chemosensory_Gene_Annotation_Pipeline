# Aphid Chemosensory Gene Annotation Pipeline

This project adapts the chemosensory gene annotation pipeline originally developed for ant genomes ([GAGA project](https://github.com/schraderL/GAGA)) to **aphid genomes**. The workflow has been reimplemented using [Nextflow](https://www.nextflow.io/) for better scalability and reproducibility.

It reuses core Perl scripts from the original GAGA repository (e.g., `get_genewise_gtf_corrected.pl`, `run_OR_classification.pl`) with minor modifications, and integrates them into an automated pipeline suitable for aphid species.

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

# Annotation of Odorant Receptor (OR) 
The annotation of the odorant receptor is conducted using 
