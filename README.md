# Aphid Chemosensory Gene Annotation Pipeline

This project adapts the chemosensory gene annotation pipeline originally developed for ant genomes ([GAGA project](https://github.com/schraderL/GAGA)) to **aphid genomes**. The workflow has been reimplemented using [Nextflow](https://www.nextflow.io/) for better scalability and reproducibility.

It reuses core Perl scripts from the original GAGA repository (e.g., `get_genewise_gtf_corrected.pl`, `run_OR_classification.pl`) with minor modifications, and integrates them into an automated pipeline suitable for aphid species.

in progress not done 

## Installation & usage 

### 1. Prerequisites
  - A working conda installation (or Docker, etc)
  - Java (17 or newer):
  - [Nextflow](https://www.nextflow.io/docs/latest/install.html)
  - [HAPpy-ABCENTH](https://github.com/biorover/HAPpy-ABCENTH)
  - [Genome tools](https://genometools.org/)
  - Blast
  - InterProScan
  - [GeMoMa v1.7.1](https://www.jstacs.de/index.php/GeMoMa) (Note that newer versions will not work in this pipeline, as the module CompareTranscripts has been modified)

### 2. Environnement Setup
Create a conda environnemet with the required tools:
```bash
conda create -n aphid-nf -c bioconda -c etetoolkit python=3 numpy ete3 wise2 mafft=7 hmmer=3 intervaltree pandas perl
conda activate aphids-nf
```

Installation of Main dependencies:
```bash
conda install -c bioconda -c conda-forge openjdk nextflow genometools-genometools blast interproscan gemoma=1.7.1
```

Installation of HAPpy-ABCENTH:
```bash
pip install HAPpy-ABCENTH
```

### 3. Running the Workflow
Example command:
```bash
nextflow run main.nf --OR true --GR true
(--genome path/to/genome.fa --outdir results/)
```
available modules:
 - `--OR` for ordorant receptor annotation (running the file `OR_annotation`)
 - `--GR` for gustatory receptor annotation (running the files `GR1.nf`, `GR2.nf`, `GR3.nf`)
 - `--OBP` (in development)
> More usage examples and parameters will be documented.

## ORs: Odorant Receptor Annotation
The annotation of the odorant receptor is performed using the [HAPpy-ABCENTH pipeline](https://github.com/biorover/HAPpy-ABCENTH). ABCENTH ("Annotation Based on Conserved Exons Noticed Through Homology") is a gene finder specifically devised for multigene families with extremely high sequence divergence but highly conserved exon structure. It is also designed to avoid gene fusion in tandem arrays.

The `OR_annotation.nf` module runs ABCENTH using HMM profiles built from manually curated OR genes with the command:
```
HAPpy --genome <target genome fasta> --ref_genome <one or more reference genome fastas> \
--annotations <one gtf per ref genome> --cutoff <p distance at which proteins are clustered, 0.45 for ant ORs> \
--search_mode exons --annotator ABCENTH
```
You can use my own HMM files in [dbOR.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/dbOR.zip) created across 44 species of aphids:

<img width="689" height="506" alt="tree_44species" src="https://github.com/user-attachments/assets/8c256946-bbc3-491e-a4db-2f07cdeb2ca7" />

Next, The OR annotation are converted to GFF3 format, and predicted protein sequences are extracted. These proteins are then:
 - Functionnally annotated using InterProScan
 - Compared to reference OR datasets (BLAST) to identify ORco and rename ORs
 - Classified as complete, partial or pseudogene

The sequence datatsets used can be found in [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip).

## GRs: Gustatory Receptor Annotation

The annotation of the gustatory receptor is performed using the [HAPpy-ABCENTH pipeline](https://github.com/biorover/HAPpy-ABCENTH). However, here we performed both ABCENTH and HAPpy method to annotate the GRs, and combine this both methods to generate a final gene model. The are also classified into complete, partial and pseudogenes.

We began by assembling a dataset based on high-quality gene models of gustatory receptors (GRs), which served as the foundation for subsequent GR annotations. This dataset includes GR sequences from 44 aphid species, labeled with the identifiers GAGA-ID: (GAGA-0001 to GAGA-0044). The corresponding protein sequences are stored in the file [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip). the file [GAGA_species_list.txt](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/GAGA_species_list.txt) contains the species names matched to each GAGA code. Additionally, we used the aphids genome annotations to construct exon HMM profiles with the HAPpy-ABCENTH tool. These profiles are compiled in the archive [db2GR.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db2GR.zip).

### 1. GR annotation with ABCENTH methods
`GR1.nf` The first pipeline performs gustatory receptor (GR) annotation using the ABCENTH tool. it requires the genome assembly and the directory containing the previously described HMM profiles. in addition, the path to the [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip) folder must be set (as specified in file [nextflow.config](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/nextflow.config)), along with the OR annotations, which help ensure that no olfactory receptors (ORs) are misidentified as GRs due to their high sequence similarity. All scripts needed to execute this pipeline are available in the [scripts](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/tree/main/script) folder.

### 2. GR annotation with HAPpy method
In this second pipeline, GRs are annotated through homology-based prediction using HAPpy, which relies on GeneWise. The required inputs include the genome assembly and the previously mentioned GR protein dataset [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip). Similar to the first pipeline, it is essential to provide the path to the [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip) directory into [nextflow.config](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/nextflow.config) and the OR annotations. All scripts necessary for running this step are located in the [scripts](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/tree/main/script) folder.

### 3. Merging both method
The third pipeline consolidates the GR models produced in the first two steps to generate a final set of gene models.

## Input Requirements 
- Genome assembly in Fasta format.
- HMM profiles: [dbOR.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/dbOR.zip) and [db2GR.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db2GR.zip).
- [db_chemo.zip](https://github.com/ichninak/Aphid_Chemosensory_Gene_Annotation_Pipeline/blob/main/db_chemo.zip).

## Output Files
Each module produces:
 - Annotated gene in GFF3 format (complete, partial and pseudogene)
 - Peptide FASTA files and nucleotide FASTA files
 - Summary of annotation


## Credits
This pipeline reuses core annotation logic and tools from the [GAGA project](https://github.com/schraderL/GAGA/tree/main), developped for ant genomes.

These components are used within this Nextflow implementation  

## TODO 
- [ ]  Finalize OBP annotation module
- [ ]  Improve modularity of Nextflow processes
- [ ]  Replace or wrap Perl scripts in native DSL processes
- [ ]  Add containerization support (e.g. Docker or Singularity)

## Contact
<ins> Maintainer</ins>: **Ishaac Chninak**

<ins> Institution</ins>: INRAE - UMR IGEEP and University of Rennes

<ins> Email</ins>: ishaac.chninak@gmail.com



#### Feel free to fork, adapt, or contribute improvements to the pipeline!
