## Introduction
This is a collection of custom scripts used for Gong et al. 2020. A Shiny app to browse the placenta transcriptome is also avaible [here](https://www.obgyn.cam.ac.uk/placentome/) and the source code is [here](https://github.com/sung/ShinyPlacentome).

## Initial mapping of RNA-Seq data
The RNA-Seq reads were trimmed (using [cutadapt](https://github.com/marcelm/cutadapt) and [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) and mapped to the primary chromosomal assemblies of the [GRCh38.p3 version](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.29/) of the human reference genome using [TopHat2](https://github.com/infphilo/tophat), a splice-aware mapper built on top of [Bowtie2](https://github.com/BenLangmead/bowtie2) short-read aligner. For more details, read [Gong et al., Epigenetics, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5989156/) and [Gong et al., JCI Insight, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124516/).

## A list of dependent software and database
- [tophat2](https://github.com/infphilo/tophat)
- [bamutils](https://genome.sph.umich.edu/wiki/BamUtil)
- [cutadapt](https://github.com/marcelm/cutadapt)
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)

## Contacts
+ [Sung Gong](https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/)
