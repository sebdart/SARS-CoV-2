# Department of Defense SARS-CoV-2 bioinformatic Protocol and Clade System

## Table of Content

This document is a life document and is updated as often as needed with more detailed information depending on specific requests from the SME's community. Please contact us if you have any question and/or how to make this protocol a better tool for all.

Last update: *20<sup>th</sup> November 2020*.

----

* [Getting the Data](#getting-the-data)
  * [From GISAID](#from-gisaid)
  * [From Genbank](#from-genbank)
* [Curing and Reformatting the Data and Metadata prior to bioinformatic modeling](#curing-and-reformatting-the-data-and-metadata-prior-to-bioinformatic-modeling)
* [Nucleotide Site Masking](#nucleotide-site-masking)
* [Tree Construction](#tree-construction)
* [Tree Dating](#tree-dating)
* [DOD SARS-CoV-2 Clade System](#dod-sars-cov-2-clade-system)
  * [Why](#why)
  * [How](#how)
  * [How does it compare with other SARS-CoV-2 clades](#how-does-it-compare-with-other-sars-cov-2-clades)

----


## Getting the Data
### From GISAID

You must have an account with [GISAID](https://www.gisaid.org/registration/register/).
If you hav one, log in to [GISAID](https://www.epicov.org/epi3/)

From there:

**Get the MetaData as prepared by NextStrain**

Click on *EPiCoV*, then *Download*, scroll down to *Genomic epidemiology*
Click on *nexstmeta*

> You will download a file named *metadata_YYYY-MM-dd_##-##.tsv.gz* with a TSV formatted metada file named *metadata.tsv*; this file is ready to be used by NextStrain's Augur bioinformatic pipeline.

**Get the FASTA Data**

Click on *EPiCoV*, then *Search*
From there, click on:
- *complete*: only considers genomic sample with at least 29,000 base pairs

- *high-coverage*: genomic samples with less than 1% Ns (undefined bases), less thann 0.05% unique Amino Acid mutations (not seen in any other sequences database) and with no insertion or deletion unless verified by submitter

- *low-coverage*: excludes all genomic samples with > 5% N

Select or search the specific sequences you want to download into a FASTA file; if you want them all select the click box on the upper left; then click the *Download* button on the lower-right and select *Sequences (FASTA)*:

> You will download a file named *gisaid_hcov-19_YYYY_MM_dd_##.fasta*

Alternatively, one can also download an easy FASTA ready for NextStain Augur bioinformatic pipeline in clicking on *EPiCoV*, *Download*, scrolling down to *Genomic epidemiology*, and clickin on *nexstfasta*.

> You will download a file named *sequences_YYYY-MM-dd_##-##.fasta.gz* with a FASTA formatted file named *sequences.fasta*


### From Genbank

[GenBank](https://www.ncbi.nlm.nih.gov/genbank/) is a publicly available genetic sequence database, a password free, fully annotated collection of all publicly available DNA sequences data repository maintained by the [US Department of Health & Human Services (HHS)](https://www.hhs.gov/)/[National Institutes of Health (NIH)](https://www.nih.gov/)/[National Center for Biotechnology Information  (NCBI)](https://www.ncbi.nlm.nih.gov/).

SARS-CoV-2 data can be easily accessed and downloaded from [GenBank NCBI virus Portal](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) The Taxonomic ID of SARS-CoV-2 is [2697049](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049).

SARS-CoV-2 Genbank's genomic data and metadata will have to go through two scripts to cure these files for easy use by phylodynamic/phylogeographic bioinformatic models (such as NextStrain's Augur), it is therefore important to follow all the instructions herewith:

On the Left-Side menu:
- *Virus*, ensure you have *2697049*, which is the GenBank TAXid of SARS-CoV-2;
- *Nucleotide Completeness*, click on *Complete*

Click on *Download* (Upper-Right Blue button)

**Get the MetaData**

*Step 1 of 3:*
Click on *CSV format* under ***Current table view result*** (Right-side)

*Step 2 of 3:*
Click on *Download All Records*

*Step 3 of 3:*
Click on *Select All*

> You will download a file named *sequences.csv*; this file will have to be processed and curred by an R-code script to be bioinformatically compliant with NextStrain's Augur requirement (or any other phylodyanmic/phylogeographic codes, such as BEAST 1 and 2).

**Get the FASTA Data**

*Step 1 of 3:*
Click on *Nucleotide* under ***Sequence data (FASTA Format)*** (Left-side)

*Step 2 of 3:*
Click on *Download All Records*

*Step 3 of 3:*
Click on *Build Custom*
*Remove* all default
*Add* in this order *Geo Location* and *Accession*

> You will download a file named *sequences.fasta*; this file will have to be processed and curred by a Linux Bash script to be bioinformatically compliant with NextStrain's Augur (or any other phylodyanmic/phylogeographic codes, such as BEAST 1 or 2).


## Curing and Reformatting the Data and Metadata prior to bioinformatic modeling

     About DOD's minimum nucleotide completeness in SARS-CoV-2 sequences

## Nucleotide Site Masking

We follow the recommendation in masking high sites

|    POS    |REF|    ALT    |FILTER |                                 JUST                                  |   GENE    |AA_POS|AA_REF|  AA_ALT   |
|-----------|---|-----------|-------|-----------------------------------------------------------------------|-----------|------|------|-----------|
|1-55       |.  |.          |mask   |seq_end                                                                |.          |.     |.     |.          |
|150        |T  |C,Y        |mask   |homoplasic<br>single_src<br>neighbour_linked                           |.          |.     |.     |.          |
|153        |T  |G,Y        |mask   |homoplasic<br>single_src<br>neighbour_linked                           |.          |.     |.     |.          |
|635        |C  |Y,T        |mask   |highly_ambiguous<br>homoplasic<br>narrow_src                           |gene-orf1ab|124   |R     |X,C        |
|1707       |C  |Y,T,A      |mask   |highly_ambiguous<br>homoplasic<br>narrow_src                           |gene-orf1ab|481   |S     |X,F,Y      |
|1895       |G  |T,K        |mask   |ambiguous                                                              |gene-orf1ab|544   |V     |L,X        |
|2091       |C  |T,Y        |mask   |highly_ambiguous<br>homoplasic<br>narrow_src                           |gene-orf1ab|609   |T     |I,X        |
|2094       |C  |T,Y        |mask   |highly_ambiguous<br>narrow_src                                         |gene-orf1ab|610   |S     |L,X        |
|2198       |G  |R,T,A      |mask   |homoplasic<br>ambiguous<br>narrow_src                                  |gene-orf1ab|645   |G     |X,C,S      |
|2604       |G  |T,K        |mask   |homoplasic<br>ambiguous<br>single_src                                  |gene-orf1ab|780   |G     |V,X        |
|3145       |G  |T          |mask   |homoplasic<br>single_src                                               |gene-orf1ab|960   |L     |F          |
|3564       |G  |T,K        |mask   |highly_ambiguous<br>highly_homoplasic<br>single_src                    |gene-orf1ab|1100  |G     |V,X        |
|3639       |G  |R,A        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|1125  |G     |X,D        |
|3778       |A  |G          |mask   |homoplasic<br>single_src                                               |gene-orf1ab|1171  |T     |T          |
|4050       |A  |C          |mask   |homoplasic<br>single_src<br>amended                                    |gene-orf1ab|1262  |N     |T          |
|5011       |A  |M,C        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|1582  |Q     |X,H        |
|5257       |A  |W,G,M      |mask   |highly_ambiguous<br>single_src                                         |gene-orf1ab|1664  |L     |X,L,X      |
|5736       |C  |T,Y        |mask   |homoplasic<br>single_src                                               |gene-orf1ab|1824  |A     |V,X        |
|5743       |G  |S,T        |mask   |highly_ambiguous<br>neighbour_linked<br>narrow_src                     |gene-orf1ab|1826  |E     |X,D        |
|5744       |T  |Y          |mask   |highly_ambiguous<br>neighbour_linked<br>narrow_src                     |gene-orf1ab|1827  |Y     |X          |
|6167       |G  |R,K        |mask   |ambiguous<br>single_src                                                |gene-orf1ab|1968  |V     |X,X        |
|6255       |C  |T          |mask   |highly_homoplasic<br>narrow_src                                        |gene-orf1ab|1997  |A     |V          |
|6869       |A  |W,T        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|2202  |T     |X,S        |
|8022       |T  |K,G,C      |mask   |highly_ambiguous<br>highly_homoplasic<br>narrow_src                    |gene-orf1ab|2586  |V     |X,G,A      |
|8026       |A  |W,G,T      |mask   |ambiguous<br>homoplasic<br>narrow_src                                  |gene-orf1ab|2587  |A     |A,A,A      |
|8790       |G  |T,K        |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-orf1ab|2842  |G     |V,X        |
|8827       |T  |W          |mask   |ambiguous<br>neighbour_linked                                          |gene-orf1ab|2854  |I     |I          |
|8828       |G  |K,A        |mask   |ambiguous<br>neighbour_linked                                          |gene-orf1ab|2855  |A     |X,T        |
|9039       |C  |Y,M,T      |mask   |ambiguous<br>single_src                                                |gene-orf1ab|2925  |A     |X,X,V      |
|10129      |T  |Y,C,W      |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-orf1ab|3288  |T     |T,T,T      |
|10239      |C  |M,A        |mask   |homoplasic<br>single_src                                               |gene-orf1ab|3325  |S     |X,Y        |
|11074      |C  |T,Y        |mask   |highly_homoplasic                                                      |gene-orf1ab|3603  |F     |F,F        |
|11083      |G  |T,K,A      |mask   |highly_homoplasic                                                      |gene-orf1ab|3606  |L     |F,X,L      |
|11535      |G  |T,K        |mask   |ambiguous<br>highly_homoplasic<br>single_src                           |gene-orf1ab|3757  |G     |V,X        |
|13402      |T  |G,K,C      |mask   |homoplasic<br>narrow_src<br>amended                                    |gene-orf1ab|4379  |Y     |*,X,Y      |
|13408      |T  |K          |mask   |homoplasic<br>narrow_src<br>amended                                    |gene-orf1ab|4381  |C     |X          |
|13476      |C  |T,Y        |mask   |highly_ambiguous<br>narrow_src                                         |gene-orf1ab|4404  |C     |V,X        |
|13571      |G  |T,K        |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-orf1ab|4436  |G     |F,X        |
|14277      |G  |T,K        |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-orf1ab|4671  |R     |V,X        |
|15435      |A  |R,G        |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-orf1ab|5057  |E     |X,R        |
|15922      |T  |Y,C        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|5219  |V     |C,C        |
|16290      |T  |K          |mask   |highly_ambiguous<br>single_src                                         |gene-orf1ab|5342  |A     |X          |
|16887      |C  |T,Y        |mask   |highly_homoplasic                                                      |gene-orf1ab|5541  |Y     |I,X        |
|19298      |A  |D,W,T,R    |mask   |highly_ambiguous<br>single_src                                         |gene-orf1ab|6345  |Y     |X,X,L,X    |
|19299      |T  |Y,G        |mask   |homoplasic<br>ambiguous<br>narrow_src                                  |gene-orf1ab|6345  |Y     |X,R        |
|19484      |C  |T,Y        |mask   |highly_ambiguous<br>amended                                            |gene-orf1ab|6407  |A     |L,L        |
|19548      |A  |R,W,G,T,D  |mask   |ambiguous<br>single_src                                                |gene-orf1ab|6428  |S     |X,X,R,L,X  |
|20056      |G  |K,A        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|6597  |E     |X,K        |
|20123      |T  |Y,C        |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-orf1ab|6620  |I     |L,L        |
|20465      |A  |G,R,W      |mask   |highly_homoplasic<br>single_src                                        |gene-orf1ab|6734  |D     |V,X,X      |
|21550      |A  |M,C        |mask   |ambiguous<br>homoplasic<br>narrow_src                                  |gene-orf1ab|7095  |N     |T,T        |
|21551      |A  |W,T,R      |mask   |ambiguous<br>homoplasic<br>narrow_src                                  |gene-orf1ab|7096  |N     |X,S,X      |
|21575      |C  |T,Y        |mask   |highly_homoplasic                                                      |gene-S     |5     |L     |F,X        |
|22335      |G  |T,K,A      |mask   |highly_ambiguous<br>amended                                            |gene-S     |258   |W     |L,X,*      |
|22516      |T  |W,D,K,C,Y  |mask   |highly_ambiguous<br>single_src                                         |gene-S     |318   |F     |X,X,X,F,F  |
|22521      |T  |K,Y        |mask   |highly_ambiguous                                                       |gene-S     |320   |V     |X,X        |
|22661      |G  |T,S,A      |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-S     |367   |V     |F,X,I      |
|22802      |C  |M,A,T,Y    |mask   |homoplasic<br>single_src<br>amended<br>interspecific_contamination     |gene-S     |414   |Q     |X,K,*,X    |
|24389      |A  |W,M,C      |mask   |highly_homoplasic<br>highly_ambiguous<br>narrow_src<br>nanopore_adapter|gene-S     |943   |S     |X,X,R      |
|24390      |G  |T,S,K,C    |mask   |highly_homoplasic<br>highly_ambiguous<br>narrow_src<br>nanopore_adapter|gene-S     |943   |S     |I,X,X,T    |
|24622      |T  |Y          |mask   |highly_ambiguous<br>single_src                                         |gene-S     |1020  |A     |A          |
|24933      |G  |K,T        |mask   |highly_ambiguous<br>highly_homoplasic<br>narrow_src                    |gene-S     |1124  |G     |X,V        |
|25202      |T  |K          |mask   |highly_ambiguous                                                       |gene-S     |1214  |W     |X          |
|25381      |A  |C,R,M      |mask   |homoplasic<br>single_src                                               |gene-S     |1273  |T     |T,T,T      |
|26549      |C  |T,Y        |mask   |homoplasic<br>single_src                                               |gene-M     |9     |T     |T,T        |
|27760      |T  |K,Y,A      |mask   |homoplasic<br>neighbour_linked                                         |.          |.     |.     |.          |
|27761      |T  |C          |mask   |homoplasic<br>neighbour_linked                                         |.          |.     |.     |.          |
|27784      |A  |W,G,T      |mask   |ambiguous<br>single_src                                                |.          |.     |.     |.          |
|28253      |C  |Y,T,G,A,S  |mask   |highly_homoplasic                                                      |gene-ORF8  |120   |F     |F,F,L,L,X  |
|28985      |G  |R,D,T,K,A  |mask   |highly_ambiguous<br>homoplasic<br>single_src                           |gene-N     |238   |G     |X,X,C,X,S  |
|29037      |C  |T,Y,M      |mask   |homoplasic<br>ambiguous<br>single_src                                  |gene-N     |255   |S     |F,X,X      |
|29039      |A  |T,W        |mask   |homoplasic<br>ambiguous<br>single_src                                  |gene-N     |256   |K     |*,X        |
|29425      |G  |T,S,A      |mask   |ambiguous<br>homoplasic<br>single_src                                  |gene-N     |384   |Q     |H,X,Q      |
|29553      |G  |A,T        |mask   |highly_homoplasic<br>single_src                                        |.          |.     |.     |.          |
|29827      |A  |G,T        |mask   |seq_end<br>highly_homoplasic<br>single_src                             |.          |.     |.     |.          |
|29830      |G  |T,A,C      |mask   |seq_end<br>highly_homoplasic<br>single_src                             |.          |.     |.     |.          |
|29804-29903|.  |.          |mask   |seq_end                                                                |.          |.     |.     |.          |


| Header          | Description                    |
|-----------------|--------------------------------|
|POS              | 1-based position of the variation on the reference |
|REF              | Reference base |
|ALT              | List of alternative alleles at the position (IUPAC ambiguity code) |
|JUST             | List of reasons for suggested exclusion (tags described in separate table) |
|GENE             | Position falls into range of this gene |
|AA_POS           | Position of amino acid residue within gene |
|REF_AA           | Reference amino acid residue |
|ALT_AA           | List of alternative amino acid residues (IUPAC ambiguity code) |


Descriptions of reasons for mask/caution are as follows:

| Tag                         | Description                                                      |
|-----------------------------|------------------------------------------------------------------|
| seq_end                     | Alignment ends are affected by low coverage and high error rates |
| ambiguous                   | Sites which show an excess of ambiguous basecalls relative to the number of alternative alleles, often emerging from a single country or sequencing laboratory |
| amended                     | Previous sequencing errors which now appear to have been fixed in the latest versions of the GISAID sequences, at least in sequences from some of the sequencing laboratories |
| highly_ambiguous            | Sites with a very high proportion of ambiguous characters, relative to the number of alternative alleles |
| highly_homoplasic           | Positions which are extremely homoplasic - it is sometimes not necessarily clear if these are hypermutable sites or sequencing artefacts |
| homoplasic                  | Homoplasic sites, with many mutation events needed to explain a relatively small alternative allele count |
| interspecific_contamination | Cases (so far only one instance) in which the known sequencing issue is due to contamination from genetic material that does not have SARS-CoV-2 origin |
| nanopore_adapter            | Cases in which the known sequencing issue is due to the adapter sequences in nanopore reads |
| narrow_src                  | Variants which are found in sequences from only a few sequencing labs (usually two or three), possibly as a consequence of the same artefact reproduced independently |
| neighbour_linked            | Proximal variants displaying near perfect linkage |
| single_src                  | Only observed in samples from a single laboratory |



## Tree Construction

### Tree Dating

## DOD SARS-CoV-2 Clade System
### Why
### How
### How does it compare with other SARS-CoV-2 clades
