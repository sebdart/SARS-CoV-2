# Department of Defense SARS-CoV-2 Bioinformatic Protocol and Clade System

## Table of Content

This document is a living document and is updated as often as needed with more detailed information depending on specific requests from the SME's community. Please contact us if you have any question and/or how to make this protocol a better tool for all.

Last update: *8<sup>th</sup> December 2020*.

----

* [Getting the Data](#getting-the-data)
  * [From GISAID](#from-gisaid)
  * [From Genbank](#from-genbank)
* [Curating and Reformatting Fasta and Metadata prior to bioinformatic modeling](#curating-and-reformatting-fasta-and-metadata-prior-to-bioinformatic-modeling)
* [Sequence Alignement](#sequence-alignement)
* [Nucleotide Site Masking and Minimum Length](#site-masking-and-minimum-length)
* [Tree Topology](#tree-topology)
* [Tree Dating](#tree-dating)
* [DOD SARS-CoV-2 Clade System](#dod-sars-cov-2-clade-system)
  * [Why](#why)
  * [How does it compare with other SARS-CoV-2 clades](#how-does-it-compare-with-other-sars-cov-2-clades)
  * [How](#how)

----


## Getting the Data
### From GISAID

Please, login to your [GISAID account](https://www.epicov.org/epi3/). If you need a new account, [register here](https://www.gisaid.org/registration/register/). 

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

Alternatively, one can also download an easy FASTA ready for NextStain Augur bioinformatic pipeline in clicking on *EPiCoV*, *Download*, scrolling down to *Genomic epidemiology*, and clicking on *nexstfasta*.

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

|         Step           |                                                                                  |
|------------------------|----------------------------------------------------------------------------------|
|    *Step 1 of 3:*      |     Click on *CSV format* under ***Current table view result*** (Right-side)     |
|    *Step 2 of 3:*      |     Click on *Download All Records*                                              |
|    *Step 3 of 3:*      |     Click on *Select All*                                                        |

> You will download a file named *sequences.csv*; this file will have to be processed and curated by an R-code script to be bioinformatically compliant with NextStrain's Augur requirement (or any other phylodyanmic/phylogeographic codes, such as BEAST 1 and 2).

**Get the FASTA Data**

|         Step           |                                                                                   |
|------------------------|------ ----------------------------------------------------------------------------|
|    *Step 1 of 3:*      |     Click on *Nucleotide* under ***Sequence data (FASTA Format)*** (Left-side)    |
|    *Step 2 of 3:*      |     Click on *Download All Records*                                               |
|    *Step 3 of 3:*      |     Click on *Build Custom*                                                       |
|                        |     *Remove* all default                                                          |
|                        |     *Add* in this order: *Geo Location* and then *Accession*                      |


> You will download a file named *sequences.fasta*; this file will have to be processed and curated by a Linux Bash script to be bioinformatically compliant with NextStrain's Augur (or any other phylodyanmic/phylogeographic codes, such as BEAST 1 or 2).

## Curating and Reformatting Fasta and Metadata prior to bioinformatic modeling

Depending on how and from where you download the raw data, some curating will be required. We wrote bash and R-scripts to set the SARS-CoV-2 Fasta data and metaData to be easily ran by phylodynamic bioinformatic codes.

A critical aspect is to ensure that the metadata file of a given sample/sequence is named in a similar manner in the FASTA file:
  - From GenBank: each sample is labelled *Country Name/Accession Number* in both the metaCata and Fasta files.
  - From GISAID: each sample is labelled *Country Name/ID/Year* in both the metaData and Fasta files.

In the Fasta file, all second duplicates by IDs (ie., two nucleotide sequences with same ID) and all second duplicates by sequences (ie., two different IDs with the exact same nucleotide sequence) are eliminated (ie., the first is kept). All nucleotide sequences are reformatted on one single line (ie., half of the total number of lines in the Fasta file is the total number of SARS-CoV-2 samples available in the Fasta file).

In the metData file from GenBank only, all incomplete dates and/or unknown Country names are eliminated. Great care is taken to ensure that all GenBank locations for each sample is written as **Continent(Region):Country:Division(State,Province):Location(County,City)**.

Further in the bioinformatic pipeline, any sample with incomplete dates from GISAID (eg., *YYYY*, *YYYY-MM*, *YYYY-XX*) will be also eliminated. An complete date is defined as *YYYY-MM-dd*.

Bash and R-codes ara available on [DoD/NGA GITHUB](https://gitlab.gs.mil/Dartevelle.Sebastien.1503290509/ncov) with a CAC card or upon request.

## Sequence Alignement

Multiple genomic sequence alignment is performed by the Multiple Alignment using Fast Fourier Transform [MAFFT algorithm v7.471](https://mafft.cbrc.jp/alignment/software/) in Linux.

Each sequences are initially grouped by pack of 250 sequences for small simulations up to 2500 sequences for larger simulations (> 10,000 data points). All MAFFT alignement is performed in SMP parallel:

```shell
$  mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 8   input > output
```

with,
* [reorder](https://mafft.cbrc.jp/alignment/software/manual/manual.html#lbAK).
* [anysymbol](https://mafft.cbrc.jp/alignment/software/anysymbol.html).
* [nomemsave](https://mafft.cbrc.jp/alignment/software/tips.html).
* [adjustdirection](https://mafft.cbrc.jp/alignment/software/adjustdirection.html).
* [thread n](https://mafft.cbrc.jp/alignment/software/multithreading.html).

An older version of [MAFFT Manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html) is available as well as [Tips](https://mafft.cbrc.jp/alignment/software/tips0.html) and [Change Log](https://mafft.cbrc.jp/alignment/software/changelog.html) with its latest updates.

All 250/2500 individual pack of sequences will be eventually aggregated together into one large fasta file of aligned sequences: *results/aligned.fasta*.

```shell
cat results/split_alignments/1.fasta results/split_alignments/0.fasta results/split_alignments/3.fasta results/split_alignments/6.fasta results/split_alignments/4.fasta results/split_alignments/8.fasta > results/aligned.fasta
```

.

## Site Masking and Minimum Length

This is an important aspect in the DoD's protocol because it leads to slightly different results in the branching and clade naming in SARS-CoV-2 phylogenetic trees when compared with other methodologies which do not aggressively mask *chalenging* sites.

We follow the recommendations from [European Molecular Biology Laboratory, UK](https://www.ebi.ac.uk/) as detailed in [their paper](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) and their [GITHUB](https://github.com/W-L/ProblematicSites_SARS-CoV2). Specifically, many mutations seen several times along the phylogenetic tree (ie., highly homoplasic) may be likely the result of contamination, recurrent sequencing errors, or hypermutability, than selection or recombination. Some homoplasic substitutions seem laboratory-specific, suggesting that they might arise from specific combinations of sample preparation, sequencing technology, and consensus calling approaches. Hence, these sites subejct too "suspect mutations" are masked (ie., ignored) in our bioinfomrtaic modeling as decribed in the Table below.

As commonly done, we systematically mask positions 1–55 and 29804–29903.

We aslo filter out any sequences that have too few resolved characters (ie., less than ***29,400*** reference bases) instead of the more common threshold of *29,000* (unless requested/specified otherwise).


| Header          | Description                    |
|-----------------|--------------------------------|
|POS              | 1-based position of the variation on the reference |
|REF              | Reference base |
|ALT              | List of alternative alleles at the position (IUPAC ambiguity code) |
|GENE             | Position falls into range of this gene |
|AA_POS           | Position of amino acid residue within gene |
|AA_REF           | Reference amino acid residue |
|AA_ALT           | List of alternative amino acid residues (IUPAC ambiguity code) |


|    POS    |REF|    ALT    |FILTER |                          JUSTIFICATION                                |   GENE    |AA_POS|AA_REF|  AA_ALT   |
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


## Tree Topology

The Phylogenetic tree inference is based on *Maximum Likelihood inference* methodology with [IQ-TREE v2.0.5](http://www.iqtree.org/) 

```shell
$ iqtree -s results/all/subsampled_alignment.fasta -m GTR -ninit 5 -n 4 -czb  -nt n
```

* -s: aligned sequences (in this case a likely subsample from all many data available from either GISAID or GenBank).

* -m: consensus [substitution model](https://en.wikipedia.org/wiki/Substitution_model) used for SARS-CoV-2, *General Time Reversible model* from [Tavare, 1986](http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf).

* -n: Specify number of iterations to stop.

* -ninit: Specify the number of initial parsimony trees; the Phylogenetic Likelihood Library (PLL) from [Flouri et al., 2015](https://academic.oup.com/sysbio/article/64/2/356/1630375) is used, which generates a random tree (randomized stepwise addition order parsimony tree) with subsequent optimized likelihood evaluation functions and topological rearrangement operations on the Tree (eg., Subtree Pruning and Regrafting-SPR).

* -czb: collapse near zero branches, so that the final tree may be multifurcating. This is useful for bootstrapping in the presence of polytomy to reduce bootstrap supports of short branches. It is a huge time saver for the next step (Tree Dating); however, the "czb" algorithm is not parallelized and rather slow.

* -nt: Specify the number of CPU cores to be used by PLL.

A complete list of all [iqtree commands is available](http://www.iqtree.org/doc/Command-Reference#general-options) as well as a detailed online [iqtree Manual](http://www.iqtree.org/doc/).

Note that iqtree, by default, saves the Tree in NEWICK format as *.treefile*, however this step is overridden by augur pipeline, which saves the output as *tree_raw.nwk*, which will be reused for the next step (ie., Tree Dating).


## Tree Dating

The Phylogenetic tree dating is based on *Maximum Likelihood inference* from [TreeTime algorithm v0.8.0](https://treetime.readthedocs.io/en/latest/). This step involves rooting to a reference sequence, working on polytomies (if any), dating, and refining the final Tree topolgy.

```shell
$ treetime --tree results/all/tree_raw.nwk --aln results/all/subsampled_alignment.fasta --dates data/metadata.tsv --reroot seq1 seq2 --keep-polytomies False --coalescent const --confidence True --covariation False --clock-filter 4 --branch-length-mode auto
```

Rooting of SARS-CoV-2 with respect to the following two sequences:

* In GISAID: `--reroot Wuhan/Hu-1/2019 Wuhan/WH01/2019`, with GISAID Accession Numbers: *EPI_ISL_402125* & *EPI_ISL_406798*
* In GenBank: `--reroot China/MN908947 China/LR757998`, with GenBank Accession Numbers: *MN908947* & *LR757998*

Note that TreeTime, by default, saves the dated Tree in NEWICK format in a specified *--outdir*, however this step is overridden by augur pipeline, which saves the output as *tree.nwk* under the *results/all/* directory.

We do not normally impose a constant rate of mutation and its standard deviation, which will be estimated from our largest samples available by [least-squares regression by TreeTime](https://treetime.readthedocs.io/en/latest/tutorials/clock.html). However, for very small simulation (ie., sample less than < 10,000 data points), we set a constant rates: `--clock-rate 0.0008 --clock-std-dev 0.0004` as commonly done in the community.

A complete list of all available TreeTime commands are [available online](https://treetime.readthedocs.io/en/latest/commands.html).

TreeTime is a powerful tool for refining a given phylogenetic Tree, its topology and dating it; however, it is not parallelized and is not adequate for very large simulations (ie., > 20,000 data points). Hence, we are planing to eventually replace TreeTime by another Tree dating algorithm: [LSD: Least-Squares methods to estimate rates and Dates from serial phylogenies](https://github.com/tothuhien/lsd2), available by default in [iqtree](http://www.iqtree.org/doc/Dating) itself and also in [R language](https://github.com/tothuhien/Rlsd2).

## DOD SARS-CoV-2 Clade System
### Why

The Department of Defense Clade system of SARS-CoV-2 phylogeny aimed to be easy to understand, universal, and reliable regardless of the size of your dataset (ie., tens of thousands data sample vs. a few thousands samples).

Very roughly, one can create a clade system based either on specific mutations at a specific position (nucleotide based) or a specific mutation at a specific position at a given time (time based).

### How does it compare with other SARS-CoV-2 clades

* Pangolin Clades (nucleotide based) can be--somehow--[visualized here](https://nextstrain.org/ncov/global?c=pangolin_lineage&d=tree,entropy,frequencies&p=full) on NextStrain website (about 3,500 samples).

* NexStrain Clades (time based) can be [visualized here](https://nextstrain.org/ncov/global?c=clade_membership&d=tree,entropy,frequencies&p=full) (about 3,500 samples).

* GISAID Clades (nucleotide based) can be [visualized here](https://nextstrain.org/ncov/global?c=GISAID_clade&d=tree,entropy,frequencies&p=full) on Nextstrain website (about 3,500 samples).

* DOD Clades (nucleotide based) can be [visualized here](https://sars-cov-2.dev.east.paas.nga.mil/sars-cov2?d=tree,entropy,frequencies&p=full) on NGA website with GISAID only data (about 15,000 samples).

* DOD Clades (nucleotide based) can be [visualized here](https://sars-cov-2.dev.east.paas.nga.mil/sars-cov2?d=tree,entropy,frequencies&p=full) on NGA website with GenBank only data (about 15,000 samples).

We eliminated the time-based clade system because the first 2020 "20A" branch is within 15 days from December 2019; depending on the time-dating model assumptions (ie., constant mutation clock throughout the year), depending on your site masking assumptions (see above), and depending on how many data points used in your model (3,564 vs. 16,427 vs 35,000 samples), the "20A" clade shifts towards the end of December 2019, which would make under that time-based logic a new clade "19C" ... Of course, all models have their own assumptions but this would create confusion---So, we decided to move to a nucleotide based type of classification. After many modeling runs over 35,000 samples, we furthermore decided to use a modified GISAID classification system for two reasons:

* First is simplicity: the DOD Clades are based on mutation from "C to T" on two different sites (241 and 3037) early on the phylogeny Tree of SARS-CoV-2; the C-Clade is based on nucleotide C241/C3037 and T-Clade is nucleotide 241T/3037T, the rest (C.1, C.2, T.1, T.2) easily flows with a similar logic (see below); it is universal no matter haw many data points you are adding; indeed,

* Second, well precisely, it is universal: the other clades systems start to get mixed up across the phylogentic tree when your dataset increase in size; at first (ie., small dataset), it looks good with "only" 3,564 data points but then it gets messier and messier as you add many more data points; the worse is when you hit many tens of thousands of data point ... One can have an idea of this from our work when we calculated the two clade systems with these larger samples (DOD vs. GISAID clades):

  * [GISAID dataset with GISAID clades (14,859 samples)](https://sars-cov-2.dev.east.paas.nga.mil/sars-cov-2?c=GISAID_clade&d=tree,entropy,frequencies&p=full); one can see with this larger dataset that the Clade "O" and "G" start to become distributed across the whole Tree, regardless of its Clade membership; it gets worse when you hit 25,000 samples and on.
  * [The same GISAID dataset with the DOD clades (14,859 samples)](https://sars-cov-2.dev.east.paas.nga.mil/sars-cov-2?d=tree,entropy,frequencies&p=full); One can visualize with this same (larger) dataset that the Clade T (and T.1, T.2) and C (and C.1, C.2) do not get inter-mixed as the dataset increases in size.

DoD remains open to all sorts of better suggestions and setups; nothing is set in stones; we remain flexible with the best path forward for a better SARS-CoV-2 Clade system.

### Where

The DOD SARS-CoV-2 clades are defined as follow (to be defined in NextStrain's *clades.tsv* file):

|   clade   | nuc sites  |   alt     |
|-----------|------------|-----------|
|   C       |   241      |	  C      |
|   C       |   3037	   |    C      |
|           |            |           |
|   T       |   241      |    T      |
|   T       |   3037	   |    T      |
|           |            |           |
|   C.1     |   241	     |    C      |
|   C.1     |   3037	   |    C      |
|   C.1     |   8782	   |    T      |
|           |            |           |
|   C.2     |   241	     |    C      |
|   C.2     |   3037	   |    C      |
|   C.2     |   1397	   |    A      |
|           |            |           |
|   T.1     |   241	     |    T      |
|   T.1     |   3037	   |    T      |
|   T.1     |   23403	   |    G      |
|   T.1     |   25563	   |    T      |
|           |            |           |
|   T.2     |   241	     |    T      |
|   T.2     |   3037	   |    T      |
|   T.2     |   23403	   |    G      |
|   T.2     |   28882	   |    A      |
