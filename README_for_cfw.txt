CFW MOUSE GENOTYPING-BY-SEQUENCING, RNA-SEQ AND COMPLEX TRAIT DATA
May 21, 2016

This Data Dryad package contains data to accompany the following
publication:

  Parker CC, Gopalakrishnan G, Carbonetto P, Gonzales NM, Leung E,
  Park YJ, Aryee E, Davis J, Blizard DA, Ackert-Bicknell CL, Lionikas
  A, Pritchard JK, Palmer AA. Genome-wide association study of
  behavioral, physiological and gene expression traits in commercially
  available outbred CFW mice. To appear in Nature Genetics.

These data were collected as part of a large study to assess the
viability of using Carworth Farms White (CFW) mice for mapping genes
and genetic loci underlying complex traits relevant to the study of
human disease and psychology.

For code implementing QTL mapping of physiological, behavioral and
gene expression phenotypes, and other analyses of these data, refer to
our code repository hosted by Github:

  http://github.com/pcarbo/cfw

To use these data for research, please cite our paper published in
Nature Genetics, as well as the Data Dryad package:

  Parker CC, Gopalakrishnan G, Carbonetto P, Gonzales NM, Leung E,
  Park YJ, Aryee E, Davis J, Blizard DA, Ackert-Bicknell CL, Lionikas
  A, Pritchard JK, Palmer AA. Data from: Genome-wide association study
  of behavioral, physiological and gene expression traits in
  commercially available outbred CFW mice. Dryad Digital
  Repository. http://dx.doi.org/<doi goes here>

CONTENTS
========

1. README.txt: this file.
2. pheno.csv: phenotype data from 1,219 CFW male mice.
3. map.txt: 92,734 SNPs genotyped in CFW mice.
4. geno.txt: Genotype data for 1,161 mice at 92,734 SNPs.
5. rna_HIP.vcf: Genotype calls from RNA-seq in the hippocampus.
6. rna_HIP.vcf: Genotype calls from RNA-seq in the striatum.
7. rna_HIP.vcf: Genotype calls from RNA-seq prefrontal cortex.
8. combined.thr.norm.fpkm: FPKM values from RNA-seq in all
   three tissues.

PHYSIOLOGICAL AND BEHAVIORAL TRAIT DATA
=======================================

File pheno.csv contains physiological and behavioural phenotype data
on 1,219 mice from the Carworth Farms White (CFW) outbred mouse stock.
All mice are males. The data are stored in comma-delimited ("csv")
format, with one line per sample. Missing entries are marked as "NA",
following the convention used in R. Use R function read.pheno in file
read.data.R from the github repository to read these data into an R
data frame.

Individual columns of the table are as follows:

  * id: Unique number assigned to each mouse.
  * round: Refers to a shipment of mice. There are at most 48 mice in
    a single shipment.
  * cageid: Each mouse is placed in a cage with 3 other mice.
  * FCbox, PPIbox, methcage: Boxes and cages used for beahvioral
    testing.
  * discard: A "yes" in this column means that there were heatlh
    concerns or other concerns with this mouse, and we advise
    excluding this sample from further analysis.
  * mixup: A "yes" in this column means the sample is flagged as a
    possible sample mixup due to mishandling or mislabeling of the
    flowcells.
  * glucoseage, methage, FCage, PPIage, sacage: Age of mouse (in days)
    when glucose levels are measured, when various behavioral tests
    are started, and when mouse is sacrificed.
  * bw0, bw1, bw2, bw3, PPIweight, sacweight: Body weights recorded at
    various time points. bw0 is the body weight measured upon arrival
    in Chicago, when tail is snipped. bw1, bw2 and bw3 are body
    weights taken during methamphetamine sensivity testing. PPIweight
    is the body weight measured during prepulse inhibition testing.
    sacweight is the body weight measured when the mouse is sacrificed.
  * BMD: Bone-mineral density.
  * TA, EDL, gastroc, plantaris, soleus: Measured weights of various
    muscles (in mg).
  * tibia: Length of the tibia bone (in mm).
  * abnormalbone: Refers to "bone health", in which 0 = healthy
    looking bone, and 1 = abnormal bone, white and swollen, which,
    might be expression of osteopetrosis. Note that this is not the
    same "abnormal bone" phenotype that we used for QTL mapping.
  * experimenters: Initials of the people who handled the mice.
  * fastglucose: Glucose levels in plasma after fasting. Measured
    using glucometer once at dissection, after fasting overnight.
    Units are milligrams per deciliter (mg/dL), which is the SI unit
    for measuring glucose in blood.
  * PreTrainD1: Average proportion of freezing on day 1 during the
    pre-training interval (30-180 seconds) before exposure to tones
    and shocks ("pre-training freezing").
  * AvToneD1: Average proportion of freezing on the first day during
    exposure to the conditioned stimulus ("freezing to tone on day 1").
  * AvContextD2: Average proportion of time freezing in the 30–180
    second interval on the second day in conditions identical to the
    first day ("freezing to same context").
  * AvAltContextD3: Average proportion of time freezing over the
    30-180 interval on the third day in an altered setting.
  * AvToneD3: Average proportion of time freezing on the third day in
    the altered setting during the 30-second intervals in which the
    tones are presented ("freezing to cue").
  * D3.180, D3.240, D3.300, D3.360: Proportion of freezing on day 3 of
    the conditioned fear tests in the 30-second intervals after
    receiving the tones at 180, 240, 300 and 360 seconds.
  * p120b1, p120b2, p120b3, p120b4, startle: During prepulse inhibition
    (PPI) tests, average startle response during the pulse-alone
    trials (with 120-dB pulses), blocks 1 through 4, and the average
    over blocks 2 and 3 ("startle").
  * pp3b1, pp3b2, pp3avg, pp6b1, pp6b2, pp6avg, pp12b1, pp12b2,
    pp12avg: Average startle response during the prepulse trials with
    3, 6 and 12 dB prepulses followed by a 120-dB pulse.
  * pp3PPIb1, pp3PPIb2, pp3PPIavg, pp6PPIb1, pp6PPIb2, pp6PPIavg
    pp12PPIb1, pp12PPIb2, pp12PPIavg, PPIavg: Average of the
    "inhibition intensity" taken as the ratio of the prepulse response
    during the 3, 6 and 12-dB prepulse trials to the pulse-alone
    startle amplitude, including an average across all amplitudes
    ("PPIavg").
  * avgnostim: Average startle response during trials when no stimulus
    is presented.
  * D1totaldist0to15, D1totaldist15to30, D1totaldist0to30
    D1totaldist0to15, D1totaldist15to30," D1totaldist0to30: "Baseline"
    measurements of locomotor response on Day 1 and Day 2 of
    methamphetamine sensitivity testing. Locomotor activity is
    measured as the distance traveled (in cm) over 15-minute and
    30-minute intervals immediately after receiving an injection of
    saline solution.
  * D1TOTDIST5, D1TOTDIST10, D1TOTDIST15, D1TOTDIST20, D1TOTDIST25,
    D1TOTDIST30, D2TOTDIST5, D2TOTDIST10, D2TOTDIST15, D2TOTDIST20,
    D2TOTDIST25, D2TOTDIST30: Baseline measurements of locomotor
    response over 5-minute intervals on days 1 and 2 of the
    methamphetamine sensitivity testing.
  * D3totaldist0to15, D3totaldist15to30, D3totaldist0to30:
    Measurements of locomotor response to methamphetamine injections,
    which is intended to produce locomotor stimulation. Locomotor
    activity is measured in the same way as on Days 1 and 2,
    immediately after administering drug.
  * D3TOTDIST5, D3TOTDIST10, D3TOTDIST15, D3TOTDIST20, D3TOTDIST25,
    D3TOTDIST30: Measurements of locomotor response to methamphetamine
    injections over 5-minute intervals on day 3 of the testing.
  * D1ctrtime0to15, D1ctrtime0to30, D2ctrtime0to15, D2ctrtime0to30,
    D3ctrtime0to15, D3ctrtime0to30: These are measures of time spent
    in the center of the arena, in seconds. On the first day of
    testing, this provides a measure of baseline response to a novel
    environment. On day 3, this provides a measure of methamphetamine
    sensitivity.
  * D1hact0to15,D1hact0to30, D1vact0to15, D1vact0to30, D2hact0to15,
    D2hact0to30, D2vact0to15, D2vact0to30, D3hact0to15, D3hact0to30,
    D3vact0to15, D3vact0to30: Measures of horizontal and vertical
    activity on days 1-3 of methamphetamine sensitivity testing. The
    day-3 measurements are used to assess methamphetamine sensitivity,
    and the day-1 and day-2 measures provide more general assessments
    of locomotor sensitivity.

GENOTYPE DATA
=============

File map.txt contains information about 92,734 single nucleotide
polymorphisms (SNPs) on chromosomes 1-19 genotyped in the CFW outbred
mouse cohort. The table is stored in space-delimited columns, with one
row per SNP. All genomic positions are based on Mouse Genome Assembly
38 from the NCBI database (mm10, December 2011).

The columns of map.txt are as follows:

  * id: Unique label assigned to each SNP. This is either the dbSNP
    database id (release 137), or if the SNP is not registered in
    dbSNP, then it is an internal id of the form "cfw-N-XXXXXX, where
    N is the chromsome number, and XXXXXX is the base-pair position of
    the SNP on the chromosome.
  * chr: Chromosome number.
  * pos: Base-pair position on chromosome.
  * ref: Reference allele. This always matches the dbSNP reference
    allele.
  * alt: Alternative allele. This usually matches one of the
    alternative alleles in the corresponding dbSNP entry.
  * quality: Imputation quality of each SNP. Here we define
    "imputation quality" as the proportion of samples in which the
    probability of the maximum-probability genotype is >0.5.

The genotype data are stored in file geno.txt. The table with
space-delimited columns contains genotype data for 1,161 mice at
92,734 SNPs. The first column ("id") is the sample id, and the second
column ("discard") indicates whether the sample should be discarded
because of flowcell samples that were mislabeled, and so we cannot be
sure of the identity of these samples.

The remaining columns give the genotypes at all SNPs. The genotypes
are represented as "dosages"; specifically, the expected number of
times the alternative allele is observed in the genotype. This will
either be an integer (0, 1 or 2), or a real number between 0 and 2
when there is some uncertainty in the estimate of the genotype.

GENOTYPES FROM RNA-SEQ DATA
===========================

Using the combined RNA-Seq reads from the three tissues, polymorphic
sites in the transcriptome were found. The file rna_HIP.vcf contains
the genotypes called at these polymorphic sites for the samples that
were used in the hippocmapus RNA-seq experiment. The genotypes are
stored using the standard VCF file format; see, for example,
http://www.1000genomes.org/wiki/Analysis/vcf4.0. Similarly, files
rna_PFC.vcf and rna_STR.vcf contain genotypes called for the samples
that were used in the prefrontal cortex and striatum RNA-seq
experiments.

NORMALIZED FPKM VALUES FOR RNA-SEQ
==================================

File combined.thr.norm.fpkm contains the thresholded, normalized fpkm
values (fragments per kilobase per million mapped reads) generated
from the RNA-seq experiments for the three brain tissues. The data
file represents a matrix, with a header for each row and column. Each
column represents a single gene-tissue combination, and each row
represents a single sample. Columns are labeled using a combination of
the Ensembl gene id and the brain tissue in which the RNA-seq
experiment was conducted; e.g., column ENSMUSG00000051951.HIP contains
expression fpkm values for gene ENSMUSG00000051951 (Xkr4) in the
hippocampus.

Only genes that had a raw (untransformed) fpkm value of greater than 1
are retained for any downstream analysis, so only these genes are
included in this file. The fpkm values for each gene are also
normalized using the quantile transform so that the empirical
distribution of fpkm values aligns with the expected distribution
under the standard normal.

CONTACT INFO
============

For questions and feedback, please contact:

Abraham Palmer
Department of Psychiatry
University of California, San Diego
9500 Gilman Drive
La Jolla, California, USA
aapalmer@ucsd.edu

CONTRIBUTORS
============

Clarissa C. Parker
Shyam Gopalakrishnan
Peter Carbonetto
Natalia M. Gonzales
Emily Leung
Yeonhee J. Park
Emmanuel Aryee
Joe Davis
David A. Blizard
Cheryl L. Ackert-Bicknell
Arimantas Lionikas
Jonathan K. Pritchard
Abraham A. Palmer

ACKNOWLEDGMENTS
===============

This project was funded by NIH R01GM097737 (A.A.P.), NIH T32DA07255
(C.C.P), NIH T32GM07197 (N.M.G.), NIH R01AR056280 (D.A.B.), NIH
R01AR060234 (C.A.B.), the Fellowship from the Human Frontiers Science
Program (P.C.), and the Howard Hughes Medical Institute (J.K.P.). The
authors wish to acknowledge technical assistance from Dana Godfrey,
Sima Lionikaite, Vikte Lionikaite, Ausra S. Lionikiene, and John
Zekos; as well as technical and intellectual input from Drs. Mark
Abney, Justin Borevitz, Karl Broman, Na Cai, Riyan Cheng, Nancy Cox,
Robert Davies, Jonathan Flint, Leo Goodstadt, Paul Grabowski, Bettina
Harr, Ellen Leffler, Richard Mott, Jerome Nicod, John Novembre, Alkes
Price, Matthew Stephens, Daniel Weeks, and Xiang Zhou.
