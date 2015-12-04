---
title: "RGeneticThesaurus Guide"
author: "Tomasz Konopka"
date: "2015-12-04"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


# RGeneticThesaurus

This vignette describes an application of package RGeneticThesaurus to analyze 
variants and mutations in a toy dataset. 



## Getting started

If you haven't already, install package RGeneticThesaurus from github. (If you need to, uncomment the following lines)


```r
## library("devtools")
## install_github("tkonopka/RGeneticThesaurus")
```

After the installation, load the package 


```r
library("RGeneticThesaurus")
```

## Toy dataset

In what follows, we will use a toy dataset that is bundled with the package. This dataset consists of a synthetic "normal" sample, a matching "tumor" sample, and all associated variant analysis files. We will have available the following:

1. known germline variants in the normal sample;
2. known somatic mutations in the tumor sample;
3. synthetic data and alignments for both the normal and tumor samples;
4. called variants in both the normal and tumor samples;
5. thesaurus annotations for variants in the tumor sample.

The first two categories are known by definition (the toy dataset is a synthetic dataset). Files under points 3 and 4 are bundled with the package for completeness, but we will not use them directly in this vignette. Files under point 5 consist of output from the java program GeneticThesaurus - they include annotated variant call files and associated thesaurus supplementary files.

The objective of the analysis in the vignette will be to recover the known germline and somatic mutations given the called variants and annotations obtained from data.

(In practice, we will look at a couple of regions on chromosome X around gene IKBKG. All variants will be written in hg38 coordinates. )

## Loading expected variants

To get started, let's load the expected germline and somatic variants. These are in files `germline.vcf` and `somatic.vcf` bundled with the RGeneticThesaurus package. We can obtain paths to these files on the local filesystem using `system.file`.


```r
file.germline = system.file("extdata", "germline.vcf", package="RGeneticThesaurus")	
file.somatic = system.file("extdata", "somatic.vcf", package="RGeneticThesaurus")	
```

Next, we can read the contents of these files


```r
germline = readVariantsFromFile(file.germline, getcolumns=5)
somatic = readVariantsFromFile(file.somatic, getcolumns=5)
```

These commands read from files in vcf format and output standard R data.frames. In contrast to `read.table`, the function `readVariantsFromFile` can read parts of a file (above, only the first five columns) and automatically filter out certain rows (by default, indels and entries marked with certain filter codes). A nice feature is that, unlike `read.table`, this function produces (empty) data frames with vcf-like formating even when the input data file contains zero rows.
      
As these new objects are standard data frames, we can view them as usual


```r
germline
```

```
##    chr  position id ref alt
## 1 chrX 154648258  .   C   T
```

```r
somatic
```

```
##    chr  position id ref alt
## 1 chrX 154552029  .   A   G
## 2 chrX 154648184  .   G   A
```

We see that the toy example contains one germline variant and two somatic mutations (btw, in hg38 coordinates).



## Loading thesaurus data

Next, let's load information on variants obtained from data, i.e. variants called from alignments of reads. 

A thesaurus annotation set consists of a file with called variants (extension `.vcf.gz`), a file with thesaurus links (extension `.vtf.gz`), and a file with allele frequency estimates (extension `.baf.tsv.gz`). The package is bundled with two sets of such analyses, one produced using a low threshold on mapping quality and one produced using a high threshold on mapping quality. (See GeneticThesaurus software documentation.)

To load these data from disk into R objects, we first need to obtain their paths using `system.file`.
   

```r
file.mmq1 = system.file("extdata", "tumor.mmq1.thesaurus.vcf.gz", package="RGeneticThesaurus")
file.mmq16 = system.file("extdata", "tumor.mmq16.thesaurus.vcf.gz", package="RGeneticThesaurus")
```

The commands look up paths to called variants for the two analyses. From these, we can load the entire set of thesaurus annotations with `readThesaurusSet`


```r
set.mmq1 = readThesaurusSet(file.mmq1)
set.mmq16 = readThesaurusSet(file.mmq16)
```

The function `readThesaurusSet` here takes as input a file in vcf format, but internally also searches for related files with the other mentioned extensions. The function then loads all the relevant information into a list. We can obtain some basic information about these object just by executing the set name, e.g.


```r
set.mmq1
```

```
## RGeneticThesaurusSet
## variants:	 5 
## links:		 4
```

We see that the set produce with mapping quality threshold 1 contains 5 variants, which altogether are annotated with 4 thesaurus links.

Later, it will be important to know that these set are actually R lists with data frames called `variants`, `links`, and `baf`. We can view each of these data frames 


```r
set.mmq1$variants[,1:7]
```

```
##    chr  position id ref alt   qual    filter
## 1 chrX 154552029  .   A   G 205.20      PASS
## 2 chrX 154556180  .   G   A  38.16 thesaurus
## 3 chrX 154556254  .   C   T  33.51 thesaurus
## 4 chrX 154648184  .   G   A  55.28 thesaurus
## 5 chrX 154648258  .   C   T  54.39 thesaurus
```

```r
set.mmq1$links
```

```
##   chr.from position.from id.from chr.to position.to id.to
## 1     chrX     154556180    <NA>   chrX   154648258  <NA>
## 2     chrX     154556254    <NA>   chrX   154648184  <NA>
## 3     chrX     154648184    <NA>   chrX   154556254  <NA>
## 4     chrX     154648258    <NA>   chrX   154556180  <NA>
```

```r
set.mmq1$baf[,c(1,2,3,4,grep("BAF", colnames(set.mmq1$baf)))]	
```

```
##    chr  position ref alt tumor.naive.BAF tumor.thesaurus.BAF
## 1 chrX 154552029   A   G          1.0000                   1
## 2 chrX 154556180   G   A          0.3810                   1
## 3 chrX 154556254   C   T          0.4118                   1
## 4 chrX 154648184   G   A          0.5652                   1
## 5 chrX 154648258   C   T          0.6316                   1
##   normal.naive.BAF normal.thesaurus.BAF
## 1             0.00                    0
## 2             0.35                    1
## 3             0.00                    0
## 4             0.00                    0
## 5             0.65                    1
```

In the last view, we see B-allele frequency (BAF) estimates for all the variants. There are two samples called `normal` and `tumor` as explained above, and there are two columns per sample. Columns labeled "naive" contain allelic frequencies obtained by counting reads at the designated locus. Columns labeled "thesaurus" contain estimates using reads aligned at relevant loci given in the links table.

*Side note:* the `baf` table also contains other columns with raw read counts. You can explore its contents just like any R data.frame.



## Changes in allele frequency
   
So far, we have just loaded data from disk into the R session. Now we can start analyzing this data.

To call somatic mutations, we can look through the B-allele frequency tables and identify loci where the BAF differs between normal and tumor samples. We can do this using the `callBafChanges` function.

Remember that we have two sets of analyses (`mmq1`, and `mmq16`) and that we have BAF estimates obtained by naive counting and by thesaurus annotations. We can therefore call mutations in at least four ways. We can use the `

A conservative approach is to use `mmq16` data and naive BAF estimates

```r
mut.mmq16.naive = callBafChanges(set.mmq16$baf, "normal.naive.BAF", "tumor.naive.BAF")
```

```
## Error in callBafChanges(set.mmq16$baf, "normal.naive.BAF", "tumor.naive.BAF"): object 'minsample' not found
```

```r
mut.mmq16.naive[,1:4]
```

```
## Error in eval(expr, envir, enclos): object 'mut.mmq16.naive' not found
```
By manual inspection, we see that this approach identifies only one of the two somatic mutations in our dataset. That's alright, but let's see other ways to analyze the data.


The other extreme is to use `mmq1` settings and thesaurus BAF estimates

```r
mut.mmq1.thesaurus = callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF")
```

```
## Error in callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF"): object 'minsample' not found
```

```r
mut.mmq1.thesaurus[,1:4]
```

```
## Error in eval(expr, envir, enclos): object 'mut.mmq1.thesaurus' not found
```
This approach gives three somatic mutation calls, which is one more than we expected. 

In the next section we will see how to evaluate these results automatically.


## Comparing variant sets

Realistic datasets may contain thousands of candidate variants, so it is important to summarize calls in an automated fashion. We will use function `compareVariants`. 

To compare our conservative somatic mutation calls with the expected ones, we perform

```r
compareVariants(somatic, mut.mmq16.naive)
```

```
## Error in compareVariants(somatic, mut.mmq16.naive): object 'mut.mmq16.naive' not found
```
This returns a list of several data frames, which split the genomic locations into true positives (TP), thesaurus true positives (TTP), false positive (FP), and false negatives (FN). The last data frame in the list contains a summary of these categories. For this analysis we see, consistent with what we observed above, that we have one true positive and one false negative. 

To compare our other mutation call set, we evaluate

```r
compareVariants(somatic, mut.mmq1.thesaurus, links=set.mmq1$links)
```

```
## Error in compareVariants(somatic, mut.mmq1.thesaurus, links = set.mmq1$links): object 'mut.mmq1.thesaurus' not found
```
Here, we provide thesaurus links to the comparison procedure. In the final summary we see that there are no false positives and no false negatives, which is a great result. Note that one of the called variants is classified as a thesaurus true positive. This means that this site was technically called incorrectly. However, the thesaurus annotation (in the `links` table) contained a warning and linked this site to another position which had a true mutation.


*Side note:* As an exercise, you can try to compare the known somatic mutations with the mmq1 mutation set, but without using the thesaurus links to help with the variant comparison. 




## Tuning the thesaurus analysis

It is worth experimenting with the settings used in `callBafChanges` as it is not obvious that the defaults are appropriate in all situations. So let's explore this `callBafChanges` function a little further. 


For example, using naive BAF estimates in our mmq1 set, under the default settings we have

```r
nrow(callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF"))
```

```
## Error in callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF"): object 'minsample' not found
```
This means there are three candidate mutation sites. Using thesaurus assisted BAF estimates we have

```r
nrow(callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF"))
```

```
## Error in callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF"): object 'minsample' not found
```
Again, three cadidate hits. 


However, if we use different thresholds in `callBafChanges`, the results may differ. Let's raise the minimum allele frequency in the tumor sample from 0.25 to 0.6. Then we have

```r
nrow(callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF", minAF=0.6))
```

```
## Error in callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF", : object 'minsample' not found
```

```r
nrow(callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF", minAF=0.6))
```

```
## Error in callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF", : object 'minsample' not found
```
Thus, raising the threshold decreased the number of hits under the naive approach but not under the thesaurus approach (The expected allele frequencies of all mutations in the synthetic datasets were 1, i.e. well above the threshold of 0.6). 

Other settings that may affect mutation calls through `callBafChanges` are `foldchange` and `maxref`. Use `help(callBafChanges)` to see information about these settings.



## Session Info


```r
sessionInfo()
```

```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] RGeneticThesaurus_0.1.0.0 knitr_1.10.5             
## 
## loaded via a namespace (and not attached):
## [1] chron_2.3-47     data.table_1.9.6 evaluate_0.7     formatR_1.2     
## [5] magrittr_1.5     markdown_0.7.7   stringi_0.5-5    stringr_1.0.0   
## [9] tools_3.1.2
```
