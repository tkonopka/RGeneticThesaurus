## ------------------------------------------------------------------------
library("devtools")
## install_github("tkonopka/RGeneticThesaurus")

## ------------------------------------------------------------------------
library("RGeneticThesaurus")

## ------------------------------------------------------------------------
file.germline = system.file("extdata", "germline.vcf", package="RGeneticThesaurus")	
file.somatic = system.file("extdata", "somatic.vcf", package="RGeneticThesaurus")	

## ------------------------------------------------------------------------
germline = readVariantsFromFile(file.germline, getcolumns=5)
somatic = readVariantsFromFile(file.somatic, getcolumns=5)

## ------------------------------------------------------------------------
germline
somatic

## ------------------------------------------------------------------------
file.mmq1 = system.file("extdata", "tumor.mmq1.thesaurus.vcf.gz", package="RGeneticThesaurus")
file.mmq16 = system.file("extdata", "tumor.mmq16.thesaurus.vcf.gz", package="RGeneticThesaurus")

## ------------------------------------------------------------------------
set.mmq1 = readThesaurusSet(file.mmq1)
set.mmq16 = readThesaurusSet(file.mmq16)

## ------------------------------------------------------------------------
set.mmq1

## ------------------------------------------------------------------------
set.mmq1$variants[,1:7]
set.mmq1$links
set.mmq1$baf[,c(1,2,3,4,grep("BAF", colnames(set.mmq1$baf)))]	

## ------------------------------------------------------------------------
mut.mmq16.naive = callBafChanges(set.mmq16$baf, "normal.naive.BAF", "tumor.naive.BAF")
mut.mmq16.naive[,1:4]

## ------------------------------------------------------------------------
mut.mmq1.thesaurus = callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF")
mut.mmq1.thesaurus[,1:4]

## ------------------------------------------------------------------------
compareVariants(somatic, mut.mmq16.naive)

## ------------------------------------------------------------------------
compareVariants(somatic, mut.mmq1.thesaurus, links=set.mmq1$links)

## ------------------------------------------------------------------------
nrow(callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF"))

## ------------------------------------------------------------------------
nrow(callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF"))

## ------------------------------------------------------------------------
nrow(callBafChanges(set.mmq1$baf, "normal.naive.BAF", "tumor.naive.BAF", minsample=0.6))
nrow(callBafChanges(set.mmq1$baf, "normal.thesaurus.BAF", "tumor.thesaurus.BAF", minsample=0.6))

## ------------------------------------------------------------------------
sessionInfo()

