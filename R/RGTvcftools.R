## 
## This file is part of the RGeneticThesaurus package.
##
## Contains helper functions for reading and manipulating vcf files and tables
##
## Author: Tomasz Konpka
##



##' Determines whether a variant is a substitution or an indel
##'
##' Given ref and alt alleles, this function returns a logical vector describing
##' whether or not the variants represent indels. The function regards entries like
##' rb="X" ab="Y" as not an indel. The function regards rb="X" and ab="Y,Z" as not an indel.
##' Other types of ref/alt combinations are reported as indels.
##'
##' @param rb vector of reference bases
##' @param ab vector of alternative bases
##'
##' @export
isIndel = function(rb, ab) {
    
    ## make sure the inputs are characters
    if (class(rb)!="character") {
        rb = as.character(rb);
    }
    if (class(ab)!="character") {
        ab = as.character(ab);
    }
    
    ## check inputs are correctly matched
    if (length(rb)!=length(ab)) {
        stop("reference and alternate base inputs have different lengths\n",call.=FALSE);
    }
    
    ## make objects like REF,ALT
    ans = paste(rb, ab, sep=",");
    ans = strsplit(ans, ",");
    ans = sapply(ans, function(x) {
        ## here first element of x is ref, subsequent are alt alleles
        ncharx = nchar(x);
        ## simple case is when all alleles have one letter
        if (sum(ncharx==1)==length(x)) {
            return(FALSE);
        }
        ## if the lengths are not all the same, that's an indel
        if (length(unique(ncharx))>1) {
            return(TRUE);
        }
        ## remaining elements must be of type e.g. AT -> TA (two singles declared together)
        return(FALSE)
        
    })
    return(as.logical(ans));  
}



##' Sort a data frame of variants by chromosome and position
##'
##' Sort a data frame by chromosome and position. (Note, this function could be faster
##' using data.table. But this implementation avoids dependencies)
##' 
##' @param vcf data frame with variants; should contain columns (chr, position)
##' @param chrorder vector of chromosome names; order of this vector will determine sort order
##' 
##' @export
sortvcf = function(vcf, chrorder=NULL) {
    
    if (is.null(chrorder)) {
        return (vcf);
    }
    
    vcf[,"chr"] = as.character(vcf[,"chr"]);
    vcfs = split(vcf, vcf[,"chr"]);
    
    ans = list();    
    for (i in 1:length(chrorder)) {        
        nowchr = as.character(chrorder[i]);
        nowvcf = vcfs[[nowchr]];        
        if (!is.null(nowvcf)) {    
            nowvcf = nowvcf[order(nowvcf[,"position"]),];      
            ans[[nowchr]] = nowvcf;
        }        
    }
    ans = do.call(rbind, ans);
    rownames(ans) = NULL;
    
    return(ans);
}




##' Read variants from disk 
##'
##' Read variants from disk and perform some pre-processing. In contrast to
##' obtaining a data.frame using read.table(), this function does not crash
##' when a vcf file does not have any valid rows (i.e. no variants detected).
##' This function also automatically allows the user to discard some columns
##' of the input table. 
##' 
##' @param ff filename
##' @param n integer, determines how many lines are read from the input to scan for the
##' vcf header. Set this to a moderatly large number, i.e. larger than the expected number
##' of header lines in the vcf file.
##' @param ignorelines vector of filter codes. Lines in the vcf that contain one
##' of these filter codes will be eliminated, i.e. not reported in the output
##' @param getcolumns number of columns to read/return. Use this to select only a few of the columns of
##' the vcf table (ignoring the info or format fields can save memory)
##' @param withindels logical. Set to TRUE to return all variants in input vcf, or
##' set to FALSE to obtain only single-nucleotide substitutions,
##' 
##' @export
readVariantsFromFile = function(ff, n=4096, ignorelines = c("thesaurushard","thesaurusmany"),
    getcolumns=10, withindels=FALSE) {
    
    ## make an empty table of variants
    empty = data.frame(matrix(0, ncol=10, nrow=0))
    colnames(empty) = c("chr", "position", "id", "ref", "alt", "qual",
                "filter", "info", "format", "sample")
    
    ## First try to read a few lines from a file connection
    
    ## try to read the top of the file to make sure it is not empty
    fcon = file(ff, open="r");
    fcontry = readLines(fcon, n=n);
    close(fcon);
    
    ## look through data and see if there are any data lines (non-header lines)
    fcontry = fcontry[substring(fcontry, 1,1)!="#"];
    if (length(fcontry)==0) {
        return(empty[,1:getcolumns])
    }
    
    ## if reached here, the file is non-empty, so try to read it and output the desired columns
    ans = read.table(ff, stringsAsFactors=F);
    
    ## in some weird cases, read.table will change "T" -> TRUE. Fix this here
    ## This step references columns 4,5 which are ref, alt
    if (class(ans[,4])=="logical") {
        ans[,4] = substring(as.character(ans[,4]),1,1);
    }
    if (class(ans[,5])=="logical") {        
        ans[,5] = substring(as.character(ans[,5]),1,1);
    }
    
    ## perhaps permanently get rid off some lines (according to filter values)
    nowhard = rep(TRUE, nrow(ans));
    for (nowignore in ignorelines) {
        nowhard[grep(nowignore,ans[,7])] = FALSE;
    }
    ans = ans[nowhard, ];
    
    ## perhaps get rid of indels
    if (!withindels) {
        ans = ans[!isIndel(ans[,4], ans[,5]),]
    }

    if (ncol(ans)>getcolumns) {
        ans = ans[,1:getcolumns];
    }
    colnames(ans)[1:getcolumns] = colnames(empty)[1:getcolumns];
    
    
    
    rownames(ans) = NULL;    
    return (ans);
}





##' Read thesaurus links from a file
##'
##' This function reads a vtf file from disk and returns a data.frame. The data.frame is
##' organized into six columns. Each row represents one thesaurus link and has an origin and
##' destination.
##'
##' @param ff filename of vtf
##' @param n number of lines. The function will read the vtf in multiple chunks, each chunk
##' containing at most these many lines.
##'
##' @export
readLinksFromFile = function(ff, n=65536) {
    
    ans = list()
    
    ## try to read the top of the file to make sure it is not empty
    fcon = file(ff, open="r");
    fcontry = readLines(fcon, n=n);
    while (length(fcontry)>0) {
        ## get rid of comment lines
        fcontry = fcontry[substring(fcontry, 1,1)!="#"];
        
        ## convert links into one table
        temp = strsplit(fcontry, "\t")
        temp = lapply(temp, function(x) {
            x1 = unlist(strsplit(x[1],":") )
            x2 = strsplit(x[2:length(x)], ":")
            x2.chr = sapply(x2, function(x) {x[1]})
            x2.pos = sapply(x2, function(x) {as.integer(x[2])})
            x2.id = sapply(x2, function(x) {x[3]} )
            x3 = data.frame(
                chr.from = x1[1], position.from = as.integer(x1[2]), id.from = x1[3],
                chr.to = x2.chr, position.to = x2.pos, id.to = x2.id, 
                stringsAsFactors=F)
            return(x3)
        })
        tempall = data.frame(data.table::rbindlist(temp))
        ans[[length(ans)+1]] = tempall;
        
        ## read more lines
        fcontry = readLines(fcon, n=n);
    }    
    close(fcon);
    
    ans = data.frame(data.table::rbindlist(ans))
    
    ## avoid returning bad format when there is no data.
    if (nrow(ans)==0) {
        ans = data.frame(chr.from="chr1", position.from=4, chr.to="chr1", position.to=5,
            id.from=NA, id.to=NA,
            stringsAsFactors=F)
        ans = ans[c(),]
    }
    ans = unique(ans)
    rownames(ans) = NULL;
    return(ans)
}





##' Read a table with baf from file
##'
##' This function is just a wrapper for read.table(). It is included in the package for
##' completeness so that all thesaurus-related files can be accessed with similar function names.
##'
##' @param ff filename
##' 
##' @export
readBafFromFile = function(ff) {
    ans = read.table(ff, header=T, stringsAsFactors=F);
    if (colnames(ans)[1] != "chr" | colnames(ans)[2] != "position") {
        stop("Contents of file does not start with chr position, is it a BAF file?\n");
    }
    ## replace all NAs or NaNs in the BAF table, if any, with 0s
    whichBAF = grep("BAF$", colnames(ans));
    for (nowcol in whichBAF) {
        ans[,nowcol] = as.numeric(ans[,nowcol]);
        ans[!is.finite(ans[,nowcol]), nowcol] = 0;
    }
    ## that's all!
    ans = unique(ans)
    rownames(ans) = NULL;
    return(ans);
}





##' Read a set of variants annotated using GeneticThesaurus
##'
##' The GeneticThesaurus software for annotating variants with thesaurus links
##' typically creates three output files - a .vcf file with called variants, a .vtf file
##' with links, and a baf.tsv file with allelic frequencies. This function here
##' reads these three file types together into a single R object (a list with the contents
##' of the three files in data frames).
##'
##' Note: the individual files are often large, so this function may take up to a few minutes
##' to complete.
##' 
##' @param variantsfile filename for a vcf file annotated with the genetic thesaurus.
##' This function assumes variantsfile ends with vcf (or vcf.gz), and then looks for files
##' with similar name with extensions vtf (or vtf.gz) and baf.tsv (or baf.tsv.gz)
##' @param n number of lines to read at a time (used by readVariantsFromFile and readLinksFromFile
##' @param ignorelines filter codes to remove from variant list
##' @param getcolumns integer; this is passed on to function readVariantsFromFile
##' @param withindels logical; determines whether variant list should include insertions/deletions
##' 
##' @export
readThesaurusSet = function(variantsfile, n=4096, ignorelines = c("thesaurushard","thesaurusmany"),
    getcolumns=10, withindels=FALSE) {    

    ## from the variantsfile input, determine filenames for files making up this thesaurus set
    if (grepl(".vcf$", variantsfile)) {
        fileroot = substring(variantsfile, 1, nchar(variantsfile)-4)
        reqfiles = paste0(fileroot, c(".vcf", ".vtf", ".baf.tsv"));
    } else if (grepl(".vcf.gz$", variantsfile)) {
        fileroot = substring(variantsfile, 1, nchar(variantsfile)-7)
        reqfiles = paste0(fileroot, c(".vcf.gz", ".vtf.gz", ".baf.tsv.gz"));
    } else {
        msg = paste0("Cannot understand format of variantsfile: ",variantsfile,"\n");
        msg = paste0(msg, "Please provide a name ending in vcf or vcf.gz\n");
        stop(msg, call.=FALSE);
    }
        
    ## check if the required files exist
    allok = sum(file.exists(reqfiles))==3;
    if (!allok) {
        for (i in 1:3) {
            if (!file.exists(reqfiles[i])) {
                cat(paste0("Cannot find file ",reqfiles[i],".\n",
                           "Note: to load a part of a thesaurus annotation set, \n",
                           "use readVariantsFromFile(), readLinksFromFile(), or readBafFromFile()\n"));
            }
        }    
        stop("Missing files\n", call. = FALSE);
    }
    
    ## if reached here, the three files exist, so load them    
    ans = list();
    ans$files = reqfiles;
    names(ans$files) = c("variants", "links", "baf");
    ## get variants from file on disk
    ans$variants = readVariantsFromFile(reqfiles[1], n=n, ignorelines=ignorelines,
        getcolumns=getcolumns, withindels=withindels);
    ## get thesaurus links and baf from files on disk
    ans$links = readLinksFromFile(reqfiles[2], n=32*n);
    ans$baf = readBafFromFile(reqfiles[3]);
    
    ## if "ignorelines" is not empty, some lines in vcf are potentially deleted from ans$variants.
    ## Here, subset the links and baf tables so that
    ## they contain only those lines that are relevant to the current ans$variants table
    if (length(ignorelines)>0 & nrow(ans$variants)>0) {
        varlabels = paste0(ans$variants[,"chr"], ":", ans$variants[,"position"]);
        ## subset the links and baf tables
        ans$baf = ans$baf[paste0(ans$baf[,"chr"],":",ans$baf[,"position"]) %in% varlabels,];
        ans$links = ans$links[paste0(ans$links[,"chr.from"],":", ans$links[,"position.from"])
            %in% varlabels, ];
        rownames(ans$baf) = NULL;
        rownames(ans$links) = NULL;
    }
    
    class(ans) = "RGeneticThesaurusSet";
    
    return(ans);    
}






##' Display information about an RGeneticThesaurus object
##' 
##' @param x an object of class RGeneticThesaurus to print
##' @param ... not used
##'
##' @export
print.RGeneticThesaurusSet = function(x, ...) {
    
    if (class(x)!="RGeneticThesaurusSet") {
        stop("Input object is not of class RGeneticThesaurusSet\n")
    }

    ## print out basic information about the thesaurus set onto screen
    cat("RGeneticThesaurusSet\n");
    cat("variants:\t",nrow(x$variants),"\n");
    cat("links:\t\t", nrow(x$links),"\n");    
}



##
## This function is DEPRECATED / not used
##' Read a table with thesaurus clusters from file
##'
##' This function is just a wrapper for read.table(). It is included in the package for
##' completeness so that all thesaurus-related files can be accessed with similar function names.
##'
##' @param ff filename
##' 
readClustersFromFile = function(ff) {
    ans = read.table(ff, header=T, stringsAsFactors=F);
    if (colnames(ans)[1] != "chr" | colnames(ans)[2] != "position" | colnames(ans)[3] != "cluster") {
        stop("Contents of file does not start with chr position cluster. Is if a cluster file?\n")
    }
    return(ans);
}






##' Compare a set of variants (hits) with a set of expected variants (known sites).
##'
##' This function compares two sets of variants and reports the number of true positives,
##' false positives and false negatives. 
##'
##' @param known data frame of known/expected variants; should contain column names (chr, position)
##' @param hits data frame of hits/observed variants, should contain column names (chr, position)
##' @param links four-column table with thesaurus links as obtained by function getLinksFromFile().
##' These links are assumed to correspond to variants given via the data frame hits.
##' 
##' @export
compareVariants = function(known, hits, links = NULL) {

    ## make a dataframe based on "hits", but with just the columns of interest
    smallhits = unique(hits[,c("chr", "position")])
    colnames(smallhits) = c("chr.from", "position.from")
    ## make a small data frame based on known (note unique to remove double liners)
    known = unique(known[,c("chr", "position")])

    ## useful vector for column labels 
    chrpos = c("chr", "position");
    ## useful dataframe with chr, position, but no rows
    emptyhits = smallhits[c(),]
    colnames(emptyhits) = chrpos;
    
    ## get strings describing hits and expected mutations
    hitcodes = paste0(hits[,"chr"], ":", hits[,"position"])
    knowncodes = paste0(known[,"chr"], ":", known[,"position"])
    
    known.found = rep(FALSE, length(knowncodes))
    names(known.found) = knowncodes

    ## append columns to vtf table indicating whether the from/to positions are present
    ## in the data frame of known sites
    if (is.null(links)) {
        links = data.frame(chr.from="chr1", position.from=0, chr.to="chr1", position.to=0, #
            site.known=FALSE, link.known=FALSE)
        links = links[c(), ]
    } else {
        links[,"site.known"] = (paste0(links[,"chr.from"],":",links[,"position.from"]) %in% knowncodes);
        links[,"link.known"] = (paste0(links[,"chr.to"],":",links[,"position.to"]) %in% knowncodes);
    }
    
    ## merge the smallhits witht vtf
    hitsvtf = merge(smallhits, links, all.x = TRUE)
    
    if (nrow(hitsvtf)>0) {
        hitsvtf[is.na(hitsvtf[,"site.known"]), "site.known"] = FALSE;
        hitsvtf[is.na(hitsvtf[,"link.known"]), "link.known"] = FALSE;    
        hitsvtf[,"site.known"] = hitsvtf[,"site.known"] |
            (paste0(hitsvtf[,"chr.from"], ":", hitsvtf[,"position.from"]) %in% knowncodes);
    }
    
    ## identify sites in the hits that match the known sites exactly
    sites.found = hitsvtf[,"site.known"]
    sites.found.codes = paste0(hitsvtf[,"chr.from"],":", hitsvtf[,"position.from"])[sites.found]
    if (length(sites.found.codes)>0) {
        known.found[sites.found.codes] = TRUE;
        TPhits = unique(hitsvtf[sites.found, paste0(chrpos,".from")]);
        colnames(TPhits) = chrpos;
    } else {
        TPhits = emptyhits;
    }
    
    ## identify sites in the hits that are linked to known sites
    links.found = hitsvtf[,"link.known"] & !hitsvtf[,"site.known"];
    links.found.codes  = paste0(hitsvtf[,"chr.to"],":", hitsvtf[,"position.to"])[links.found];
    if (length(links.found.codes)>0) {
        known.found[links.found.codes] = TRUE;
        TTPhits = unique(hitsvtf[links.found, paste0(chrpos, ".from")])
        colnames(TTPhits) = chrpos;        
    } else {
        TTPhits = emptyhits;
    }
    
    ## find those rows that are not and do not link to known
    goodhits = hitsvtf[hitsvtf[,"site.known"] | hitsvtf[,"link.known"],]
    goodhits = unique(paste0(goodhits[,"chr.from"],":", goodhits[,"position.from"]))
    badhits = hitsvtf[!hitsvtf[,"site.known"] & !hitsvtf[,"link.known"],]
    badhits = unique(paste0(badhits[,"chr.from"],":",badhits[,"position.from"]))    
    badhits = unique(badhits[!(badhits %in% goodhits)])
    badhits = badhits[badhits!=":"]; ## badhits becomes ";" when the data.frame above is empty
    
    ## helper function to convert between a vector of codes (e.g. chr4:339920)
    ## into a dataframe with chr and position columns (chr4 330020)
    code2df = function(codevec) {
        if (length(codevec)>0) {            
            temp = strsplit(codevec, ":")
            temp1 = sapply(temp, function(x) {x[1]});
            temp2 = sapply(temp, function(x) {x[2]});
            ans = data.frame(chr=temp1, position=as.integer(temp2), stringsAsFactors=F);
            return(unique(ans));
        }
        ## if vector does not have at least one element, return an empty data frame
        temp = data.frame(chr="aa", position=0);
        return(temp[c(),])
    }
    
    ## convert vectors of codes into data frames
    FNhits = code2df(names(known.found[!known.found]));
    FPhits = code2df(badhits);
    
    return(
        list(TP=TPhits, TTP=TTPhits, FP=FPhits, FN=FNhits,
             summary = data.frame(
                 TP=nrow(TPhits), TTP=nrow(TTPhits),
                 FP=nrow(FPhits), FN=nrow(FNhits)))
        )
    
}





##' Call sites showing BAF changes
##'
##' This function performs very crude mutation calling based on BAFs. The calling is
##' depends only on fold change and minimum AF thresholds.
##' 
##' @param baftable data.frame with variants and allelic frequencies
##' @param cref identifier for column with reference BAF (e.g. normal tissue)
##' @param csample identifier for column with matched sample BAF (e.g. tumor tissue)
##' @param foldchange minimum fold.change threshold to call mutation
##' @param maxref numeric, maximum allowed BAF in the reference sample
##' @param minAF numeric, minimum BAF required in the sample
##'
##' @export
callBafChanges = function(baftable, cref, csample,
    foldchange = 1.2, maxref = 0.05, minAF=0.15) {
    
    ## sanity checks
    ## calculations assume for fold change is greater than 1
    if (foldchange<1) {
        stop("foldchange threshold must be greater than one\n");        
    }
    ## check that the columns are presented and formated as expected
    if (!(cref %in% colnames(baftable))) {
        stop(paste0("Column ",cref, " is not present in the input table\n"))
    }
    if (!(csample %in% colnames(baftable))) {
        stop(paste0("Column ",csample, " is not present in the input table\n"))
    }    
    if (class(baftable[,cref]) != "numeric" ) {
        stop("cref column is not numeric\n");
    }
    if (class(baftable[,csample]) != "numeric") {
        stop("csample column is not numeric\n");
    }
    
    ## require a certain fold change in allelic frequency, also certain absolute thresholds
    hit1 = (baftable[,cref] > foldchange*baftable[,csample]) &
        (baftable[,csample]<maxref) & (baftable[,cref]>minAF)
    hit2 = (baftable[,csample] > foldchange*baftable[,cref]) &
        (baftable[,cref]<maxref) & (baftable[,csample]>minAF)
    
    ## report the entries in the input table that meet the criteria
    return(baftable[(hit1 | hit2), ,drop=FALSE])        
}





##' Call sites de novo mutations from family trio data
##'
##' This function performs crude mutation calling based on BAF. It searches for variants
##' in a child sample that are not present in either parent. The calculation depends
##' on changes in AF between child and parents.
##'
##' @param baftable data.frame with variants and allelic frequency
##' @param cref1, string, code refererring to one of the matched controls, i.e. one of the parents.
##'                      The function assumes the baftable contains columns such as
##'                      cref1.thesaurus.BAF, etc as output by GeneticThesaurus software
##' @param cref2 string, similar to cref1
##' @param csample string, code referring to the child sample. The function assumes the baftable
##'                  contains columns such as csample.thesaurus.BAF, etc. as output
##'                  by GeneticThesaurus software.
##' @param thesaurus boolean, set TRUE to use thesaurus adjusted BAF and coverage estimates
##' @param foldchange numeric, required fold change in allelic frequency
##' @param maxref numeric, maximal allelic frequency allowed in the reference samples
##' @param minAF numeric, minimum allelic frequency required in the child sample
##' @param mincov numeric, minimum coverage required in the reference samples
##' @param exceptcovchr1 vector of chromosome names that should be excluded from the
##'                       minimum coverage requirement (e.g. set chrY for a female sample)
##' @param exceptcovchr2 vector, similar to exceptcovchr1
##'
##' @export
callTrioDeNovo = function(baftable, cref1="mother", cref2="father", csample="child",
    thesaurus=FALSE, foldchange=1.2, maxref=0.05, minAF=0.15, mincov=10,
    exceptcovchr1=c(), exceptcovchr2=c()) {
    
    ## obtain a true/false vector selecting items from baftable as candidates
       
    ## use the input data to determine which columns of the baftable should be used in the analysis
    if (thesaurus) {
        ## Define columns with thesaurus-adjusted BAF estimates
        cref1baf = paste0(cref1,".thesaurus.BAF");
        cref2baf = paste0(cref2,".thesaurus.BAF");
        csamplebaf = paste0(csample, ".thesaurus.BAF");
        
        ## Define columns with thesaurus-adjusted coverage estimates
        ## (this is complicated because of a bug in GeneticThesaurus (now fixed)
        ## with spelling error in columns)
        cref1cov = paste0(cref1,".thesarus.cov");
        cref2cov = paste0(cref2,".thesarus.cov");
        if (!(cref1cov %in% colnames(baftable))) {
            cref1cov = paste0(cref1,".thesaurus.cov");
        }
        if (!(cref2cov %in% colnames(baftable))) {
            cref2cov = paste0(cref2,".thesaurus.cov");
        }
        
        ## for thesaurus analysis, require naive and adjusted estimates to be different
        cthescov = paste0(csample,".thesarus.cov");
        if (!(cthescov %in% colnames(baftable))) {
            cthescov = paste0(csample,".thesaurus.cov");
        }        
        hitcovdiff =
            (baftable[,paste0(csample,".naive.cov")] != baftable[,cthescov]) |
                (baftable[,"thesaurus.synonyms"]==0);
        
    } else {
        ## Define columns for naive BAF estimates
        cref1baf = paste0(cref1,".naive.BAF");
        cref2baf = paste0(cref2,".naive.BAF");
        csamplebaf = paste0(csample, ".naive.BAF");

        ## Define columns for coverage
        cref1cov = paste0(cref1,".naive.cov");
        cref2cov = paste0(cref2,".naive.cov");

        ## thesaurus-adjusted calculation has an extra condition, here avoid
        ## that condition by allowing all candidates to pass that requirement
        hitcovdiff = rep(TRUE, nrow(baftable));
    }

    ## select hits based on AF changes (comared to reference 1)
    hit1 = (baftable[,cref1baf] > foldchange*baftable[,csamplebaf]) &
        (baftable[,csamplebaf]<maxref) & (baftable[,cref1baf]>minAF)
    hit2 = (baftable[,csamplebaf] > foldchange*baftable[,cref1baf]) &
        (baftable[,cref1baf]<maxref) & (baftable[,csamplebaf]>minAF)
    ## select hits based on AF changes (comared to reference 2)
    hit3 = (baftable[,cref2baf] > foldchange*baftable[,csamplebaf]) &
        (baftable[,csamplebaf]<maxref) & (baftable[,cref2baf]>minAF)
    hit4 = (baftable[,csamplebaf] > foldchange*baftable[,cref2baf]) &
        (baftable[,cref2baf]<maxref) & (baftable[,csamplebaf]>minAF)
    
    ## when user asks to ignore some chromosomes in coverage thresholds,
    ## implement this by temporarily adjusting coverage numbers
    baftable[baftable[,1] %in% exceptcovchr1, cref1cov] = mincov+1;
    baftable[baftable[,1] %in% exceptcovchr2, cref2cov] = mincov+1;
    ## select hits based on minimum coverage in reference samples
    hitcov = (baftable[,cref1cov]>mincov) & (baftable[,cref2cov]>mincov)

    ## merge all the criteria together
    selectvec = (hit1 | hit2) & (hit3 | hit4) & hitcov & hitcovdiff;
    
    ## report the entries in the input table that meet the criteria
    return(baftable[selectvec, ,drop=FALSE])        
}


