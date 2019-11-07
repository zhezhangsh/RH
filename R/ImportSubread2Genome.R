## Import results from PacBio subreads to whole genome alignment using blasr
###############################################################################
# If format==1, tab-delimited text file with the following columns
# 1. library ID
# 2. read ID
# 3. subread start (actual base - 1)
# 4. subread end
# 5. query start
# 6. query end
# 7. reference start
# 8. reference end
# 9. reference chromosome
#10. reference strand: "1/-1"
#11. blasr score (negative integer)
#12. primary or secondary alignment: "p/s"
###############################################################################
ImportSubread2Genome <- function(fnm, format=1) {
  format <- as.integer(format[1]);
  
  if (format==1) {
    aln <- read.delim(fnm, sep='\t', header=FALSE, stringsAsFactors=FALSE, 
                      colClasses=c(rep('character', 2), rep('integer', 6), 'character', rep('integer', 2), 'character'));
    aln[, 12] <- as.integer(aln[, 12]=='p');
    aln <- aln[, c(1:10, 12, 11)];
  }
  
  colnames(aln) <- c('library', 'read', 'sub_start', 'sub_end', 'query_start', 'query_end', 'ref_start', 'ref_end', 
                     'chromosome', 'strand', 'primary', 'score');
  
  aln;
}