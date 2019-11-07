ParseBlasrSam <- function(fname, primary.only=TRUE) {
  require(GenomicAlignments);
  
  aln <- read.table(fname, sep='\t', comment.char = '@', header = FALSE, stringsAsFactors = FALSE);
  if (primary.only) aln <- aln[aln[, 2]==0 | aln[, 2]==16, , drop=FALSE];
  
  as <- as.integer(sub('AS:i:', '', aln[, 17]));
  qw <- as.integer(sub('qe:i:', '', aln[, 14]));
  
  rs <- aln[, 4];
  re <- cigarWidthAlongReferenceSpace(aln[, 6], flag=aln[, 2]) + rs - 1;
  
  aq <- cigarRangesAlongQuerySpace(aln[, 6], flag=aln[, 2], ops=setdiff(CIGAR_OPS, c('H', 'S')), with.ops = TRUE);
  qs <- sapply(aq, function(q) min(start(q)));
  qe <- sapply(aq, function(q) max(end(q)));
  
  str <- rep('+', nrow(aln));
  str[aln[,2]==16 | aln[,2]==272] <- '-';
  
  ids <- strsplit(aln[, 1], '/'); 
  lib <- sapply(ids, function(i) i[1]);
  zmw <- sapply(ids, function(i) i[2]);
  loc <- strsplit(sapply(ids, function(i) i[3]), '_');
  stt <- as.integer(sapply(loc, function(l) l[1]));
  end <- as.integer(sapply(loc, function(l) l[2]));
  
  aln <- data.frame(subread=aln[, 1], read=zmw, sample=lib, reference=aln[, 3], strand=str, score=as, sstart=stt,  send=end,
                    qstart=qs, qend=qe, rstart=rs, rend=re, sequence=aln[, 10], stringsAsFactors = FALSE);
  rownames(aln) <- 1:nrow(aln); 
  
  if (nrow(aln[aln$strand=='-', , drop=FALSE])>0) {
    aln$sequence[aln$strand=='-'] <- as.character(reverseComplement(DNAStringSet(aln$sequence[aln$strand=='-'])));
    len <- aln[aln$strand=='-', 'send'] - aln[aln$strand=='-', 'sstart'] + 1;
    aln[aln$strand=='-', c('qstart', 'qend')] <- cbind(len-aln[aln$strand=='-', 'qend'], len-aln[aln$strand=='-', 'qstart']);
  }
  aln$sequence <- substr(aln$sequence, aln$qstart, aln$qend);
  
  aln; 
}



