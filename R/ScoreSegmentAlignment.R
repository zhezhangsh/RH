# Calculate alignment scores of segments of a sequence to a reference sequence
ScoreSegmentAlignment <- function(seq1, seq0, window.size, step.size) {
  require('Biostrings');
  require('GenomicRanges');
  
  len <- nchar(seq1);
  loc <- seq(1, len, step.size);
  ind <- cbind(loc-window.size, loc+window.size);
  ind[, 1] <- pmax(1, ind[, 1]);
  ind[, 2] <- pmin(len, ind[, 2]);
  
  seq <- DNAStringSet(apply(ind, 1, function(i) subseq(seq1, i[1], i[2])));
  pw  <- pairwiseAlignment(seq, seq0, type='overlap');
  sc  <- score(pw);
  lc  <- Views(pw);
  
  out <- cbind(location=loc, start=ind[, 1], end=ind[, 2], 
               start_ref=start(lc), end_ref=end(lc), score=sc);
  out;
};
