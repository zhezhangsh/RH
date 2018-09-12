LigateSubread <- function(seqs, hole=NA, ind1=NA, ind2=NA) {
  
  if (identical(hole, NA) | identical(ind1, NA) | identical(ind2, NA)) {
    ids <- strsplit(names(seqs), '/');
    ind <- sapply(ids, function(i) rev(i)[1]);
    ind <- strsplit(ind, '_');
    
    hole <- sapply(ids, function(i) rev(i)[2]);
    ind1 <- as.integer(sapply(ind, function(i) i[1]));
    ind2 <- as.integer(sapply(ind, function(i) i[2]));
  };
  
  ids <- sort(unique(hole));
  
  seq <- sapply(ids, function(i) { 
    ind <- which(hole==i);
    stt <- ind1[ind];
    end <- ind2[ind];
    sub <- seqs[ind]; 
    
    ord <- order(stt); 
    stt <- stt[ord];
    end <- end[ord];
    sub <- sub[ord];
    
    seq <- lapply(1:length(stt), function(j) {
      n <- stt[j] - c(0, end)[j];
      c(strrep('N', n), sub[j]); 
    });
    seq <- unlist(seq, use.names = FALSE);
    paste(seq, collapse='');
  });
  names(seq) <- ids;
  
  seq;
}