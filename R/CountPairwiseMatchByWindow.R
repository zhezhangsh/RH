CountPairwiseMatchByWindow <- function(seq1, seq2, size=1000, start=NA, end=NA) {
  
  mp   <- c('A'=1, 'C'=2, 'G'=3, 'T'=4, '-'=0);
  s1   <- strsplit(as.character(seq1), '')[[1]];
  s2   <- strsplit(as.character(seq2), '')[[1]];
  mtch <- cbind(Seq1=mp[s1], Seq2=mp[s2]);

  if (identical(NA, start)) start <- 1;
  if (identical(NA, end)) end <- nrow(mtch);
  
  ind1 <- which(mtch[, 1]!=0);
  ind2 <- which(mtch[, 2]!=0);
  
  countMatchPercent <- function(aligned, ind, ind.pos, window, column) {
    i1  <- ind[max(1, ind.pos-window)];
    i2  <- ind[min(length(ind), ind.pos+window)];
    sub <- aligned[i1:i2, , drop=FALSE];
    n1  <- nrow(sub);
    n2  <- nrow(sub[sub[, column]!=0, , drop=FALSE]);
    n3  <- nrow(sub[sub[, 1]==sub[, 2], , drop=FALSE]);
    c(n1, n2, n3); 
  }; 
  
  cnt1 <- sapply(1:length(ind1), function(i) countMatchPercent(mtch, ind1, i, ceiling(size/2), 1));
  cnt2 <- sapply(1:length(ind2), function(i) countMatchPercent(mtch, ind2, i, ceiling(size/2), 1));
  
  cnt1 <- cbind(Position=ind1, Len_Alignment=cnt1[1, ], Len_Base=cnt1[2, ], Len_Match=cnt1[3, ]);
  cnt2 <- cbind(Position=ind2, Len_Alignment=cnt2[1, ], Len_Base=cnt2[2, ], Len_Match=cnt2[3, ]);
  
  list(seq1=cnt1, seq2=cnt2);
}