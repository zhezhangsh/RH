SummarizeSubreadAlignment <- function(gr, seq) {
  require(GenomicRanges);
  
  cov <- cv0 <- coverage(gr)[[1]];
  cov[cov>1] <- 1;
  cvd <- GRanges(seqlevels(gr), IRanges(start(cov), end(cov)));
  cvd <- cvd[runValue(cov)>0];
  smm <- sapply(cvd, function(c) {
    sub <- cv0[start(c):end(c)];
    c(mean(sub), max(sub))
  })
  elementMetadata(cvd) <- data.frame(Length=width(cvd), Mean=round(smm[1, ], 3), Max=smm[2, ]);
    
  gro <- lapply(1:length(cvd), function(i) {
    c <- cvd[i];
    g <- gr[countOverlaps(gr, c)>0];
    g <- g[order(g$as)];
    g <- g[!duplicated(g$ss)];
    g <- g[order(g$ss)];
    g[strand(g) != '*'];
  });
  str <- lapply(gro, function(g) as.vector(unique(strand(g))));
  cvd$Strand400 <- sapply(str, function(s) if (identical(s, '+')) 1 else if (identical(s, '-')) -1 else 0);
  cvd$Count400 <- sapply(gro, length);
  cvd$MinSS <- sapply(gro, function(g) min(g$ss));
  cvd$MaxSE <- sapply(gro, function(g) max(g$se));
}


nonz <- lapply(gr0, function(gr) {
  cov <- coverage(gr)[[1]];
  cov[cov>1] <- 1;
  cvd <- GRanges(seqlevels(gr), IRanges(start(cov), end(cov)));
  cvd[runValue(cov)>0];
})
