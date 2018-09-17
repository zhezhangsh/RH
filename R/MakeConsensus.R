WeightConsensus <- function(mtrx, method=c('rank', 'r', 'count')) {
  ltt <- unique(as.vector(mtrx));
  tbl <- sum <- matrix(0, nr=nrow(mtrx), nc=length(ltt), dimnames = list(rownames(mtrx), ltt));
  ttl <- rep(0, nrow(mtrx));
  rng <- apply(mtrx, 2, function(a) range(which(a!='-')));
  rng <- c(max(rng[1, ]), min(rng[2, ]));
  sub <- mtrx[rng[1]:rng[2], , drop=FALSE];
  for (i in 1:ncol(mtrx)) for (j in 1:length(ltt)) sum[mtrx[, i]==ltt[j], j] <- sum[mtrx[, i]==ltt[j], j] + 1;
  cnt <- apply(sub, 2, function(s) {
    n <- rep(0, length(s));
    for (i in 1:length(ltt)) n[s==ltt[i]] <- sum[rng[1]:rng[2], , drop=FALSE][s==ltt[i], ltt[i]];
    n - 1;
  });
  if (tolower(method)[1] == 'rank') weight <- rank(-colMeans(cnt)) else
    if (tolower(method)[1] == 'r') weight <- colMeans(cnt)/(ncol(cnt)-1) else
      if (tolower(method)[1] == 'count') weight <- colMeans(cnt) else
        weight <- rep(1, ncol(cnt));

  for (i in 1:ncol(mtrx)) {
    x <- mtrx[, i];
    rng <- range(which(x!='-'));
    ttl[rng[1]:rng[2]] <- ttl[rng[1]:rng[2]] + 1;
    for (j in 1:length(ltt)) {
      ind <- which(x==ltt[j]);
      ind <- ind[ind>=rng[1] & ind<=rng[2]];
      tbl[ind, j] <- tbl[ind, j] + weight[i];
    };
  };
  rownames(tbl) <- 1:nrow(tbl);

  out <- cbind(tbl, Count=ttl, Gap=apply(mtrx, 1, function(x) length(x[x=='-'])));

  list(weighted=out, weight=weight, method=method);
};

CallConsensus <- function(weighted) {
  tbl <- weighted[, c('-', 'A', 'C', 'G', 'T')];
  cll <- apply(tbl[, 1:5], 1, function(x) {
    mx <- which(x==max(x));
    if (length(mx)>1) sample(mx, 1) else mx;
  });
  cll <- colnames(tbl)[1:5][cll];
  names(cll) <- rownames(weighted);
  cll;
}

AlignConsensus <- function(con, ref) {
  require(S4Vectors);
  require(Biostrings);

  pw <- pairwiseAlignment(con, ref, type='global-local');

  ptn <- as.character(pattern(pw));
  sbj <- as.character(subject(pw));
  prs <- data.frame(P=strsplit(ptn, '')[[1]], S=strsplit(sbj, '')[[1]], stringsAsFactors = FALSE);
  ind <- rep(NA, nrow(prs));
  ind[prs[,1]!='-'] <- 1:length(prs[prs[,1]!='-', 1]);
  rownames(prs) <- 1:nrow(prs);
  prs$index <- ind;

  ind0 <- prs[prs[,1]==prs[,2] & prs[,1]!='-', 3];
  ind1 <- prs[prs[,1]!='-' & prs[,2]!='-' & prs[,1]!=prs[,2], 3];
  ind2 <- prs[prs[,2]=='-', 3];
  ind3 <- prs[which(prs[,1]=='-')-1, 3];
  ind3 <- ind3[!is.na(ind3)];
  ind4 <- Rle(as.integer(prs[,1]=='-'));
  ind4 <- runLength(ind4)[runValue(ind4)==1];

  cnt <- c(match=nmatch(pw), mismatch=nmismatch(pw), insertion=nindel(pw)@insertion[1,],
           deletion=nindel(pw)@deletion[1,]);
  idl <- nindel(pw);

  list(alignment=pw,
       index=list(match=ind0, mismatch=ind1, insertion=ind2, deletion=ind3, deletion.width=ind4),
       count=c(match=nmatch(pw), mismatch=nmismatch(pw), insertion=idl@insertion[1,], deletion=idl@deletion[1,]));
}


######## TODO
OptimizeConsensus <- function(freq, count, gap, ref) {
  CallCon <- function(x) {
    mx <- which(x==max(x));
    if (length(mx)>1) sample(mx, 1) else mx;
  };

  lvl1 <- unique(count);
}


