MakeConsensus <- function(seq, ref=NA, method=c("Muscle", "ClustalW", "ClustalOmega"),
                          weight=c('r', 'rank', 'count')) {
  require(msa);
  require(Biostrings);

  runMsa <- function(seq, method) {
    if (tolower(method)[1]=='clustalw') {
      msaClustalW(seq, order='input');
    } else if (tolower(method)[1]=='clustalomega') {
      msaClustalOmega(seq, order='input');
    } else {
      msaMuscle(seq, order='input');
    };
  };

  # Run MSA and remove sequence(s) not agree with the others
  removed <- c();
  pss <- FALSE;
  while(!pss & length(seq)>2) {
    seq <- seq[!(names(seq) %in% removed)];

    # Run consensus algorithm
    msa <- runMsa(seq, method);

    # Base by base
    aln <- do.call('cbind', strsplit(as.character(msa), ''));
    rownames(aln) <- 1:nrow(aln);

    # Weigth sequences
    wgt <- WeightConsensus(aln, method='count');

    if (min(wgt[[2]])<(0.5*(length(seq)-1))) {
      removed <- c(removed, names(seq)[wgt[[2]]==min(wgt[[2]])]);
    } else pss <- TRUE;
  };

  # Remove sequences hurts the consensus by reducing number of bases
  pss <- FALSE;
  while(!pss & length(seq)>2) {
    # Base by base
    aln <- do.call('cbind', strsplit(as.character(msa), ''));
    rownames(aln) <- 1:nrow(aln);

    wgt <- WeightConsensus(aln, method=weight);

    # Call consensus
    rg1 <- apply(aln, 2, function(a) range(which(a!='-')));
    rg0 <- c(max(rg1[1, ]), min(rg1[2, ]));
    cll <- CallConsensus(wgt$weighted[rg0[1]:rg0[2], , drop=FALSE]);
    cll <- cll[cll!='-'];

    sm0 <- sum(apply(cbind(cll, aln[names(cll), ]), 1, function(x) length(x[x==x[1]]))-1);
    sm1 <- sapply(1:ncol(aln), function(i) {
      a <- aln[, -i];
      w <- WeightConsensus(a, method=weight);
      r <- apply(a, 2, function(a) range(which(a!='-')));
      r <- c(max(r[1, ]), min(r[2, ]));
      c <- CallConsensus(w$weighted[r[1]:r[2], , drop=FALSE]);
      c <- c[c!='-'];
      n <- apply(cbind(c, a[names(c), ]), 1, function(x) length(x[x==x[1]]))-1;
      sum(n);
    });

    if (max(sm1) < sm0) pss <- TRUE else {
      removed <- c(removed, names(seq)[sm1==max(sm1)][1]);
      seq <- seq[!(names(seq) %in% removed)];

      # Run consensus algorithm
      msa <- runMsa(seq, method);
    }
  };

  list(consensus=paste(cll, collapse=''), range=rg0, base=aln, alignment=msa, weight=wgt,
       removed=removed, input=list(sequence=seq, reference=ref, method=method, weight=weight));
};

WeightConsensus <- function(mtrx, method=c('r', 'rank', 'count')) {
  ltt <- c('-', 'A', 'C', 'G', 'T')
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

AlignConsensus2Reference2 <- function(cons, ref) {
  nms <- names(cons);
  if (is.null(nms) | length(nms)!=length(cons)) nms <- paste('seq', 1:length(cons));
  if (class(cons) != 'DNAStringSet') cons <- DNAStringSet(cons);
  names(cons) <- nms;

  aln <- lapply(names(cons), function(nm) {
    cat(nm, '\n');
    con <- cons[[nm]];
    pw1 <- AlignConsensus2Reference(con, ref);
    pw2 <- AlignConsensus2Reference(reverseComplement(con), ref);

    if (score(pw1$alignment) >= score(pw2$alignment)) {
      pw <- pw1;
      pw$strand <- 1;
    } else {
      pw <- pw2;
      pw$strand <- -1;
    }
    pw;
  });
  names(aln) <- nms;

  stat <- t(sapply(aln, function(x) x[[3]]));
  stat <- cbind(score=sapply(aln, function(a) score(a$alignment)), strand=sapply(aln, function(a) a$strand), stat);

  list(summary=stat, consensus=cons, alignment=lapply(aln, function(a) a$alignment),
       index=lapply(aln, function(a) a$index))
}

AlignConsensus2Reference <- function(con, ref) {
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


