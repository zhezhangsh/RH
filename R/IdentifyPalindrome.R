IdentifyPalindrome <- function(seq, size, step, rate=0.75, lmin=1500, name=names(seq)) {
  require(Biostrings);
  require(GenomicRanges);

  if (isSingleString(seq)) seq <- DNAString(seq);
  if (identical(NA, name) | identical(NULL, name) | name=='') name <- paste('Len', nchar(seq), sep='');

  ind <- seq(1, length(seq)-size, step);
  sub <- DNAStringSet(lapply(ind, function(i) subseq(seq, i, width=size)));
  names(sub) <- ind;

  frq <- alphabetFrequency(sub);
  sub <- sub[frq[, 'N']/size <= (1-rate)];

  aln <- lapply(names(sub), function(nm) { print(nm);
    mp1 <- matchPattern(sub[[nm]], seq, max.mismatch = round(size*(1-rate)),
                        with.indels = TRUE);
    mp1 <- mp1[start(mp1) != as.integer(nm)];
    mp2 <- matchPattern(reverseComplement(sub[[nm]]), seq, max.mismatch = round(size*(1-rate)),
                        with.indels = TRUE);

    if (length(mp1) > 0) gr1 <- GRanges(name, IRanges(start(mp1), end(mp1)), strand='+') else
      gr1 <- GRanges();
    if (length(mp2) > 0) gr2 <- GRanges(name, IRanges(start(mp2), end(mp2)), strand='-') else
      gr2 <- GRanges();
    c(gr1, gr2);
  });
  names(aln) <- names(sub);

  # ol1 <- sapply(1:(length(aln)-1), function(i)
  #   sum(pmin(1, countOverlaps(aln[[i]], aln[[i+1]], maxgap = size*(1-rate)))));
  # ol2 <- sapply(1:(length(aln)-1), function(i)
  #   sum(pmin(1, countOverlaps(aln[[i+1]], aln[[i]], maxgap = size*(1-rate)))));

  ttl <- Reduce('c', aln);
  ttl$subread <- rep(names(sub), sapply(aln, length));

  cov <- coverage(ttl, width=nchar(seq))[[1]];

  # Segmentation
  smt <- as.vector(runmean(cov, k=2*round(lmin/4)+1, endrule = 'constant'));
  names(smt) <- 1:length(smt);
  inc <- rep(0, length(smt));
  for (i in 1:length(sub)) inc[as.numeric(names(sub)[i]):(as.numeric(names(sub)[i])+size)] <- 1;
  smt[inc==0] <- NA;
  ind <- which(inc>0);
  bmn <- sapply(ind, function(i) smt[i] == min(smt[max(1, i-lmin):min(length(smt), i+lmin)], na.rm=TRUE));
  rle <- Rle(0, length(smt));
  rle[ind[!bmn]] <- 1;
  loc <- cbind(start(rle), end(rle), runValue(rle));
  loc <- loc[loc[, 3]>0, ];
  loc <- loc[(loc[,2]-loc[,1]+1) >= lmin, , drop=FALSE];

  # vld <- SelectPalindromeSegment(seq, loc[, 1], loc[, 2], 'all');

  out <- list(sequence=seq, alignment=aln, coverage=cov, location=loc);
  out;
};

OptimizePairwiseAlignment <- function(s1, s2, window.size=50, min.similarity=0.8) {
  require(S4Vectors);
  require(Biostrings);

  window.size <- round(window.size);

  pw1 <- pairwiseAlignment(s1, s2);
  pw2 <- pairwiseAlignment(s1, reverseComplement(s2));
  if (score(pw1) > score(pw2)) {
    pw <- pw1;
    di <- 1;
  } else {
    pw <- pw2;
    di <- -1;
  };

  m1 <- strsplit(as.character(pattern(pw)), '')[[1]];
  m2 <- strsplit(as.character(subject(pw)), '')[[1]];

  agr <- Rle(as.integer(m1==m2));
  rnm <- runmean(agr, k=2*window.size+1, endrule = 'constant');
  ind <- which(rnm>=min.similarity);
  if (length(ind)==0) rng <- range(which(rnm==max(rnm))) else rng <- range(ind);

  t1 <- m1[rng[1]:rng[2]];
  t2 <- m2[rng[1]:rng[2]];
  t1 <- paste(t1[t1!='-'], collapse='');
  t2 <- paste(t2[t2!='-'], collapse='');

  stt1 <- m1[1:rng[1]];
  stt1 <- length(stt1[stt1!='-']);
  end1 <- stt1 + nchar(t1) - 1;

  stt2 <- m2[1:rng[1]];
  stt2 <- length(stt2[stt2!='-']);
  end2 <- stt2 + nchar(t2) - 1;

  pw0 <- pairwiseAlignment(t1, t2);
  smm <- summary(pw0);

  out <- c(score(pw0), rng[2]-rng[1]+1, di, stt1, end1, stt2, end2,
           smm@nmatch, smm@nmismatch, smm@ninsertion[2], smm@ndeletion[2]);
  names(out) <- c('Score', 'Length', 'Strand', 'Start1', 'End1', 'Start2', 'End2',
                     'Match', 'Mismatch', 'Insertion', 'Deletion');

  list(stat=out, alignment=pw0, similarity=rnm, trimmed=list(seq1=t1, seq2=t2));
};

SelectPalindromeSegment <- function(seq, stt, end, similarity=0.6) {
  require(Biostrings);
  require(msa);

  # Order segments by location
  ord <- order(stt);
  stt <- stt[ord];
  end <- end[ord];

  # Retrieve segment sequences from full sequences
  seg <- sapply(1:length(stt), function(i) subseq(seq, stt[i], end[i]));
  seg <- DNAStringSet(seg);
  names(seg) <- stt;

  # Best alignment of segment pairs
  cmb <- t(combn(length(seg), 2));
  opt <- t(apply(cmb, 1, function(i) OptimizePairwiseAlignment(seg[[i[1]]], seg[[i[2]]])));
  prs <- t(sapply(opt, function(x) x[[1]]));
  tbl <- cbind(Seg1=cmb[, 1], Seg2=cmb[, 2], prs);

  scs <- tbl[, 'Match'] - tbl[, 'Mismatch'] - tbl[, 'Insertion'] - tbl[, 'Deletion'];
  scs <- tbl[, 'Strand']*scs/max(scs);
  crr <- matrix(NA, nc=length(seg), nr=length(seg), dimnames=list(1:length(seg), 1:length(seg)));
  for (i in 1:nrow(tbl)) crr[tbl[i,1], tbl[i,2]] <- crr[tbl[i,2], tbl[i,1]] <- scs[i];

  mxr <- apply(crr, 1, function(r) max(abs(r), na.rm=TRUE));
  cff <- min(max(mxr), similarity);
  crr <- crr[mxr>=cff, mxr>=cff];
  mid <- apply(crr, 1, function(r) median(abs(r), na.rm=TRUE));
  while (length(mid[mid<=0.5]) & nrow(crr)>2) {
    crr <- crr[mid>min(mid), mid>min(mid)];
    mid <- apply(crr, 1, function(r) median(abs(r), na.rm=TRUE));
  };

  sgn <- crr;
  sgn[sgn<0] <- -1;
  sgn[sgn>0] <- 1;
  hcl <- hclust(as.dist(1-sgn));
  ctr <- cutree(hcl, h=0);
  cls <- lapply(1:max(ctr), function(i) rownames(crr)[ctr==i]);
  cls <- cls[order(sapply(cls, length))];

  sel <- seg[as.numeric(rownames(crr))];
  names(sel) <- rownames(crr);
  if (length(cls)>1) sel[cls[[2]]] <- reverseComplement(sel[cls[[2]]]);

  loc <- cbind(start=stt, end=end);
  rownames(loc) <- 1:nrow(loc);

  list(segment=sel, location=loc, pair=tbl, full=seq);
};

# BestAlign <- function(s1, s2) {
#   pw1 <- pairwiseAlignment(s1, s2);
#   pw2 <- pairwiseAlignment(s1, reverseComplement(s2));
#   if (score(pw1) > score(pw2)) {
#     pw <- pw1;
#     di <- 1;
#   } else {
#     pw <- pw2;
#     di <- -1;
#   };
#
#   al1 <- strsplit(as.character(subject(pw)), '')[[1]];
#   al2 <- strsplit(as.character(pattern(pw)), '')[[1]];
#   rle <- Rle(as.integer(al1==al2));
#   val <- runValue(rle);
#   len <- runLength(rle);
#   stt <- start(rle);
#   end <- end(rle);
#   ind <- which(val>0 & len>3);
#
#   stt1 <- al1[1:stt[ind[1]]]
#   stt2 <- al2[1:stt[ind[1]]]
#   stt1 <- length(stt1[stt1!='-']);
#   stt2 <- length(stt2[stt2!='-']);
#
#   end1 <- al1[end[ind[length(ind)]]:length(al1)];
#   end2 <- al2[end[ind[length(ind)]]:length(al2)];
#   end1 <- length(end1[end1!='-']);
#   end2 <- length(end2[end2!='-']);
#   end1 <- nchar(s1)-end1+1;
#   end2 <- nchar(s2)-end2+1;
#
#   sb1 <- gsub('-', '', subject(pw));
#   sb2 <- gsub('-', '', pattern(pw));
#   pw0 <- pairwiseAlignment(substr(sb1, stt1, end1), substr(sb2, stt2, end2));
#   smm <- summary(pw0);
#
#   c(score(pw0), di, stt1, end1, stt2, end2, smm@nmatch, smm@nmismatch, smm@ninsertion[2], smm@ndeletion[2])
# }
