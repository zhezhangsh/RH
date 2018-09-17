IdentifyPalindrome <- function(seq, size, step, rate, lmin, name='') {
  require(Biostrings);
  require(GenomicRanges);

  if (isSingleString(seq)) seq <- DNAString(seq);
  if (identical(NA, name) | name=='') name <- paste('Len', nchar(seq), sep='');

  ind <- seq(1, nchar(seq)-size, step);
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

  cov <- coverage(ttl)[[1]];
  cv1 <- coverage(ttl[strand(ttl)=='+'])[[1]];
  cv2 <- coverage(ttl[strand(ttl)=='-'])[[1]];

  # Segmentation
  smt <- as.vector(runmean(cov, k=2*size+1, endrule = 'constant'));
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

  vld <- ValidatePalindrome(seq, loc[, 1], loc[, 2], 'all');

  out <- list(sequence=seq, alignment=aln, coverage=cov, location=loc, validation=vld$summary, trimmed=vld$trimmed);
  out;
};

ValidatePalindrome <- function(seq, stt, end, pair=c('all', 'next')) {
  require(Biostrings);

  BestAlign <- function(s1, s2) {
    pw1 <- pairwiseAlignment(s1, s2);
    pw2 <- pairwiseAlignment(s1, reverseComplement(s2));
    if (score(pw1) > score(pw2)) {
      pw <- pw1;
      di <- 1;
    } else {
      pw <- pw2;
      di <- -1;
    };

    al1 <- strsplit(as.character(subject(pw)), '')[[1]];
    al2 <- strsplit(as.character(pattern(pw)), '')[[1]];
    rle <- Rle(as.integer(al1==al2));
    val <- runValue(rle);
    len <- runLength(rle);
    stt <- start(rle);
    end <- end(rle);
    ind <- which(val>0 & len>3);

    stt1 <- al1[1:stt[ind[1]]]
    stt2 <- al2[1:stt[ind[1]]]
    stt1 <- length(stt1[stt1!='-']);
    stt2 <- length(stt2[stt2!='-']);

    end1 <- al1[end[ind[length(ind)]]:length(al1)];
    end2 <- al2[end[ind[length(ind)]]:length(al2)];
    end1 <- length(end1[end1!='-']);
    end2 <- length(end2[end2!='-']);
    end1 <- nchar(s1)-end1+1;
    end2 <- nchar(s2)-end2+1;

    sb1 <- gsub('-', '', subject(pw));
    sb2 <- gsub('-', '', pattern(pw));
    pw0 <- pairwiseAlignment(substr(sb1, stt1, end1), substr(sb2, stt2, end2));
    smm <- summary(pw0);

    c(score(pw0), di, stt1, end1, stt2, end2, smm@nmatch, smm@nmismatch, smm@ninsertion[2], smm@ndeletion[2])
  };

  ord <- order(stt);
  stt <- stt[ord];
  end <- end[ord];

  seg <- sapply(1:length(stt), function(i) subseq(seq, stt[i], end[i]));
  seg <- DNAStringSet(seg);
  names(seg) <- stt;

  if (tolower(pair[1])=='next') cmb <- cbind(1:(length(stt)-1), 2:length(stt)) else
    cmb <- t(combn(length(seg), 2));
  prs <- t(apply(cmb, 1, function(i) BestAlign(seg[[i[1]]], seg[[i[2]]])));
  tbl <- cbind(cmb, prs);
  colnames(tbl) <- c('Seg1', 'Seg2', 'Score', 'Strand', 'Start1', 'End1', 'Start2', 'End2',
                     'Match', 'Mismatch', 'Insertion', 'Deletion');

  scs <- split(rep(tbl[, 3], 2), c(tbl[, 1], tbl[, 2]));
  str <- split(rep(tbl[, 4], 2), c(tbl[, 1], tbl[, 2]));
  stt <- split(c(tbl[, 5], tbl[, 7]), c(tbl[, 1], tbl[, 2]));
  end <- split(c(tbl[, 6], tbl[, 8]), c(tbl[, 1], tbl[, 2]));
  max <- sapply(stt, max);
  min <- sapply(end, min);
  trm <- lapply(1:length(seg), function(i) subseq(seg[[i]], max[i], min[i]));
  trm <- DNAStringSet(trm);
  elementMetadata(trm) <- data.frame(name=names(seg), start=max, end=min, stringsAsFactors = FALSE);

  list(summary=tbl, trimmed=trm);
};
