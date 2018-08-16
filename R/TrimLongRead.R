# Trim a long read into a set of shorter reads
TrimLongRead <- function(id, seq, qual, length, step=length, min.length=-Inf, output='fastq', filename='trimmed', thread=1) {
  # id      A character vector of read ids
  # seq     A character vector of read sequences, same length as <id>
  # qual    A character vector of read quality scores, same length as <id> and each element has matching length to <seq>
  # length  Length of trimmed reads
  # step    Size of stepwise trimming, to create overlapping trimmed reads
  # output  Output format
  # thread  Number of threads for parallele computing

  if (tolower(output[1])=='fastq') fnm <- paste(filename, '.fastq', sep='') else fnm <- '';
  if (file.exists(fnm)) file.remove(fnm);

  newID <- function(i, loc) {
    w <- rev(gregexpr('/', i)[[1]])[1];
    a <- substr(i, 1, w-1);
    b <- strsplit(substr(i, w+1, nchar(i)), '_')[[1]];
    c <- as.integer(b[1]);
    d <- c + loc[, 1] - 1;
    e <- c + loc[, 2];
    paste(a, '/', d, '_', e, sep='');
  }

  trimLongRead <- function(x) {
    i <- x[[1]];
    s <- x[[2]];
    q <- x[[3]];
    l <- x[[4]];
    p <- x[[5]];
    n <- nchar(s);

    x <- seq(1, n, p);
    y <- seq(l, n+l, p);
    z <- cbind(from=x[1:length(y)], to=y);
    z[z>n] <- n;

    ii <- newID(i, z);
    ss <- apply(z, 1, function(z) substr(s, z[1], z[2]));
    qq <- apply(z, 1, function(z) substr(q, z[1], z[2]));
    oo <- data.frame(from=z[, 1], to=z[, 2], seq=ss, qual=qq, stringsAsFactors = FALSE);
    rownames(oo) <- ii;

    oo;
  };
  writeFastq <- function(fnm, trimmed, append=FALSE) {
    lns <- rbind(paste('@', rownames(trimmed), sep=''), trimmed$seq, rep('+', nrow(trimmed)), trimmed$qual);
    lns <- as.vector(lns);
    if (append) write(lns, fnm, append = TRUE) else writeLines(lns, fnm);
  }
  if (thread <= 1) {
    trimmed <- lapply(1:length(id), function(i) {
      trimmed <- trimLongRead(list(id[i], seq[i], qual[i], length, step));
      trimmed <- trimmed[trimmed$to-trimmed$from+1>=min.length, , drop=FALSE];
      if (nrow(trimmed) > 0) {
        if (tolower(output[1])=='fastq') writeFastq(fnm, trimmed, TRUE);
      }
    });
  } else {
    x <- lapply(1:length(id), function(i) list(id[i], seq[i], qual[i], length, step));

    require(snow);
    cl <- makeCluster(thread, type='SOCK');
    trimmed <- clusterApplyLB(cl, x, trimLongRead);
    stopCluster(cl);

    trimmed <- do.call('rbind', trimmed);
    trimmed <- trimmed[trimmed$to-trimmed$from+1>=min.length, , drop=FALSE];

    if (tolower(output[1])=='fastq') writeFastq(fnm, trimmed);
  };
}
