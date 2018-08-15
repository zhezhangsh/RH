# Trim a long read into a set of shorter reads
TrimLongRead <- function(id, seq, qual, length, step=length, min.length=-Inf, output='fastq', filename='trimmed', thread=1) {
  # id      A character vector of read ids
  # seq     A character vector of read sequences, same length as <id>
  # qual    A character vector of read quality scores, same length as <id> and each element has matching length to <seq>
  # length  Length of trimmed reads
  # step    Size of stepwise trimming, to create overlapping trimmed reads
  # output  Output format
  # thread  Number of threads for parallele computing
  
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
    
    ii <- paste(i, z[, 1], sep='::');
    ss <- apply(z, 1, function(z) substr(s, z[1], z[2]));
    qq <- apply(z, 1, function(z) substr(q, z[1], z[2]));
    oo <- data.frame(from=z[, 1], to=z[, 2], seq=ss, qual=qq, stringsAsFactors = FALSE);
    rownames(oo) <- ii;
    
    oo;
  }; 
  
  if (thread <= 1) {
    trimmed <- lapply(1:length(id), function(i) trimLongRead(list(id[i], seq[i], qual[i], length, size)));
  } else {
    x <- lapply(1:length(id), function(i) list(id[i], seq[i], qual[i], length, size));
    
    require(snow);
    cl <- makeCluster(thread, type='SOCK');
    trimmed <- clusterApplyLB(cl, x, trimLongRead); 
    stopCluster(cl);
  }; 
  
  trimmed <- do.call('rbind', trimmed);
  trimmed <- trimmed[trimmed$to-trimmed$from+1>=min.length, , drop=FALSE];
  
  if (tolower(output[1])=='fastq') {
    fnm <- paste(filename, '.fastq', sep='');
    lns <- rbind(paste('@', rownames(trimmed), sep=''), trimmed$seq, rep('+', nrow(trimmed)), trimmed$qual);
    lns <- as.vector(lns);
    writeLines(lns, fnm);
  } else {}
}