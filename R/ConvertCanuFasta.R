ConvertCanuFasta <- function(fin, fout, name) {
  l <- readLines(fin);
  i <- grep('^>', l);
  h <- l[i];
  
  id <- sub('^>', '', sapply(strsplit(h, ' '), function(x) x[1]));
  ind <- cbind(i+1, c(i[-1]-1, length(l)));
  seq <- apply(ind, 1, function(i) paste(l[i[1]:i[2]], collapse=''));
  
  l1 <- paste0('>', name, '/', sub('^tig', '', id), '/', 0, '_', nchar(seq));
  
  writeLines(as.vector(rbind(l1, seq)), fout);
}