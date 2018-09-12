ExtractSamField <- function(fnm, path='', cmn=c(1, 2, 3, 4, 5, 6, 14, 15, 17)) {
  
  field.name <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual", 
                  read_group="rg", num_pass="np", query_end="qe", query_start="qs", 
                  hole_number="zm", align_score="as", "nm"); 
  
  fld <- field.name[cmn];
  out <- paste(sub('.sam$', '', fnm, ignore.case = TRUE));
  
  lns <- sapply(cmn, function(c) {
    fc <- paste('-f', c, sep='');
    ff <- paste(out, '_', field.name[c], '.txt', sep='');
    if (!identical(path[1], NA) & path[1]!='') ff <- paste(sub('/$', '', path[1]), out, sep='/');
    paste("grep -v '^@'", fnm, "| cut", fc, ">",  ff);
  });
  
  lns;
};

