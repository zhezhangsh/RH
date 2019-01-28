# From a set of subread alignment, identify a segment of ungapped alignment of continous subreads and the genomic loci they are aligned to.
# An alignment is a pair of subread and a location it was aligned to
IdentifyAlignmentSegment <- function(align, align.score, sub.start, sub.end, score.best=-1500, score.least=-1200, 
                                     min.number=10, min.length=1000, max.gap=0, ignore.strand=FALSE) {
  # align         GRanges object containing information about subread alignment
  # align.score   Vector of alignment scores, with matching length and order to parameter <align>
  # sub.start     Starting index of subread in the full read; must be a vector non-negative integer with matching length and order to parameter <align>
  # sub.end       Starting index of subread in the full read; must be a vector non-negative integer with matching length and order to parameter <align>
  # score.least   The minimum requirement of all alignment scores to keep subread alignments
  # score.best    The best alignment score require for a segment
  # min.number    Minimum number of alignments for a segment to be reported
  # min.length    Minimum length of bases for a segment to be reported
  # max.gap       Max length of gap allowed between subreads and the locations they aligned to; if 0, subreads and locations must be continous
  # ignore.strand If FALSE, segments will be identified separately for the 2 strands
  
  require(GenomicRanges);
  
  min.number  <- max(1, min.number[1]);
  min.length  <- max(1, min.length[1]);
  max.gap     <- max(0, max.gap[1]);
  align.score <- align.score[1:length(align)];
  sub.start   <- sub.start[1:length(align)];
  sub.end     <- sub.end[1:length(align)];
  
  align$as <- align.score;
  align$ss <- sub.start;
  align$se <- sub.end;
  
  chr <- as.vector(seqnames(align[1]));
  aln <- align[align$as <= score.least];
  if (is.null(names(aln))) names(aln) <- 1:length(aln);
  if (ignore.strand) aln <- list('*'=aln) else aln <- split(aln, as.vector(strand(aln)));
  aln <- aln[sapply(aln, length)>=min.number & sapply(aln, function(a) min(a$as, na.rm=TRUE))<=score.best];
  
  if (length(aln) == 0) NULL else {
    segs <- lapply(aln, function(a0) {
      a1 <- a0[order(start(a0))];
      end(a1) <- end(a1) + max.gap;
      cv <- coverage(a1)[[chr]];
      cv[cv>1] <- 1;
      gr <- GRanges(chr, IRanges(start(cv), end(cv)))[runValue(cv)==1];
      
      sg <- lapply(1:length(gr), function(i) { 
        a2 <- a1[countOverlaps(a1, gr[i])>0];
        a2 <- a2[order(a2$ss)];
        gr1 <- GRanges(chr, IRanges(a2$ss, a2$se+max.gap));
        names(gr1) <- names(a2);
        cv <- coverage(gr1)[[chr]];
        cv[cv>1] <- 1;
        rg <- GRanges(chr, IRanges(start(cv), end(cv)))[runValue(cv)==1];
        gr2 <- lapply(1:length(rg), function(i) {
          g <- gr1[countOverlaps(gr1, rg[i])>0];
          a2[names(g)];
        });
        gr2 <- gr2[sapply(gr2, function(g) min(g$as)) <= score.best];
        
        if (length(gr2) == 0) NULL else {
          stt1 <- sapply(gr2, function(g) min(start(g)));
          end1 <- sapply(gr2, function(g) max(end(g)));
          stt2 <- sapply(gr2, function(g) min(g$ss));
          end2 <- sapply(gr2, function(g) max(g$se));
          
          gr3 <- GRanges(chr, IRanges(stt1, end1)); 
          gr3$subread <- sapply(gr2, length);
          gr3$ss      <- stt2;
          gr3$se      <- end2;
          gr3$total   <- sapply(gr2, function(g) sum(g$as));
          gr3$mean    <- round(sapply(gr2, function(g) weighted.mean(g$as, g$se-g$ss+1)));
          gr3$best   <- sapply(gr2, function(g) min(g$as));
          
          gr3;
        }
      });
      sg <- sg[sapply(sg, length)>0];
      do.call('c', sg);
    });
    str <- rep(names(segs), sapply(segs, length));
    seg <- Reduce('c', segs); 
    if (length(seg)>0) {
      strand(seg) <- str;
      seg <- seg[order(start(seg))];
      seg <- seg[order(seg$ss)];
    };
    seg <- seg[width(seg)>=min.length & seg$subread>=min.number]
    if (length(seg) == 0) NULL else seg;
  }
}