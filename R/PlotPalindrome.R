PlotPalindrome <- function(start0, end0, start1, end1, strand, score=NA,
                           gene='', read='', type=c('read', 'gene')) {
  
  
  if (tolower(type[1])=='read') {
    xl <- 1000*c(floor(min(start1)/1000), ceiling(max(end1/1000)));
    xr <- xl[2]-xl[1];

    if (length(score) != length(strand)) col <- c('#E74C3C', '#3498DB')[(3-strand)/2] else {
      s <- -score;
      x <- round(100*(score-min(score))/(max(score)-min(score)));
      col <- gplots::colorpanel(101, '#E74C3C', '#3498DB')[x+1];
    };

    par(mar=c(0,0,2,0));
    plot(0, type='n', ylim=c(-2, 4), xlim=c(xl[1]-0.1*xr, xl[2]+0.1*xr), axes=FALSE, xlab='', ylab='');
    text(xl, 0, pos=c(2, 4), label=format(xl, big.mark = ','));
    rect(xl[1], -0.1, xl[2], 0.1, col='grey', border='black')
    points(seq(xl[1]+1000, xl[2]-1000, 1000), rep(0, xr/1000-1), col='white', pch='|', cex=0.25);

    title(main=paste(paste('Gene', gene), paste('Read', read), sep=' on '), cex.main=1.5);

    for (i in 1:length(start0)) {
      if (strand[i]==1) {
        arrows(start1[i], .5, end1[i], .5, length=0.025, lwd=1.5, col=col[i]);
        text(start1[i]/2+end1[i]/2, 1, srt=80, pos=3, font=2, col=col[i], cex=0.4,
             label=paste(format(start0[i], big.mark=','), format(end0[i], big.mark=','), sep=' - '));
      } else {
        arrows(end1[i], -.5, start1[i], -.5, length=0.025, lwd=1.5, col=col[i]);
        text(end1[i]/2+start1[i]/2, -1, srt=280, pos=1, font=2, col=col[i], cex=0.4,
             label=paste(format(end0[i], big.mark=','), format(start0[i], big.mark=','), sep=' - '));
      }
    }

    if (length(score) == length(strand))
      legend('topright', bty='n', lty=1, horiz = TRUE, lwd=2,
             col=gplots::colorpanel(101, '#E74C3C', '#3498DB')[c(1, 101)], legend=round(range(score)));

  } else if (tolower(type[1])=='gene') {
    xl <- 1000*c(floor(min(start0)/1000), ceiling(max(end0/1000)));
    xr <- xl[2]-xl[1];
    if (length(score) != length(strand)) col <- c('#E74C3C', '#3498DB')[(3-strand)/2] else {
      s <- -score;
      x <- round(100*(score-min(score))/(max(score)-min(score)));
      col <- gplots::colorpanel(101, '#E74C3C', '#3498DB')[x+1];
    };

    par(mar=c(0,0,2,0));
    plot(0, type='n', ylim=c(-length(strand), length(strand)/20), xlim=c(xl[1]-0.1*xr, xl[2]+0.1*xr), axes=FALSE, xlab='', ylab='');
    text(xl, 0, pos=c(2, 4), label=format(xl, big.mark = ','));
    rect(xl[1], -length(strand)/200, xl[2], length(strand)/200, col='grey', border='black')
    points(seq(xl[1]+1000, xl[2]-1000, 1000), rep(0, xr/1000-1), col='white', pch='|', cex=.25);
    title(main=paste(paste('Read', read), paste('Gene', gene), sep=' on '), cex.main=1.5);

    for (i in 1:length(start0)) {
      if (strand[i]==1) {
        arrows(start0[i], -i, end0[i], -i, length=0.05, lwd=1.5, col=col[i]);
        text(start0[i]/2+end0[i]/2, -i+1/3, font=2, col=col[i], cex=0.4,
             label=paste(format(start1[i], big.mark=','), format(end1[i], big.mark=','), sep=' - '));
      } else {
        arrows(end0[i], -i, start0[i], -i, length=0.05, lwd=1.5, col=col[i]);
        text(end0[i]/2+start0[i]/2, -i+1/3, font=2, col=col[i], cex=0.4,
             label=paste(format(end1[i], big.mark=','), format(start1[i], big.mark=','), sep=' - '));
      }
    }

    if (length(score) == length(strand))
      legend('topright', bty='n', lty=1, horiz = TRUE, lwd=2,
             col=gplots::colorpanel(101, '#E74C3C', '#3498DB')[c(1, 101)], legend=round(range(score)));

  }
}
