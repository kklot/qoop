#' @importFrom grDevices adjustcolor dev.new palette
#' @importFrom graphics abline axis barplot boxplot layout legend lines 
#' @importFrom graphics mtext par plot points polygon rect segments text
showPch <- function() {
    dev.new()
    plot(0:25, pch=0:25, col=1:25)
    text(0:25, labels=0:25, pos=3, xpd=T)
    text(0:25, labels=0:25, pos=3, xpd=T)
}

solarized <- c("#002b36", "#dc322f", "#b58900", "#268bd2", "#859900", "#6c71c4", "#d33682", "#2aa198")
harvard <- c(cod.gray="#0b0b09", vivid.burgundy="#961b36", medium.champagne="#f6eeab", tulip.tree="#eba938", monte.carlo="#7dc4ba", red.damask="#d96043", wattle="#d4d849", wheat="#f7deb2", light.blue="#b2dbde", fern.frond="#64821c", swans.down="#d8e9dc", chocolate="#d46619", allports="#206b87", red.damask="#d76340", boston.blue="#4683a8", aluminum="#989897", silver="#c2c2c1")

#' @importFrom grDevices col2rgb rgb
AddAlpha <- function (plotclr, alpha = 0.5, verbose = 0) {
    tmp <- col2rgb(plotclr, alpha = alpha)
    tmp[4, ] = round(alpha * 255)
    for (i in 1:ncol(tmp)) {
        plotclr[i] = rgb(tmp[1, i], tmp[2, i], tmp[3, i], tmp[4, 
            i], maxColorValue = 255)
    }
    return(plotclr)
}

Kgrid <- function(bg = NA, cols = "gray93" ) {
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg,
         border = NA)
    xaxp <- par("xaxp"); yaxp <- par("yaxp")
    Vvec <- seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3])
    Hvec <- seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3])
    vvec <- seq(xaxp[1] + diff(Vvec)[1]/2, xaxp[2], by = abs(diff(Vvec)[1]))
    hvec <- seq(yaxp[1] + diff(Hvec)[1]/2, yaxp[2], by = abs(diff(Hvec)[1]))
    abline(v=Vvec, h = Hvec, lty=1, col = cols, lwd = 1)
    abline(v=vvec, h = hvec, lty=1, col = cols, lwd = .5)
}
vline <- function(v = 0, ...) abline(v = v, lty = 3, ...)
hline <- function(h = 0, ...) abline(h = h, lty = 3, ...)
put <- function(n.row, n.col, mar = NULL,...) {
  if (is.null(mar)) par(mfrow = c(n.row, n.col), mar = c(5,4,2,1)+0.1, ...)
  else par(mfrow = c(n.row, n.col), mar=mar, ...)
}
blankplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2) 
  }
}
lineplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2)
  }
  lines(...)
}
pointplot <- function(..., autoax = TRUE) {
  plot(..., type = "n", axes = FALSE)
  Kgrid()
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2)
  }
  points(...)
}
Kaxis <- function(side = 1, col='gray93', colticks='dimgray', ...) {
    axis(side, col=col, col.ticks=colticks, ...)
}

Kolygon <- function(x, y, ylow=NULL, col='gray70', alpha=0.4, border='white',...) {
  if (is.null(border)) border <- col
  if (missing(x)) x <- 1:length(y)
  if (missing(y)) {
    y <- x
    x <- seq(length(y))
  }
  xx <- c(x, rev(x))
  if (is.null(ylow)) yy <- c(y, rep(0, length(y)))
    else yy <- c(y, rev(ylow))
  polygon(xx, yy, col = AddAlpha(col, alpha), border = border, ...)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices  colorRampPalette
genCols <- function(n, pallete = "Set1") {
  getCols <- colorRampPalette(brewer.pal(8, pallete))
  return(getCols(n))
}

area_plot <- function(..., add=FALSE, autoax = TRUE) {
  if (!add) {
    plot(..., type = "n", axes = FALSE)
    Kgrid()
  }
  Kolygon(...)
  if (autoax) {
    Kaxis(1); Kaxis(2, las = 2) 
  }
}

# plot demographic.
# -----------------------------------------------------------------------------
plot.qoop <- function(mod, start_year=1971, min_age=0, bin_year=5, bin_age=5, max_age=91, colors=solarized, ...) {
  browser()
  mod <- mod$data[,,1,]
  if (interactive()) dev.new(width=9, height=6)
  put(1, 2, c(5,4,2,2))
  palette(colors)
  mains  <- c("Male", "Female")
     yl  <- dim(mod)[3]
  y_name <- start_year:(start_year+yl-1)
  n_age  <- dim(mod)[1]
  cols   <- rep(1:floor(n_age/5), each=bin_age, length.out = n_age)
  at_age <- which(!duplicated(cols))
  at_yr  <- which(!duplicated(1:yl%/%bin_year))

  sapply(1:dim(mod)[2], function(y) {
    blankplot(1:yl, numeric(yl), xlab="Years", ylab="Pop. size",
              ylim = c(0, max(mod)),
              main=mains[y], autoax=FALSE)
    Kaxis(1, at=at_yr, labels=y_name[at_yr])
    Kaxis(2)
    sapply(1:n_age, function(x) lines(mod[x, y, ], col=cols[x],...))
    text(yl+1, mod[at_age,y,yl], (min_age:max_age)[at_age], cex=.7)
  })
  invisible()
}