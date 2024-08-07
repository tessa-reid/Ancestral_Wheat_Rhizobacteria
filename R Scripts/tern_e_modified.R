#tern_e function from: https://github.com/BulgarelliD-Lab/Microbiota_mapping/blob/main/QRMC-3HS_Fig1_SFig3/tern_e.R
#tern_e_modifed: modified for aesthetics including point shape and axis labels.

tern_e_modified <- function (x, scale = 1, dimnames = NULL, dimnames_position = c("corner", "edge", "none"), dimnames_color = "black", id = NULL, id_color = "black", coordinates = FALSE, grid = TRUE, grid_color = "gray", labels = c("inside", "outside", "none"), labels_color = "black", border = "black", bg = "white", pch = 19, cex = 1, prop_size = FALSE, shape = 19, col = "red", newpage = TRUE, pop = TRUE, ...) {
  labels <- match.arg(labels)
  if (grid == TRUE)
    grid <- "dotted"
  if (coordinates)
    id <- paste("(", round(x[, 1] * scale, 1), ",", round(x[, 2] * scale, 1), ",", round(x[, 3] * scale, 1), ")", sep = "")
  dimnames_position <- match.arg(dimnames_position)
  if (is.null(dimnames) && dimnames_position != "none")
    dimnames <- colnames(x)
  if (is.logical(prop_size) && prop_size)
    prop_size <- 3
  if (ncol(x) != 3)
    stop("Need a matrix with 3 columns")
  if (any(x < 0))
    stop("X must be non-negative")
  s <- rowSums(x)
  if (any(s <= 0))
    stop("each row of X must have a positive sum")
  x <- x/s
  top <- sqrt(3)/2
  if (newpage)
    grid.newpage()
  xlim <- c(-0.03, 1.03)
  ylim <- c(-1, top)
  pushViewport(viewport(width = unit(1, "snpc")))
  #if (!is.null(main))
    #grid.text(main, y = 0.9, gp = gpar(fontsize = 18, fontstyle = 1))
  pushViewport(viewport(width = 0.8, height = 0.8, xscale = xlim, yscale = ylim, name = "plot"))
  eps <- 0.01
  grid.polygon(c(0, 0.5, 1), c(0, top, 0), gp = gpar(fill = bg, col = border), ...)
  if (dimnames_position == "corner") {
    grid.text(x = c(0, 1, 0.5), y = c(-0.04, -0.04, top + 0.04), label = dimnames,gp = gpar(fontsize = 12))
  }
  if (dimnames_position == "edge") {
    shift <- eps * if (labels == "outside")
      8
    else 0
    grid.text(x = 0.25 - 2 * eps - shift, y = 0.5 * top + shift, label = dimnames[2], rot = 60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.75 + 3 * eps + shift, y = 0.5 * top + shift, label = dimnames[1], rot = -60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.5, y = -0.02 - shift, label = dimnames[3], gp = gpar(col = dimnames_color))
  }
  if (is.character(grid))
    for (i in 1:4 * 0.2) {
      grid.lines(c(1 - i, (1 - i)/2), c(0, 1 - i) * top, gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(1 - i, 1 - i + i/2), c(0, i) * top, gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(i/2, 1 - i + i/2), c(i, i) * top, gp = gpar(lty = grid, col = grid_color))
      if (labels == "inside") {
        grid.text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 * top, label = i * scale, gp = gpar(col = labels_color), rot = 120)
        grid.text(x = 1 - i + i/4 + eps, y = i/2 * top - eps, label = (1 - i) * scale, gp = gpar(col = labels_color), rot = -120)
        grid.text(x = 0.5, y = i * top + eps, label = i * scale, gp = gpar(col = labels_color))
      }
      if (labels == "outside") {
        grid.text(x = (1 - i)/2 - 6 * eps, y = (1 - i) * top, label = (1 - i) * scale, rot = -40, gp = gpar(col = labels_color))
        grid.text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 - i) * top + 5 * eps, label = i * scale, rot = 60, gp = gpar(col = labels_color))
        grid.text(x = i + eps, y = -0.05, label = (1 - i) * scale, rot = 60, gp = gpar(col = labels_color))
      }
    }
  xp <- x[, 2] + x[, 3]/2
  yp <- x[, 3] * top
  size = unit(if (prop_size)
    prop_size * rowSums(x*x*s) / max( rowSums(x*x*s) )
    #prop_size * sqrt(s/max(s))
    else cex, "lines")
  grid.points(xp, yp, pch = shape, gp = gpar(col = col), default.units = "snpc", 
              size = size, ...)
  if (!is.null(id))
    grid.text(x = xp, y = unit(yp - 0.015, "snpc") - 0.5 * 
                size, label = as.character(id), gp = gpar(col = id_color, 
                                                          cex = cex))
  if (pop)
    popViewport(2)
  else upViewport(2)
}
