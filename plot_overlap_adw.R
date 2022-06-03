
plot_overlap_adw <- function (object, niches = c(1, 2), niche_col = c("blue", 
                                                  "red"), data = TRUE, data_col = c("blue", "red"), 
          background = FALSE, background_type, proportion = 0.3, background_col = viridis::viridis, 
          change_labels = FALSE, xlab = "", ylab = "", 
          zlab = "", legend = TRUE) 
{
  if (missing(object)) {
    stop("Argument 'object' is necessary to perform the analysis.")
  }
  if (length(niche_col) < length(object@ellipsoids)) {
    message("Number of niches to plot exceeds number of colors, using automatic selection.")
    data_col <- rainbow(length(object@ellipsoids))
    niche_col <- rainbow(length(object@ellipsoids))
  }
  if (background == TRUE & missing(background_type)) {
    stop("Argument 'background_type' needs to be defined if background = TRUE.")
  }
  var_names <- object@variable_names
  iter <- niches[1]:niches[2]
  backs <- paste0("Niche_", niches[1], "_vs_", 
                  niches[2])
  if (length(var_names) > 2) {
    if (data == TRUE) {
      points <- lapply(iter, function(x) {
        sp_data <- object@data[[x]][, var_names]
        if (x == niches[1]) {
          if (change_labels == TRUE) {
            rgl::plot3d(sp_data[, 1:3], col = data_col[x], 
                        size = 8, xlab = xlab, ylab = ylab, zlab = zlab)
          }
          else {
            rgl::plot3d(sp_data[, 1:3], col = data_col[x], 
                        size = 8)
          }
        }
        else {
          rgl::plot3d(sp_data[, 1:3], col = data_col[x], 
                      size = 8, add = TRUE)
        }
      })
    }
    if (background == TRUE) {
      if (background_type == "full") {
        mh_sort <- object@full_background[[backs]]
      }
      if (background_type == "back_union") {
        mh_sort <- object@union_background[[backs]]
      }
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 
      ]
      mh_sort <- mh_sort[sample(ceiling(nrow(mh_sort) * 
                                          proportion)), ]
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 
                         1:3]
      if (data == TRUE) {
        rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)), 
                    add = TRUE)
      }
      else {
        if (change_labels == TRUE) {
          rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)), 
                      xlab = xlab, ylab = ylab, zlab = zlab)
        }
        else {
          rgl::plot3d(mh_sort, col = background_col(nrow(mh_sort)))
        }
      }
    }
    ellipsoides <- lapply(iter, function(x) {
      centroid <- object@ellipsoids[[x]]@centroid[1:3]
      cov_mat <- object@ellipsoids[[x]]@covariance_matrix[1:3, 
                                                          1:3]
      level <- object@ellipsoids[[x]]@level/100
      ell <- rgl::ellipse3d(cov_mat, centre = centroid, 
                            level = level)
      if (change_labels == TRUE) {
        rgl::wire3d(ell, col = niche_col[x], alpha = 1, 
                    xlab = xlab, ylab = ylab, zlab = zlab, lit = F)
      }
      else {
        rgl::wire3d(ell, col = niche_col[x], alpha = 1, lit = F)
      }
    })
    if (legend == TRUE) {
      rgl::legend3d("topright", legend = paste("Niche", 
                                               niches[1:2]), lty = 1, col = niche_col, inset = 0.02, 
                    bty = "n")
    }
  }
  else {
    el1 <- lapply(iter, function(x) {
      centroid <- object@ellipsoids[[x]]@centroid[1:2]
      cov_mat <- object@ellipsoids[[x]]@covariance_matrix[1:2, 
                                                          1:2]
      level <- object@ellipsoids[[x]]@level/100
      ellipse::ellipse(x = cov_mat, centre = centroid, 
                       level = level)
    })
    xlim <- range(unlist(lapply(el1, function(x) {
      x[, 1]
    })))
    ylim <- range(unlist(lapply(el1, function(x) {
      x[, 2]
    })))
    par(mar = c(4, 4, 1, 1))
    if (background == TRUE) {
      if (background_type == "full") {
        mh_sort <- object@full_background[[backs]]
      }
      if (background_type == "back_union") {
        mh_sort <- object@union_background[[backs]]
      }
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 
      ]
      mh_sort <- mh_sort[sample(ceiling(nrow(mh_sort) * 
                                          proportion)), ]
      mh_sort <- mh_sort[order(mh_sort[, "Niche_1_S"]), 
                         1:2]
      if (change_labels == TRUE) {
        plot(mh_sort, col = background_col(nrow(mh_sort)), 
             xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
      }
      else {
        plot(mh_sort, col = background_col(nrow(mh_sort)), 
             xlim = xlim, ylim = ylim, xlab = var_names[1], 
             ylab = var_names[2])
      }
    }
    if (data == TRUE) {
      points <- lapply(iter, function(x) {
        sp_data <- object@data[[x]][, var_names]
        if (x == niches[1]) {
          if (background == TRUE) {
            points(sp_data[, 1:2], pch = 19, col = data_col[x])
          }
          else {
            if (change_labels == TRUE) {
              plot(sp_data[, 1:2], pch = 19, col = data_col[x], 
                   xlim = xlim, ylim = ylim, xlab = xlab, 
                   ylab = ylab)
            }
            else {
              plot(sp_data[, 1:2], pch = 19, col = data_col[x], 
                   xlim = xlim, ylim = ylim, xlab = var_names[1], 
                   ylab = var_names[2])
            }
          }
        }
        else {
          points(sp_data[, 1:2], pch = 19, col = data_col[x])
        }
      })
    }
    ellipsoides <- lapply(iter, function(x) {
      lines(el1[[x]], col = niche_col[x], lwd = 1.5)
    })
    if (legend == TRUE) {
      legend("topright", legend = paste("Niche", 
                                        niches[1:2]), lty = 1, col = niche_col, bty = "n", 
             horiz = TRUE)
    }
  }
}

  