require(toolboxH)

venndiagram = function (x = x, y = y, z = z, w = w, unique = T, title = "Venn Diagram", 
          labels = c("x", "y", "z", "w"), lines = 1, lcol = 1, tcol = 1, 
          diacol = 1, plot = T, type = "3",fill_colors =fill_colors,  printsub = TRUE, ...) 
{
  if (unique == T) {
    x <- unique(x)
    x <- as.vector(stats::na.omit(x))
    y <- unique(y)
    y <- as.vector(stats::na.omit(y))
    if (!missing("z")) {
      z <- unique(z)
      z <- as.vector(stats::na.omit(z))
    }
    if (!missing("w")) {
      w <- unique(w)
      w <- as.vector(stats::na.omit(w))
    }
  }
  if (!type %in% c("2", "2map", "3", "3map", "4", "4map", 
                   "4el", "4elmap")) {
    return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
  }
  if (type == "2") {
    q1 <- x[x %in% y]
    q2 <- x[!x %in% y]
    q3 <- y[!y %in% x]
    qlist <- list(q1 = q1, q2 = q2, q3 = q3)
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query = names(count), count = as.vector(count))
    olDF <- data.frame(x = c(5, 3.1, 7), y = c(6.1, 6.1, 
                                               6.1), count = countDF$count)
    if (printsub == TRUE) {
      mysub <- paste(paste("N unique: xy =", length(unique(c(x, 
                                                             y)))), paste("; x =", length(unique(x))), paste("; y =", 
                                                                                                             length(unique(y))), sep = "")
    }
    else {
      mysub <- ""
    }
    if (plot == T) {
      # graphics::symbols(x = c(4, 6), y = c(6, 6), circles = c(2, 
      #                                                         2), xlim = c(0, 10), ylim = c(0, 10), inches = F, 
      #                   main = title, sub = mysub, xlab = "", ylab = "", 
      #                   xaxt = "n", yaxt = "n", bty = "n", fg = lines, 
      #                   ...)
      graphics::symbols(x = c(4, 6), y = c(6, 6), circles = c(2, 2), 
                        xlim = c(0, 10), ylim = c(0, 10), inches = F, 
                        main = title, sub = mysub, xlab = "", ylab = "", 
                        xaxt = "n", yaxt = "n", bty = "n", fg = lines, 
                        bg = if(!is.null(fill_colors)) fill_colors[1:2] else NA,
                          ...)
      
      graphics::text(olDF$x, olDF$y, olDF$count, col = tcol, 
                     ...)
      graphics::text(c(2, 8), c(8.8, 8.8), labels[1:2], 
                     col = lcol, ...)
    }
    return(qlist)
  }
  if (type == "2map") {
    olDFdebug <- data.frame(x = c(5, 3.1, 7), y = c(6.1, 
                                                    6.1, 6.1), count = paste("q", 1:3, sep = ""), ...)
    graphics::symbols(x = c(4, 6), y = c(6, 6), circles = c(2, 
                                                            2), xlim = c(0, 10), ylim = c(0, 10), inches = F, 
                      main = "Mapping Venn Diagram", xlab = "", ylab = "", 
                      xaxt = "n", yaxt = "n", bty = "n", fg = lines, ...)
    graphics::text(olDFdebug$x, olDFdebug$y, olDFdebug$count, 
                   col = tcol, ...)
    graphics::text(c(2, 8), c(8.8, 8.8), paste(labels[1:2], 
                                               "=", c("x", "y")), col = lcol, ...)
  }
  if (type == "3") {
    q1 <- x[x %in% y & x %in% z]
    q2 <- x[x %in% z]
    q2 <- q2[!q2 %in% y]
    q3 <- y[y %in% z]
    q3 <- q3[!q3 %in% x]
    q4 <- x[x %in% y]
    q4 <- q4[!q4 %in% z]
    q5 <- x[!x %in% y]
    q5 <- q5[!q5 %in% z]
    q6 <- y[!y %in% z]
    q6 <- q6[!q6 %in% x]
    q7 <- z[!z %in% x]
    q7 <- q7[!q7 %in% y]
    qlist <- list(q1 = q1, q2 = q2, q3 = q3, q4 = q4, q5 = q5, 
                  q6 = q6, q7 = q7)
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query = names(count), count = as.vector(count))
    olDF <- data.frame(x = c(5, 3.8, 6.3, 5, 3, 7, 5), y = c(5.6, 
                                                             4.6, 4.6, 6.9, 6.5, 6.5, 3), count = countDF$count)
    if (printsub == TRUE) {
      mysub <- paste(paste("N unique: xyz =", length(unique(c(x, 
                                                              y, z)))), paste("; x =", length(unique(x))), 
                     paste("; y =", length(unique(y))), paste("; z =", 
                                                              length(unique(z))), sep = "")
    }
    else {
      mysub <- ""
    }
    if (plot == T) {
      # Draw circles individually to allow different line colors
      graphics::symbols(x = 4, y = 6, circles = 2, xlim = c(0, 10), ylim = c(0, 10), 
                        inches = F, main = title, sub = mysub, xlab = "", ylab = "", 
                        xaxt = "n", yaxt = "n", bty = "n", 
                        fg = if(length(lines) >= 1) lines[1] else lines[1],
                        bg = if(!is.null(fill_colors) && length(fill_colors) >= 1) fill_colors[1] else NA,
                        add = FALSE, ...)
      graphics::symbols(x = 6, y = 6, circles = 2, inches = F,
                        fg = if(length(lines) >= 2) lines[2] else lines[1],
                        bg = if(!is.null(fill_colors) && length(fill_colors) >= 2) fill_colors[2] else NA,
                        add = TRUE, ...)
      graphics::symbols(x = 5, y = 4, circles = 2, inches = F,
                        fg = if(length(lines) >= 3) lines[3] else lines[1],
                        bg = if(!is.null(fill_colors) && length(fill_colors) >= 3) fill_colors[3] else NA,
                        add = TRUE, ...)
      
      # graphics::symbols(x = c(4, 6, 5), y = c(6, 6, 4), 
      #                   circles = c(2, 2, 2), xlim = c(0, 10), ylim = c(0, 
      #                                                                   10), inches = F, main = title, sub = mysub, 
      #                   xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
      #                   bty = "n", fg = lines, ...)
      # 
      graphics::text(olDF$x, olDF$y, olDF$count, col = tcol, 
                     ...)
      graphics::text(c(2, 8, 6), c(8.8, 8.8, 1.1), labels[1:3], 
                     col = lcol, ...)
    }
    return(qlist)
  }
  if (type == "3map") {
    olDFdebug <- data.frame(x = c(5, 3.8, 6.3, 5, 3, 7, 
                                  5), y = c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3), count = paste("q", 
                                                                                            1:7, sep = ""), ...)
    graphics::symbols(x = c(4, 6, 5), y = c(6, 6, 4), circles = c(2, 
                                                                  2, 2), xlim = c(0, 10), ylim = c(0, 10), inches = F, 
                      main = "Mapping Venn Diagram", xlab = "", ylab = "", 
                      xaxt = "n", yaxt = "n", bty = "n", fg = lines, ...)
    graphics::text(olDFdebug$x, olDFdebug$y, olDFdebug$count, 
                   col = tcol, ...)
    graphics::text(c(2, 8, 6), c(8.8, 8.8, 1.1), paste(labels[1:3], 
                                                       "=", c("x", "y", "z")), col = lcol, ...)
  }
  if (type == "4" | type == "4el" | type == "4elmap") {
    xy <- x[x %in% y]
    xz <- x[x %in% z]
    xw <- x[x %in% w]
    yz <- y[y %in% z]
    yw <- y[y %in% w]
    zw <- z[z %in% w]
    q1 <- xy[xy %in% zw]
    q2 <- xw[xw %in% z]
    q2 <- q2[!q2 %in% y]
    q3 <- yz[yz %in% w]
    q3 <- q3[!q3 %in% x]
    q4 <- yz[yz %in% x]
    q4 <- q4[!q4 %in% w]
    q5 <- xw[xw %in% y]
    q5 <- q5[!q5 %in% z]
    q6 <- xy[!xy %in% z]
    q6 <- q6[!q6 %in% w]
    q7 <- zw[!zw %in% x]
    q7 <- q7[!q7 %in% y]
    q8 <- xz[!xz %in% y]
    q8 <- q8[!q8 %in% w]
    q9 <- yw[!yw %in% x]
    q9 <- q9[!q9 %in% z]
    q10 <- x[!x %in% c(y, z, w)]
    q11 <- y[!y %in% c(x, z, w)]
    q12 <- z[!z %in% c(x, y, w)]
    q13 <- w[!w %in% c(x, y, z)]
    q14 <- xw[!xw %in% y]
    q14 <- q14[!q14 %in% z]
    q15 <- yz[!yz %in% x]
    q15 <- q15[!q15 %in% w]
    qlist <- list(q1 = q1, q2 = q2, q3 = q3, q4 = q4, q5 = q5, 
                  q6 = q6, q7 = q7, q8 = q8, q9 = q9, q10 = q10, q11 = q11, 
                  q12 = q12, q13 = q13, q14 = q14, q15 = q15)
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query = names(count), count = as.vector(count))
    olDF <- data.frame(x = c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 
                             4.8, 3, 6.5, 3, 6.5, 3, 6.5, 4.8, 4.8), y = c(5.2, 
                                                                           4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 
                                                                           3.2, 3.2, 1, 0.4), count = countDF$count)
    if (printsub == TRUE) {
      mysub <- paste(paste("N unique: xyzw =", length(unique(c(x, 
                                                               y, z, w)))), paste("; x =", length(unique(x))), 
                     paste("; y =", length(unique(y))), paste("; z =", 
                                                              length(unique(z))), paste("; w =", length(unique(w))), 
                     sep = "")
    }
    else {
      mysub <- ""
    }
    if (plot == T & type == "4") {
      graphics::symbols(x = c(4, 5.5, 4, 5.5), y = c(6, 
                                                     6, 4.5, 4.5), circles = c(2, 2, 2, 2), xlim = c(0, 
                                                                                                     10), ylim = c(0, 10), inches = F, main = title, 
                        sub = mysub, xlab = "", ylab = "", xaxt = "n", 
                        yaxt = "n", bty = "n", fg = lines, ...)
      graphics::text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], 
                     col = tcol, ...)
      graphics::text(c(2, 7.5, 2, 7.5), c(8.3, 8.3, 2, 
                                          2), labels, col = lcol, ...)
      graphics::text(c(3.8, 3.8), c(1, 0.4), c(paste("Only in ", 
                                                     labels[1], " & ", labels[4], ": ", olDF$count[14], 
                                                     sep = ""), paste("Only in ", labels[2], " & ", 
                                                                      labels[3], ": ", olDF$count[15], sep = "")), 
                     col = diacol, ...)
    }
    if (plot == T & (type == "4el" | type == "4elmap")) {
      olDF <- data.frame(x = c(5, 4.2, 6.4, 3.6, 5.8, 
                               2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5, 5), 
                         y = c(2.8, 1.4, 4, 4, 1.4, 5.9, 5.9, 2.2, 2.2, 
                               4.8, 7.2, 7.2, 4.8, 0.7, 6), count = countDF$count)
      plotellipse <- function(center = c(1, 1), radius = c(1, 
                                                           2), rotate = 1, segments = 360, xlab = "", ylab = "", 
                              ...) {
        angles <- (0:segments) * 2 * pi/segments
        rotate <- rotate * pi/180
        ellipse <- cbind(radius[1] * cos(angles), radius[2] * 
                           sin(angles))
        ellipse <- cbind(ellipse[, 1] * cos(rotate) + 
                           ellipse[, 2] * sin(rotate), ellipse[, 2] * 
                           cos(rotate) - ellipse[, 1] * sin(rotate))
        ellipse <- cbind(center[1] + ellipse[, 1], center[2] + 
                           ellipse[, 2])
        graphics::plot(ellipse, type = "l", xlim = c(0, 
                                                     10), ylim = c(0, 10), xlab = "", ylab = "", 
                       ...)
      }
      ellipseVenn <- function(lines = lines, olDF, title = title, 
                              labels = labels, sub = mysub, main, lcol = lcol, 
                              tcex = 1.3, ...) {
        graphics::split.screen(c(1, 1))
        plotellipse(center = c(3.5, 3.6), radius = c(2, 
                                                     4), rotate = -35, segments = 360, xlab = "", 
                    ylab = "", col = lines[1], axes = FALSE, main = title, 
                    sub = mysub, ...)
        graphics::screen(1, new = FALSE)
        plotellipse(center = c(4.7, 4.4), radius = c(2, 
                                                     4), rotate = -35, segments = 360, xlab = "", 
                    ylab = "", col = lines[2], axes = FALSE, ...)
        graphics::screen(1, new = FALSE)
        plotellipse(center = c(5.3, 4.4), radius = c(2, 
                                                     4), rotate = 35, segments = 360, xlab = "", 
                    ylab = "", col = lines[3], axes = FALSE, ...)
        graphics::screen(1, new = FALSE)
        plotellipse(center = c(6.5, 3.6), radius = c(2, 
                                                     4), rotate = 35, segments = 360, xlab = "", 
                    ylab = "", col = lines[4], axes = FALSE, ...)
        graphics::text(olDF[1:15, 1], olDF[1:15, 2], 
                       olDF[1:15, 3], col = tcol, ...)
        graphics::text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 
                                                8.3, 8.3, 7.3), labels, col = lcol, ...)
        graphics::close.screen(all = TRUE)
      }
      if (type == "4el") {
        ellipseVenn(olDF = olDF, lcol = lcol, lines = lines, 
                    labels = labels, title = title, ...)
      }
      if (type == "4elmap") {
        olDFdebug <- data.frame(x = c(5, 4.2, 6.4, 3.6, 
                                      5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 
                                      5, 5), y = c(2.8, 1.4, 4, 4, 1.4, 5.9, 5.9, 
                                                   2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6), count = paste("q", 
                                                                                                        1:15, sep = ""), ...)
        ellipseVenn(olDF = olDFdebug, lcol = lcol, lines = lines, 
                    labels = paste(labels, "=", c("x", "y", "z", 
                                                  "w")), title = "Mapping Venn Diagram", ...)
      }
    }
    return(qlist)
  }
  if (type == "4map") {
    olDFdebug <- data.frame(x = c(4.8, 3.9, 5.7, 3.9, 5.7, 
                                  4.8, 4.8, 3, 6.5, 3, 6.5, 3, 6.5, 4.8, 4.8), y = c(5.2, 
                                                                                     4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 
                                                                                     3.2, 3.2, 1, 0.4), count = paste("q", 1:15, sep = ""), 
                            ...)
    graphics::symbols(x = c(4, 5.5, 4, 5.5), y = c(6, 6, 
                                                   4.5, 4.5), circles = c(2, 2, 2, 2), xlim = c(0, 
                                                                                                10), ylim = c(0, 10), inches = F, main = "Mapping Venn Diagram", 
                      xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", 
                      fg = lines, ...)
    graphics::text(olDFdebug$x[1:13], olDFdebug$y[1:13], 
                   olDFdebug$count[1:13], col = tcol, ...)
    graphics::text(c(2, 7.5, 2, 7.5), c(8.3, 8.3, 2, 2), 
                   paste(labels, "=", c("x", "y", "z", "w")), col = lcol, 
                   ...)
    graphics::text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", 
                                                       labels[1], " & ", labels[4], ": ", olDFdebug$count[14], 
                                                       sep = ""), paste("Only in ", labels[2], " & ", labels[3], 
                                                                        ": ", olDFdebug$count[15], sep = "")), col = tcol, 
                   ...)
  }
}


venn2 = function (x1, y1, mytitle = "2-Way Venn Diagram", mylabels = NA, 
          plotte = T, venntype = "2", 
          fill_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 0.5),
          line_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 0),
          text_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 1)
          ) 
{
  if (all(is.vector(x1) | is.factor(x1), is.vector(y1) | is.factor(y1)) == 
      F) 
    stop("All input data must be vectors...")
  if (is.na(mylabels[1])) 
    mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)))
  qlist <- venndiagram(x = x1, y = y1, unique = T, title = mytitle, 
                       labels = mylabels, plot = plotte, 
                       lines = line_colors, 
                       lcol = text_colors, 
                       fill_colors = fill_colors,
                       tcol = c(1, 1, 1), lwd = 3, cex = 1.3, printsub = T, 
                       type = venntype)
  qlist
}


venn3 = function (x1, y1, z1, mytitle = "3-Way Venn Diagram", mylabels = NA, 
          plotte = T, venntype = "3",
          fill_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 0.5),
          line_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 0),
          text_colors = alpha(colour = c("steelblue", "salmon", "grey33"),alpha = 1)) 
{
  if (all(is.vector(x1) | is.factor(x1), is.vector(y1) | is.factor(y1), 
          is.vector(z1) | is.factor(z1)) == F) 
    stop("All input data must be vectors...")
  if (is.na(mylabels[1])) 
    mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)), 
                 deparse(substitute(z1)))
  qlist <- venndiagram(x = x1, y = y1, z = z1, unique = T, 
                       title = mytitle, labels = mylabels, plot = plotte, 
                       lines = line_colors, 
                       lcol = text_colors, 
                       fill_colors = fill_colors,
                       tcol = c(1, 1, 1, 1, 1, 1, 1), 
                       lwd = 3, cex = 1.3, printsub = T, type = venntype)
  qlist
}

# minimal examples ----
# customize colors with 
a=1:4
b=3:5
venn2(a,b)

mycolors2 = c("brown", "orange", "pink")

venn2(a,b, fill_colors = alpha(mycolors2, 0.2), line_colors = mycolors2, text_colors = mycolors2)

c=4:6
venn3(a,b,c)
mycolors3 = c("brown", "orange", "pink")

venn3(a,b,c, fill_colors = alpha(mycolors3, 0.3), line_colors = mycolors3, text_colors = mycolors3)

