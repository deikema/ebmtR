##' Plot survival curves
##'
##' @title Plot survival curves
##' @param model Survival or prodlim object
survplot <- function(fit, cause = 1, byval = 12, tstart = 0, xlab = '', order = 0, stacked = F, title = '', colormap = 'Set1', bw = F, legendpos = 'topright') {
  
  if (length(cause) == 2) par(mfrow = c(1, 2))

  a <- 0.1
  xlab <- 'Months since Tx'
  tstart <- 0
  tmax <- max(fit$time)
  atrisk.at <- seq(tstart, tmax, byval)
  varnames <- gsub(paste0(fit$call$formula[3], '='), '', names(fit$strata))
  if (length(varnames) == 0) varnames <- ''
  leftmargin <- ifelse(length(varnames) != 1, ifelse(max(nchar(varnames)) >= 10, 10, 5), 5)
  bottommargin <- ifelse(length(fit$strata) == 0, 1, length(fit$strata))
  par(mar = c(1 + bottommargin,leftmargin, 2, 1), mgp = c(2.5, 0.75, 0))
  ngroups <- length(fit$strata)
  
  # replace lower and upper NA (only happens when surv == 0) by 0
  fit$lower[is.na(fit$lower)] <- 0
  fit$upper[is.na(fit$upper)] <- 0
  
  for (k in 1:length(cause)) {
    .cause <- ifelse(length(cause) == 2, k, cause) 
    
    cols <- RColorBrewer::brewer.pal(9, colormap)
    if (ngroups == 0 & length(cause) == 1 | bw == T) { 
      cols <- '#000000'
      linetype = 1:(ngroups + 1)
    }
    
    
    if (is.null(fit$strata)) fit$strata <- length(fit$time)
    ind <- c(0, cumsum(fit$strata))
    
    if (fit$type == 'right') {
      
      ylab <- 'Survival'
      plot(fit, conf.int = F, bty = 'l', axes = F, col = cols, main = title[k], lty = linetype)
      lower <- fit$lower
      upper <- fit$upper
      
      if (ngroups > 1) text(x = 0, y = 0.05, labels = test.logrank(fit), pos = 4)
      
    } else if (fit$type == 'mright' & !stacked) {
      ylab <- 'Cumulative incidence'
      lower <- fit$lower[, .cause]
      upper <- fit$upper[, .cause]
      # replace na nan etc by 0 because it just means they are out of plotting range
      lower[!sapply(lower, is.finite)] <- 0
      upper[!sapply(upper, is.finite)] <- 0
      # base plotting
      plot(fit[, .cause], conf.int = FALSE, bty = 'l', axes = FALSE, col = cols, ylim = c(0, 1), main = title[k])
      if (k == 1) fit$n.risk <- fit$n.risk[, ncol(fit$n.risk)]
      
      if (ngroups > 1) text(x = 0, y = 0.95, labels = test.gray(fit, .cause), pos = 4)
      
    } else if (fit$type == 'mright' & stacked) {
      ylab <- 'Cumulative incidence'
      
      inc <- fit$pstate[, 1:(ncol(fit$pstate) - 1)]
      
      if (order[1] != 0) {
        if (length(order) != (ncol(fit$pstate) - 1)) { order <- 1:(ncol(fit$pstate) - 1) }
        
        inc <- inc[, order]
        inc.cs <- t(apply(inc, 1, cumsum))
        plot(y = inc.cs[, 1], x = fit$time, bty = 'l', axes = FALSE, ylim = c(0, 1), 
             type = 's', col = cols[1], ylab = '', xlab = '', main = title[k], lty = linetype)
        inc.cs.int <- cbind(rep(0, length(fit$time)), inc.cs[, 1:length(order)])
        inc.cs.int <- rbind(rep(0, ncol(inc.cs.int)), inc.cs.int)
        
        for (j in 1:(length(order))) {
          u <- inc.cs.int[, j + 1]
          l <- inc.cs.int[, j]
          t <- c(0, fit$time)
          x <- rep(c(t, rev(t)), each = 2)
          y_u <- c(u[1], rep(u, each = 2)[1:((length(u)*2) - 1)])
          y_l <- c(rev(rep(l, each = 2))[2:(length(l)*2)], l[1])
          
          lines(y = inc.cs.int[, j + 1], x = t, type = 's', col = cols[j], lty = linetype[j])
          polygon(x, c(y_u, y_l), col = adjustcolor(cols[j], a), border = NA)
        }
        atrisk <- fit$n.risk[, ncol(fit$n.risk)]
        atrisk.n <- atrisk[c(1, findInterval(atrisk.at[-1], fit$time, all.inside = T) + 1)]
        mtext(text = atrisk.n, side = 1, line = 1, at = atrisk.at, cex = 0.9)
        
      }
    }
    
    axis(side = 2, at = seq(0, 1, by = 0.25), labels = paste0(seq(0, 100, 25), ' %') , las = 1, pos = 0)
    axis(side = 1, at = atrisk.at, xaxs = 'i', pos = 0)
    
    if (!stacked) {
      for (i in 1:ngroups) {
        
        x.grp <- fit$time[(ind[i] + 1):ind[i + 1]]
        
        yu.grp <- upper[(ind[i] + 1):ind[i + 1]]
        yl.grp <- lower[(ind[i] + 1):ind[i + 1]]
        x <- rep(c(x.grp, rev(x.grp)), each = 2)
        y_u <- c(yu.grp[1], rep(yu.grp, each = 2)[1:((length(yu.grp)*2) - 1)])
        y_l <- c(rev(rep(yl.grp, each = 2))[2:(length(yl.grp)*2)], yl.grp[1])
        polygon(x, c(y_u, y_l), col = adjustcolor(cols[i], a), border = NA)
        
        n.risk.grp <- fit$n.risk[(ind[i] + 1):ind[i + 1]]
        atrisk.grp <- n.risk.grp[c(1, findInterval(atrisk.at[-1], x.grp, all.inside = TRUE) + 1)]
        
        mtext(text = paste0(varnames[i], '      '), line = i, side = 1, at = -1, cex = 0.9, col = cols[i], adj = 1)
        mtext(text = atrisk.grp, side = 1, line = i, at = atrisk.at, cex = 0.9)
      }
    }
    mtext(text = ylab, at = 0.5, side = 2, line = 3, cex = 1.2)
    text(x = tmax, y = 0.03, labels = xlab, cex = 0.8, adj = 1)
  }

  if (bw) legend(legendpos, lty = linetype, legend = varnames)
  
  par(mfrow = c(1, 1))
}