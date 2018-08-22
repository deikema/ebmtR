# MEDIAN FOLLOWUP [checked]
med <- function(fit, digits = 1, cause = 1, timescale = 'months', caption = '', internal = F) {
     ngroups <- ifelse(is.null(fit$strata), 1, length(fit$strata))
     r1 <- rep(NA, ngroups)
     medfu.str <- rep(NA, ngroups)
     
     if (ngroups > 1) {
          subsets <- c(0, cumsum(fit$strata))
          r2 <- names(fit$strata)
     } else {
          subsets <- c(0, length(fit$time))
          r2 <- NULL
     }
     for (i in 1:ngroups) {
          if (fit$type == 'right') {
               est <- fit$surv[(subsets[i] + 1):subsets[i + 1]]
               lower <- fit$lower[(subsets[i] + 1):subsets[i + 1]]
               upper <- fit$upper[(subsets[i] + 1):subsets[i + 1]]
               
               lower[is.na(lower)] <- 0
               upper[is.na(upper)] <- 0
               
               time.tmp <- fit$time[(subsets[i] + 1):subsets[i + 1]]
               
               if (min(est) > 0.5) {
                    medfu.str[i] <- 'Not reached' 
               } else {
                    medfu.str[i] <- paste0(round(time.tmp[min(which(est < 0.5))], digits), ' (', 
                                           round(time.tmp[min(which(lower < 0.5))], digits), ' - ',
                                           round(time.tmp[min(which(upper < 0.5))], digits), ')')
               }
          } else {
               est <- fit$pstate[(subsets[i] + 1):subsets[i + 1], cause]
               lower <- fit$lower[(subsets[i] + 1):subsets[i + 1], cause]
               upper <- fit$upper[(subsets[i] + 1):subsets[i + 1], cause]
               
               lower[is.na(lower)] <- 0
               upper[is.na(upper)] <- 0
               
               time.tmp <- fit$time[(subsets[i] + 1):subsets[i + 1]]
               
               if (min(est) > 0.5) {
                    medfu.str[i] <- 'Not reached' 
               } else {
                    medfu.str[i] <- paste0(round(time.tmp[min(which(est > 0.5))], digits), ' (', 
                                           round(time.tmp[min(which(upper > 0.5))], digits), ' - ',
                                           round(time.tmp[min(which(lower > 0.5))], digits), ')')
               }
          }
     }
     r1[1] <- 'Median time (95% CI): '
     df <- data.frame(r1 = cbind(r1, r2), medfu = medfu.str, row.names = NULL)
     options(knitr.kable.NA = '')
     if (internal) {
          return(df$medfu)
     } else {
          return(knitr::kable(df, col.names = c(rep('', ncol(r1)), timescale), caption = caption))
     }
}

# TABLE FUNCTIONS
#EXAMINE [checked]
examine <- function(cont, cat = '', colname = '', caption = '', xcol = '', digits = 1) {
     if (cat[1] == '' & length(cat) == 1) cat <- factor(rep(colname, length(cont)))
     tb <- do.call(data.frame, aggregate(cont, list(cat), function(x) { round(quantile(x, na.rm = T), digits)}))
     n <- table(droplevels(cat[!is.na(cont)], useNA = 'no'))
     m <- rep(NA, nlevels(droplevels(cat)))
     tb <- cbind(n, m, tb[, 2:ncol(tb)])
     total <- c('Total', 
                length(cont[!is.na(cat)]), 
                sum(is.na(cont) | is.na(cat)), 
                round(quantile(cont, na.rm = T), digits ))
     tb[, 3] <- table(droplevels(cat[is.na(cont)], useNA = 'no'))
     if (nlevels(cat) > 1) tb <- suppressWarnings(rbind(total, tb))
     if (xcol != '') rownames(tb)[1] <- xcol
     cn <- c(colname, 'N', 'Missing', 'Min','Q1','Median','Q3','Max')
     return(knitr::kable(tb, caption = caption, digits = digits, col.names = cn))
}

# CONDITIONING: ADD REVISED TABLE


# FREQUENCY TABLE
freqtab <- function(x, caption = '', digits = 1, order = '', valid.only = F, percent = T, ...) {
     inp <- match.call()
     varargin <- list(...)
     labs <- ''
     unpack.list(varargin)
     
     vname <- unlist(strsplit(as.character(inp[2]), '\\$'))[2]
     caption <- ifelse(caption == '' & is.data.frame(labs), 
                       as.character(labs[which(rownames(labs) == vname), ]), caption)
     
     if (class(x) == 'factor') x <- droplevels(x)
     perc <- 'Percent'
     vperc <- 'Valid percent'
     
     if (!percent) perc <- vperc <- NULL
     if (valid.only) {
          percent <- TRUE
          vperc <- 'Valid percent'
     }
     freqs <- table(x, useNA = 'ifany')
     props.incl <- prop.table(freqs)*100
     props.excl <- prop.table(table(x, useNA = 'no'))*100
     if (sum(is.na(x)) > 0) props.excl <- c(props.excl, NA)
     rnames <- names(freqs)
     rnames[is.na(rnames)] <- 'Missing'
     rnames <- c(rnames, 'Total')
     tb <- data.frame(rbind(cbind(freqs, props.incl, props.excl),
                            c(sum(freqs), 100, 100 )), row.names = rnames)
     tb.missing <- NULL
     .omitted.rows <- 1
     if (order == 'desc') {
          if (rownames(tb)[length(rownames(tb)) - 1] == 'Missing') {
               tb.missing <- tb[dim(tb)[1] - 1, ]
               .omitted.rows <- 2
          }
          tb.total <- tb[dim(tb)[1], ]
          tb <- tb[order(tb$freqs[1:(dim(tb)[1] - .omitted.rows)], decreasing = T), ]
          tb <- rbind.data.frame(tb, tb.missing, tb.total)
     }
     tb <- cbind(rownames(tb), tb)
     rownames(tb) <- NULL
     if (!percent) {
          tb <- tb[, 1:2] 
          valid.only <- F
     }
     if (valid.only) tb <- tb[, c(1, 2, 4)]
     
     cnames <- c('', 'Frequency', perc, vperc)
     options(knitr.kable.NA = '')
     return(knitr::kable(tb, digits = digits, row.names = NA, col.names = cnames, caption = caption))
}

# stat tests (to use: test = TRUE) perhaps manual test selection would be better?
# cat variables
# chisq
# cont vars:
# 2 subgroups: mann whitney u
# >2 subgroups: kruskal wallis test

# cat.names for custom categorical names.
# cont.names for custom cont. names. 

# the thing gets confused when all variables have an equal number of categories. If this happens, it will but the columns side by side, whereas if at least 1 column has a different number of categories, this will not happen.
tableone <- function(data, by = '', margin = 2, cat.names = '', cont.names = '', caption = '', digits = 1, test = F, df = F) {
     isfactor <- sapply(data, is.factor)
     df.cat  <- data.frame(data[, isfactor])
     for (i in 1:ncol(df.cat)) df.cat[, i] <- droplevels(df.cat[, i])
     
     lvls <- lapply(data, levels)
     cont <- sapply(lvls, is.null)
     df.cnt <- data.frame(data[, cont])
     
     byvar <- factor(rep(1, nrow(df.cat)))
     
     cat.tab <- NULL
     if (ncol(df.cat) > 0) {
          # categorical stuff
          if (cat.names[1] == '') cat.names <- names(df.cat)
          vnames <- trimws(rep(cat.names, sapply(lvls[!cont], length)), 'right')
          gnames <- trimws(unlist(lvls[!cont]), 'right')
          
          if (by != '') {
               byvar <- data[, by]
               
               perc <- vals <- NULL
               df.cat$temp <- factor(1)
               
               vals <- do.call(rbind, apply(df.cat, 2, function(x) { table(x, byvar, useNA = 'no') }))
               perc <- do.call(rbind, apply(df.cat, 2, function(x) {
                    round(prop.table(table(x, byvar, useNA = 'no'), margin = margin)*100, digits = digits) }))
               
               mat <- matrix(NA, nrow = nrow(vals), ncol = ncol(vals)*2)
               mat[, seq(1, dim(vals)[2]*2, 2)] <- vals
               mat[, seq(2, dim(vals)[2]*2, 2)] <- perc
               miss <- sapply(df.cat[, -ncol(df.cat)], function(x) {sum(is.na(x))})
               
               if (test) {
                    cat.p <- sapply(df.cat, function(x) {
                      if (length(levels(x)) == 1) {
                        return(NA) 
                      } else {
                        return(round(chisq.test(x, y = byvar)$p.value, 3)) 
                      }
                    })
                    mat <- data.frame(mat, p = NA)
                    mat$p[!duplicated(vnames)] <- cat.p
               }
               
               mat <- mat[-nrow(mat), ]
          } else {
               test <- FALSE
               # add temp variable to get rid of the weird behavior when all factors have an equal number of levels
               df.cat$temp <- factor(1)
               miss <- sapply(df.cat[, -ncol(df.cat)], function(x) {sum(is.na(x))})
               vals <- unlist(sapply(df.cat, function(x) {table(x, useNA = 'no' ) }))
               
               perc <- unlist(sapply(df.cat, function(x) { round(prop.table(table(x, useNA = 'no'))*100, 
                                                                 digits = digits) }))
               mat <- data.frame(X1 = vals, X2 = perc)
               
               # get rid of temp.
               mat <- mat[-nrow(mat), ]
          }
          
          missing.vect <- rep(miss, sapply(lvls[!cont], length))
          dupnames <- duplicated(data.frame(vnames, missing.vect))
          
          vnames[dupnames] <- NA
          missing.vect[dupnames] <- NA
          missing.vect[missing.vect == 0] <- NA
          
          cat.tab <- data.frame(vnames, gnames, missing.vect, mat, row.names = NULL)
     }
     cnt.tab <- NULL
     if (ncol(df.cnt) > 0) {
          if (cont.names[1] == '') cont.names <- names(df.cnt)
          if (by != '') byvar <- data[, by] 
          
          cnt.tmp <- sapply(df.cnt, function(x) {
               tapply(x, byvar,  function(x) {
                    qt <- round(quantile(x, na.rm = TRUE, probs = c(0.5, 0, 1)), digits)
                    paste0(qt[1], ' (', qt[2], '-', qt[3], ')')})
          })
          
          if (by != '') { cnt.tmp <- data.frame(t(cnt.tmp)) 
          } else {cnt.tmp <- data.frame(cnt.tmp) }
          
          miss.cnt <- sapply(df.cnt, function(x) {sum(is.na(x))})
          miss.cnt[miss.cnt == 0] <- NA
          cnt.tmp2 <- data.frame(matrix(NA, nrow = nrow(cnt.tmp), ncol = ncol(cnt.tmp)*2) )
          
          freq.cnt <- sapply(df.cnt, function(x) {tapply(x, byvar, function(y) {sum(!is.na(y))})})
          
          if (by != '') freq.cnt <- t(freq.cnt)
          
          cnt.tmp2[, seq(1, ncol(cnt.tmp2), 2) ] <- freq.cnt
          cnt.tmp2[, seq(2, ncol(cnt.tmp2), 2) ] <- cnt.tmp
          
          cnt.tab <- data.frame(vnames = cont.names,
                                gnames = 'Median (range)',
                                missing.vect = miss.cnt,
                                cnt.tmp2,
                                row.names = NULL)
          if (test) {
               if (nlevels(byvar) == 2) {
                    cnt.p <- data.frame(p = sapply(df.cnt, function(x) { round(wilcox.test(x ~ byvar)$p.value, 3)}))
               } else {
                    cnt.p <- data.frame(p = sapply(df.cnt, function(x) { round(kruskal.test(sd$age_TX1 ~ byvar)$p.value, 3)}))
               }
               cnt.tab <- cbind(cnt.tab, cnt.p)
          }  
          table1 <- rbind(cat.tab, cnt.tab)
     } else {
          cont.names <- NULL
          table1 <- cat.tab
     }
     
     if (by != '') table1 <- table1[c(rep(cat.names, sapply(lvls[!cont], length)), cont.names) != by, ]
     col.names <- c('','', 'N Missing', rep(c('N', '%'), nlevels(byvar)))
     if (nlevels(byvar) > 1) {
       row1 <- rep(NA, ncol(table1))
       row1[col.names == 'N'] <- levels(byvar)
       table1 <- rbind(row1, table1)
     }
     if (test) { 
          col.names <- c(col.names, 'p')
          table1[1, ncol(table1)] <- NA
     }
     
     if (df) {
          return(table1)
     } else {
          options(knitr.kable.NA = '')
          return(knitr::kable(table1, row.names = FALSE, col.names = col.names, caption = caption))
     }
}


# FGR does not allow for interactions?
# so in the case of any interactions, manually compute the variable yourself beforehand...
# dont use FGR yet. issues with interactions


# TBD
# 1. variable group names for interactions
# 2. indicator for tt variables
# 3. does it work with long form data?

# COVARIATE TABLE FOR COXPH
# REPLACEMENT FOR COXOUT:-> INCLUDES SUPPORT FOR RISKREGRESSION (NOT CRR)
# THIS IS A DEV VERSION (WORKS IF YOU CAN IGNORE TERRIBLE CODING)
# ALSO, FOR CSH, BE SURE TO FIRST CODE THE STATUS VAR AS 0-1, FOR WHATEVER CAUSE YOU ARE USING, OTHERWISE IT DOES NOT WORK CORRECTLY (IT EXPECTS THE OUTCOME INDICATOR AS A VARIABLE NAME, IT CANT HANDLE A CONDITIONAL)
mvatab <- function(fit, caption = '', showevents = TRUE, digits = 2) {
     # get data frame used in model
     data <- get(as.character(fit$call$data))
     # get the outcome variables (differentiate between Hist or Surv object)
     outcome <- unlist(strsplit(gsub(pattern = paste0(substr(fit$call$formula[2], 1, 4), '\\(|\\)'), '', 
                                     as.character(fit$call$formula[2])), ', '))
     
     # specify the type of hazard ratio based on the entered model
     estimate <- ifelse(is.null(fit$cause), 'HR', 'sHR')
     
     # flag if dealing with rsadd model:
     relsurvobj <- ifelse(class(fit)[1] == 'rsadd', 1, 0) 
     
     # extract the covariate names
     tmp <- unlist(strsplit(as.character(fit$call$formula[3]), ' +'))
     tmp <- tmp[seq(1,length(tmp), by = 2)]
     
     # find tt if it exists and extract the variable name and store the index
     # for now just keep it and treat it as a contnouos variable
     # find index
     tt_index <- grep(pattern = 'tt\\(', tmp)
     # replace
     tmp <- gsub(pattern = '(tt\\()|\\)', '', tmp, fixed = F)
     # keep variable name
     tt_var <- tmp[tt_index]
     
     # subset data to only get included variable names
     data <- data[, c(tmp, outcome)]
     # and only use complete cases
     data <- data[complete.cases(data), ]
     
     # remove the ratetable and following items from covariate list
     if (relsurvobj) tmp <- tmp[1:(pmatch('ratetable', tmp) - 1)]
     # remove potential strata, cluster or frailty terms from table
     other.covs <- grep(pattern = '[(strata)|(cluster)|(frailty)]*([a-zA-Z]).+)', tmp)
     if (length(other.covs > 0)) tmp <- tmp[-other.covs]
     
     # compute hazard ratios
     tmp.tab <- summary(fit)
     # only use hr, confidence interval and p-values
     if (!relsurvobj) {
          tmp.tab2 <- cbind(round(tmp.tab$conf.int[, c(1, 3, 4)], 2), 
                            round(tmp.tab$coef[, ncol(tmp.tab$coef)], 3))
          tmp.tab$nevent <- sum(data[, outcome[2]])
     } else {
          tmp.tab2 <- cbind(round(exp(tmp.tab$coefficients[, 1]), digits),
                            round(exp(tmp.tab$coefficients[, 1] - 1.96*tmp.tab$coefficients[, 2]), digits), 
                            round(exp(tmp.tab$coefficients[, 1] + 1.96*tmp.tab$coefficients[, 2]), digits),
                            round(tmp.tab$coefficients[, 4], 3))
          # we also need the number of events and sample size for later
          tmp.tab$nevent <- fit$nevent
          tmp.tab$n <- fit$n
     }
     
     vnames.rep <- rep(names(fit$assign), lapply(fit$assign, function(x) length(x) + 1))
     # find which variables are continuous and which are categorical
     # categorical: 
     catvars <- !is.na(match(vnames.rep, names(fit$xlevels)))
     intvars <- grep(':', vnames.rep)
     
     # finally replace pvalues <0.001 by <0.001 string 
     tmp.tab2[tmp.tab2[, 4] < 0.001, 4] <- '<0.001'
     
     # paste final hazard ratios
     tb <- cbind(paste0(tmp.tab2[, 1], '(', tmp.tab2[, 2], '-', tmp.tab2[, 3], ')'), tmp.tab2[, 4])
     
     # start with table
     # get number of levels for each included covariate
     
     # create table of NA, we will fill in later
     mat <- matrix(NA, nrow = length(vnames.rep), ncol = 6)
     mat[ , 1] <- vnames.rep
     
     mat[duplicated(vnames.rep), c(5,6)] <- tb
     
     freqs <- unlist(sapply(data[, names(fit$xlevels)], function(x) table(x, useNA = 'no')))
     grp.events <- unlist(sapply(data[data[, outcome[2]] == 1, names(fit$xlevels)], function(x) table(x, useNA = 'no')))
     
     mat[catvars, 2] <- unlist(fit$xlevels)
     mat[catvars, 3] <- freqs
     mat[catvars, 4] <- grp.events
     
     # remove reference rows for either the interactions or cont variables
     mat <- mat[-which(!catvars & !duplicated(mat[, 1])), ]
     mat[duplicated(mat[, 1]), 1] <- NA
     
     mat <- rbind(c('Total', NA, tmp.tab$n, tmp.tab$nevent, NA, NA), mat)
     
     colnames.str <- c('','group','n','n events',paste0(estimate, ' (95% CI)'),'p')
     
     options(knitr.kable.NA = '')
     return(knitr::kable(mat, col.names = colnames.str, caption = caption, row.names = FALSE))
}


# below: add option

##' Summarize survival estimates from survival objects
##' 
##' @title Summarize survival estimates
##' @param fit Survfit object, works with single and competing risks outcomes. if a list is passed, survtabs for each covariate are stacked.
##' @param times Vector of survival times at which estimates are given
##' @param cause Competing risk cause. The default is 1. If the status indicator in the survfit object is a factor, the first level is assumed to be censored, and cause = 1 corresponds to the first non-censored status type, 2 for the second etc.  
##' @param digits Number of decimal points in estimates
##' @param caption Table title
##' @param outcomes Vector of outcome names, will be included in the final table. If no oucomes are given, the corresponding table cells will be empty. 
##' 
survtab <- function(fit, times, cause = 1, digits = 0, caption = '', outcomes = '', median = F, test = T) {
     tb.final = NULL
     modelNum <- 1
     if (class(fit) == 'list') { 
          fitlist <- fit
          modelNum <- length(fitlist)
     }
     for (k in 1:modelNum) {
          if (modelNum > 1) fit <- fitlist[[k]]
          
          
          med.tit <- NULL
          med.str <- NULL
          pname <- NULL
          gt <- NULL
          finaltab <- NULL
          
          if (length(outcomes) < length(cause)) outcomes <- rep('', length(cause))
          
          ngroups <- ifelse(is.null(fit$strata), 1, length(fit$strata))
          # add leading 0 for reasons
          if (times[1] != 0) times <- c(0, times)
          if (ngroups == 1) test <- F
          
          # force single cause if survival model, to avoid mistakes
          if (fit$type == 'right') cause <- 1
          for (i in 1:length(cause)) {
               # extract table (use extend)
               tb <- summary(fit, times = times, extend = T)
               # replace any na(n)'s by 0
               if (any(is.na(tb$upper))) tb$upper[is.na(tb$upper)] <- 0
               if (any(is.na(tb$lower))) tb$lower[is.na(tb$lower)] <- 0
               
               .cause <- ifelse(length(cause) > 1, i, cause) 
               # for competing risks
               if (fit$type == 'mright') {
                    
                    if (median) {
                         med.str <- as.character(med(fit, internal = T, cause = .cause))
                         med.tit <- 'median (95% CI)'
                    }
                    # if the model is without strata
                    if (ngroups == 1) {
                         
                         # for each cause, get the estimate and lower / upper bounds
                         tb2 <- round(cbind(tb$pstate[, .cause]*100, tb$lower[, .cause]*100, tb$upper[, .cause]*100), digits)
                         tab <- cbind(paste0(tb2[, 1], '% (', tb2[ , 2], '-', tb2[, 3], '%)'))
     
                         tab <- c('Cum Inc','', t(tab[-1]), pname, med.str)
                         finaltab <- c(finaltab, tab)
                         
                         # if the model is with covariates do the following:
                    } else {
                         varname.mat = matrix(unlist(strsplit(names(fit$strata), '=')), ncol = 2, byrow = TRUE)
                         
                         tb2 <- round(cbind(tb$pstate[, .cause]*100, tb$lower[, .cause]*100, tb$upper[, .cause]*100), digits)
                         tab <- cbind(paste0(tb2[, 1], '% (', tb2[ , 2], '-', tb2[, 3], '%)'))
                         tab <- matrix(tab, nrow = ngroups, byrow = T)
                         tab <- cbind(rep('', nrow(tab)), tab)
                         
                         if (test) {
                              gt <- rep(NA, ngroups)
                              gt[1] <- test.gray(fit, .cause)
                              pname <- 'p'
                         }
                         
                         # replace [,1] with varnames
                         tab[,1:2] <- varname.mat
                         tab[duplicated(tab[, 1]),1] <- NA
                         
                         tab <- cbind(tab, gt, med.str)
                         finaltab <- cbind(finaltab, tab)
                    }
               }
               # for single outcomes
               if (fit$type == 'right') {
                    if (median) {
                         med.str <- med(fit, internal = T)
                         med.tit <- 'median (95% CI)'
                    }
                    tb2 <- round(cbind(tb$surv*100, tb$lower*100, tb$upper*100), digits)
                    tab <- cbind(paste0(tb2[, 1], '% (', tb2[ , 2], '-', tb2[, 3], '%)'))
                    
                    # if the model is without strata
                    if (ngroups == 1) {
                         tab <- c('Survival', '', tab)
                         
                         # remove time = 0
                         tab <- c(tab[-2], as.character(med.str))
                         # convert to dataframe [note, transpose to force into rows
                         finaltab <- data.frame(t(tab))
                         column_names <- c(outcomes[1], paste0('Time: ', as.character(times)[-1]), med.tit)
                         
                         # if the model is with covariates do the following:
                    } else {
                         if (test) pval <- test.logrank(fit)
                         
                         varname.mat = matrix(unlist(strsplit(names(fit$strata), '=')), ncol = 2, byrow = TRUE)
                         
                         # reshape tab
                         tab <- matrix(tab, nrow = ngroups, byrow = T)
                         tab <- cbind(rep('', nrow(tab)), tab)
                         
                         # replace with row names
                         tab[ ,1:2] <- varname.mat
                         tab[duplicated(tab[, 1]),1] <-  NA
                         
                         if (test) {
                              tab <- cbind(tab, rep(NA, nrow(tab)))
                              tab[1, ncol(tab)] <- pval
                              pname <- 'p'
                         }
                         # insert into finaltab matrix
                         finaltab <- cbind(finaltab, tab)
                         finaltab <- cbind(finaltab, as.character(med.str))
                         finaltab <- data.frame(finaltab)
                    }
               }
          }
          if (ngroups == 1) {
               finaltab <- data.frame(t(finaltab))
          } else {
               finaltab <- data.frame(finaltab)
          }
          tb.final <- rbind(tb.final, finaltab)
     }
     # column names for the final table
     column_names <- rep(c('', ' ', paste0('Time: ', as.character(times)[-1]), pname, med.tit), length(cause))
     column_names[column_names == ''] <- outcomes
     
     options(knitr.kable.NA = '')
     return(knitr::kable(tb.final, col.names = column_names, row.names = FALSE, caption = caption))
}