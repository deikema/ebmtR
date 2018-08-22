# LIST VARIABLE N FUNCTION ARGUMENTS
unpack.list <- function(object) for (.x in names(object)) assign(value = object[[.x]], x = .x, envir = parent.frame())

# FIND LAST INDEX
last <- function(x) { return( x[length(x)] ) }

# COMPARISON IN VECTORS WITH NA'S
'%==%' <- function(a, b) a == b & !is.na(a) & !is.na(b)

# LOGRANK
test.logrank <- function(fit) {
     .logrank <- survdiff(formula(fit$call$formula), get(as.character(fit$call$data)))
     .pval <- round(1 - pchisq(.logrank$chisq, length(fit$strata)), 4) 
     .lt <- ifelse(.pval < 0.001, 'p<0.001', paste0('p=', .pval))
     return(.lt)
}

# GRAY
test.gray <- function(fit, .cause) {
     .dat <- as.character(fit$call$data)
     .out <- get(.dat)[, c(unlist(strsplit(as.character(fit$call$formula[2]),'[,(=)]|\\s'))[c(2,4)], 
                          as.character(fit$call$formula[3]))]
     .gray.test <- round(cuminc(.out[, 1], .out[, 2], .out[, 3], cencode = 0)$Tests[.cause, 2], 4)
     .gt <- ifelse(.gray.test < 0.001, 'p<0.001', paste0('p=', .gray.test))   
     return(.gt)
}

# HELPFUL THINGS
# SPSS DATES: SPSS USES SECONDS SINCE START OF GREG CALENDAR (although technically it was implemented on feb 24, this is not what spss does)
# as.POSIXlt(x, origin="1582/10/14")
