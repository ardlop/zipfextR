pkgname <- "zipfextR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('zipfextR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("moezipf")
### * moezipf

flush(stderr()); flush(stdout())

### Name: moezipf
### Title: The Marshal-Olkin Extended Zipf Distribution (MOEZipf).
### Aliases: moezipf dmoezipf dmoezipf pmoezipf qmoezipf rmoezipf

### ** Examples

dmoezipf(1:10, 2.5, 1.3)
pmoezipf(1:10, 2.5, 1.3)
qmoezipf(0.56, 2.5, 1.3)
rmoezipf(10, 2.5, 1.3)




cleanEx()
nameEx("moezipfFit")
### * moezipfFit

flush(stderr()); flush(stdout())

### Name: moezipfFit
### Title: MOEZipf parameters estimation.
### Aliases: moezipfFit residuals.moezipfR fitted.moezipfR coef.moezipfR
###   plot.moezipfR print.moezipfR summary.moezipfR logLik.moezipfR
###   AIC.moezipfR BIC.moezipfR

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
data <- as.data.frame(table(data))
data[,1] <- as.numeric(data[,1])
initValues <- moezipf_getInitialValues(data)
obj <- moezipfFit(data, init_alpha = initValues$init_alpha, init_beta = initValues$init_beta)



cleanEx()
nameEx("moezipfMean")
### * moezipfMean

flush(stderr()); flush(stdout())

### Name: moezipfMean
### Title: Expected value.
### Aliases: moezipfMean

### ** Examples

moezipfMean(2.5, 1.3)
moezipfMean(2.5, 1.3, 10^(-3))



cleanEx()
nameEx("moezipfMoments")
### * moezipfMoments

flush(stderr()); flush(stdout())

### Name: moezipfMoments
### Title: Distribution Moments.
### Aliases: moezipfMoments

### ** Examples

moezipfMoments(3, 4.5, 1.3)
moezipfMoments(3, 4.5, 1.3,  1*10^(-3))



cleanEx()
nameEx("moezipfVariance")
### * moezipfVariance

flush(stderr()); flush(stdout())

### Name: moezipfVariance
### Title: Variance of the MOEZipf distribution.
### Aliases: moezipfVariance

### ** Examples

moezipfVariance(3.5, 1.3)



cleanEx()
nameEx("moezipf_getInitialValues")
### * moezipf_getInitialValues

flush(stderr()); flush(stdout())

### Name: moezipf_getInitialValues
### Title: Calculates initial values for the alpha and beta parameters.
### Aliases: moezipf_getInitialValues

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
data <- as.data.frame(table(data))
initials <- moezipf_getInitialValues(data)



cleanEx()
nameEx("zipfpe")
### * zipfpe

flush(stderr()); flush(stdout())

### Name: zipfpe
### Title: The Zipf-Poisson Extreme Distribution (Zipf-PE).
### Aliases: zipfpe dzipfpe dzipfpe pzipfpe qzipfpe rzipfpe

### ** Examples

dzipfpe(1:10, 2.5, -1.5)
pzipfpe(1:10, 2.5, -1.5)
qzipfpe(0.56, 2.5, 1.3)
rzipfpe(10, 2.5, 1.3)




cleanEx()
nameEx("zipfpeFit")
### * zipfpeFit

flush(stderr()); flush(stdout())

### Name: zipfpeFit
### Title: Zipf-PE parameters estimation.
### Aliases: zipfpeFit residuals.zipfpeR fitted.zipfpeR coef.zipfpeR
###   plot.zipfpeR print.zipfpeR summary.zipfpeR logLik.zipfpeR AIC.zipfpeR
###   BIC.zipfpeR

### ** Examples

data <- rzipfpe(100, 2.5, 1.3)
data <- as.data.frame(table(data))
data[,1] <- as.numeric(data[,1])
obj <- zipfpeFit(data, 1.1, 0.1)



cleanEx()
nameEx("zipfpeMean")
### * zipfpeMean

flush(stderr()); flush(stdout())

### Name: zipfpeMean
### Title: Expected value of the Zipf-PE distribution.
### Aliases: zipfpeMean

### ** Examples

zipfpeMean(2.5, 1.3)
zipfpeMean(2.5, 1.3, 10^(-3))



cleanEx()
nameEx("zipfpeMoments")
### * zipfpeMoments

flush(stderr()); flush(stdout())

### Name: zipfpeMoments
### Title: Distribution Moments.
### Aliases: zipfpeMoments

### ** Examples

zipfpeMoments(3, 4.5, 1.3)
zipfpeMoments(3, 4.5, 1.3,  1*10^(-3))



cleanEx()
nameEx("zipfpeVariance")
### * zipfpeVariance

flush(stderr()); flush(stdout())

### Name: zipfpeVariance
### Title: Variance of the Zipf-PE distribution.
### Aliases: zipfpeVariance

### ** Examples

zipfpeVariance(3.5, 1.3)



cleanEx()
nameEx("zipfpssFit")
### * zipfpssFit

flush(stderr()); flush(stdout())

### Name: zipfpssFit
### Title: Zipf-PSS parameters estimation.
### Aliases: zipfpssFit residuals.zipfpssR fitted.zipfpssR coef.zipfpssR
###   plot.zipfpssR print.zipfpssR summary.zipfpssR logLik.zipfpssR
###   AIC.zipfpssR BIC.zipfpssR

### ** Examples

data <- rzipfpss(100, 2.5, 1.3)
data <- as.data.frame(table(data))
data[,1] <- as.numeric(data[,1])
obj <- zipfpssFit(data, 1.1, 0.1)



cleanEx()
nameEx("zipfpssMean")
### * zipfpssMean

flush(stderr()); flush(stdout())

### Name: zipfpssMean
### Title: Expected value of the Zipf-PSS distribution.
### Aliases: zipfpssMean

### ** Examples

zipfpssMean(2.5, 1.3)
zipfpssMean(2.5, 1.3, TRUE)



cleanEx()
nameEx("zipfpssMoments")
### * zipfpssMoments

flush(stderr()); flush(stdout())

### Name: zipfpssMoments
### Title: Distribution Moments.
### Aliases: zipfpssMoments

### ** Examples

zipfpssMoments(1, 2.5, 2.3)
zipfpssMoments(1, 2.5, 2.3, TRUE)



cleanEx()
nameEx("zipfpssVariance")
### * zipfpssVariance

flush(stderr()); flush(stdout())

### Name: zipfpssVariance
### Title: Variance of the Zipf-PSS distribution.
### Aliases: zipfpssVariance

### ** Examples

zipfpssVariance(4.5, 2.3)
zipfpssVariance(4.5, 2.3, TRUE)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
