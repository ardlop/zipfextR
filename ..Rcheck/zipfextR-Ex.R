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
### Aliases: moezipfFit

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
data <- zipfExtR_getDataMatrix(data)
obj <- moezipfFit(data, 1.001, 0.001)



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
### Title: Variance.
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
data <- zipfExtR_getDataMatrix(data)
initials <- moezipf_getInitialValues(data)



cleanEx()
nameEx("zipfExtR_getDataMatrix")
### * zipfExtR_getDataMatrix

flush(stderr()); flush(stdout())

### Name: zipfExtR_getDataMatrix
### Title: Convert a sample vector to a frequency matrix.
### Aliases: zipfExtR_getDataMatrix

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
zipfExtR_getDataMatrix(data)




cleanEx()
nameEx("zpe")
### * zpe

flush(stderr()); flush(stdout())

### Name: zpe
### Title: The Zipf-Poisson Extreme Distribution (ZPE).
### Aliases: zpe dzpe dzpe pzpe qzpe rzpe

### ** Examples

dzpe(1:10, 2.5, -1.5)
pzpe(1:10, 2.5, -1.5)
qzpe(0.56, 2.5, 1.3)
rzpe(10, 2.5, 1.3)




cleanEx()
nameEx("zpeFit")
### * zpeFit

flush(stderr()); flush(stdout())

### Name: zpeFit
### Title: ZPE parameters estimation.
### Aliases: zpeFit

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
data <- zipfExtR_getDataMatrix(data)
obj <- zpeFit(data, 1.001, 0.001)



cleanEx()
nameEx("zpeMoments")
### * zpeMoments

flush(stderr()); flush(stdout())

### Name: zpeMoments
### Title: Distribution Moments.
### Aliases: zpeMoments

### ** Examples

moezipfMoments(3, 4.5, 1.3)
moezipfMoments(3, 4.5, 1.3,  1*10^(-3))



cleanEx()
nameEx("zpssFit")
### * zpssFit

flush(stderr()); flush(stdout())

### Name: zpssFit
### Title: ZPSS parameters estimation.
### Aliases: zpssFit

### ** Examples

data <- rmoezipf(100, 2.5, 1.3)
data <- zipfExtR_getDataMatrix(data)
obj <- zpssFit(data, 1.001, 0.001)



cleanEx()
nameEx("zpssMean")
### * zpssMean

flush(stderr()); flush(stdout())

### Name: zpssMean
### Title: Expected value of the Z-PSS distribution.
### Aliases: zpssMean

### ** Examples

zpssMean(2.5, 1.3)
zpssMean(2.5, 1.3, TRUE)



cleanEx()
nameEx("zpssVariance")
### * zpssVariance

flush(stderr()); flush(stdout())

### Name: zpssVariance
### Title: Variance of the Z-PSS distribution.
### Aliases: zpssVariance

### ** Examples

zpssVariance(4.5, 2.3)
zpssVariance(4.5, 2.3, TRUE)



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
