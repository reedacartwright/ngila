library(gsl)

x <- seq(1,10.1,0.1)
y <- log((x-1)*zeta(x))
y[1] <- 0
f <- splinefun(x,y)
z <- get("z", envir = environment(f))

fmt <- "%.17e"

m <- matrix(sprintf(fmt,z$y),ncol=4,byrow=T)
cat("\tstatic double y[92] = {\n\t\t")
cat(paste(apply(m,1,function(x) paste(x,collapse=", ")), collapse=",\n\t\t"))
cat("\n\t};\n")

m <- matrix(sprintf(fmt,z$b),ncol=4,byrow=T)
cat("\tstatic double b[92] = {\n\t\t")
cat(paste(apply(m,1,function(x) paste(x,collapse=", ")), collapse=",\n\t\t"))
cat("\n\t};\n")

m <- matrix(sprintf(fmt,z$c),ncol=4,byrow=T)
cat("\tstatic double c[92] = {\n\t\t")
cat(paste(apply(m,1,function(x) paste(x,collapse=", ")), collapse=",\n\t\t"))
cat("\n\t};\n")

m <- matrix(sprintf(fmt,z$d),ncol=4,byrow=T)
cat("\tstatic double d[92] = {\n\t\t")
cat(paste(apply(m,1,function(x) paste(x,collapse=", ")), collapse=",\n\t\t"))
cat("\n\t};\n")

