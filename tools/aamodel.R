# order of nucleotides in the file
aa <- unlist(strsplit("ARNDCQEGHILKMFPSTWYV",""))

# load table
v <- scan("lg_LG.PAML.txt")
u <- v[1:190]
f <- v[191:210]
s <- matrix(0,20,20)
s[upper.tri(s)] <- u
s <- t(s)
s[upper.tri(s)] <- u

# create subst matrix
m <- t(s*f)
diag(m) <- -apply(m,1,sum)
m <- t(m*sqrt(f))/sqrt(f)

# decompose
em <- eigen(m,symmetric=TRUE)
#em$values[abs(em$values) < 2*.Machine$double.eps] <- 0

v <- c(sqrt(f),em$values,t(em$vectors))
v <- matrix(sprintf("% .17e", v),ncol=4,byrow=T)
cat(paste(apply(v,1,function(x) paste(x,collapse=", ")), collapse=",\n"))
cat("\n")
