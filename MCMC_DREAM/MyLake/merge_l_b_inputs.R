## R --vanilla --args LAE_input_basin2_land.txt LAE_input_basin1_to_basin2.txt LAE_input_basin2.txt < merge\ land\ and\ basin\ inputs.R 

args <- commandArgs(trailingOnly=TRUE) ## after --args
fromlandfn <- args[1]
fromupperbasinfn <- args[2]
outfn <- args[3]

f <- file(fromlandfn, 'rt')
firsttwolines <- readLines(f, 2)  # two strings without EOL
data1 <- read.table(f, sep='\t', na.strings=c('NaN', ''))[,1:20]
close(f)

data2 <- read.table(fromupperbasinfn, sep=',')

if (sum(is.nan(unlist(c(data1))) > 0)) {
  stop('error in land inpit see script')
}
#if (sum(is.nan(unlist(c(data2))) > 0)) {
#  stop('error in upper basin see script')
#}


v1 <- data1[,11]
v2 <- data2[,1] # from upper basin
print(dim(data1))
print(dim(data2))
readline("Inspect lenght, press Enter") 
## print(dim(diag(v1)))

temp <- array(NA, dim = c(nrow(data1), 9))
temp[, 1] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 1)] + v2 * data2[, (1 + 1)])
temp[, 2] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 2)] + v2 * data2[, (1 + 2)])
temp[, 3] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 3)] + v2 * data2[, (1 + 3)])
temp[, 4] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 4)] + v2 * data2[, (1 + 4)])
temp[, 5] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 5)] + v2 * data2[, (1 + 5)])
temp[, 6] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 6)] + v2 * data2[, (1 + 6)])
temp[, 7] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 7)] + v2 * data2[, (1 + 7)])
temp[, 8] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 8)] + v2 * data2[, (1 + 8)])
temp[, 9] <- (1 / (v1 + v2)) * (v1 * data1[, (11 + 9)] + v2 * data2[, (1 + 9)])
## temp <- diag(1 /(v1 + v2)) %*% (diag(v1) %*% as.matrix(data1[, 12:18]) +
##                                 diag(v2) %*% as.matrix(data2[, 2:8]))
# volume weighted concentration
out <- data.frame(data1[, 1:10], v1 + v2, temp)
fn <- outfn
f <- file(fn, 'wt')
cat(firsttwolines, sep='\n', file=f)
cat('\n', file=f)
write.table(out, file=f, # appends anyway to a file connection
            sep='\t', na='NaN', row.names=FALSE, col.names=FALSE)
close(f)

