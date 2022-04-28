


library(factoextra)


args <- commandArgs(trailingOnly = TRUE)
S = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))  # can be symetric matrix
nbClusters = as.numeric(args[2])
unitig_names = rownames(S)
print(dim(S))


S = dist(S)
hc <- hclust(S, method = "ward.D2")
clusters = cutree(hc, k = nbClusters)

df = data.frame(Name=unitig_names, Colour=clusters)
print(df)

write.csv(df, paste0(args[1], ".csv"), row.names = FALSE, quote=FALSE)


