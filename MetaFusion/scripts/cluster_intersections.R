library(RBGL)
library(data.table)

#Command args
args <- commandArgs(TRUE)
fid_intersection_file <- as.character(args[1])
fid_clusters_file <- as.character(args[2])

message("Reading FID.intersections.tsv file")
FIDs=fread(fid_intersection_file, header = TRUE, showProgress = TRUE)
FIDs$V1 = NULL

# For some reason pairToPair output has duplicates, need to remove them
FIDs.nodup <- FIDs[!duplicated(FIDs ), ]
message(nrow(FIDs.nodup), " records to build graph")


# Build graph
message("Building graph")
g <- ftM2graphNEL(as.matrix(FIDs.nodup), W=NULL, V=NULL, edgemode="directed")
edgemode(g) <- "undirected"

# "connections" is the object containing FID clusters
connections <- connectedComp(g)
connections <- lapply(connections, function(x){paste0(x, collapse=",")})
connections <- data.frame(FIDs=matrix(unlist(connections)),stringsAsFactors=FALSE)

#write tsv
fwrite(connections, file=fid_clusters_file, quote=FALSE, sep='\t', row.names=FALSE)
