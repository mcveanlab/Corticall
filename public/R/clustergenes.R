dm.file = "/Users/kiran/repositories/INDIANA/test.matrix";
groups.out = "/Users/kiran/repositories/INDIANA/cluster";

if (!interactive()) {
    args = commandArgs(TRUE);

    dm.file = args[1];
    groups.out = args[2];
}

suppressMessages(library(R.utils));

#message("Gathering table information...");
d.tmp = read.table(dm.file, header=TRUE, stringsAsFactors=FALSE, nrows=1);
cs = c("character", rep("numeric", ncol(d.tmp) - 1));
rows = countLines(dm.file);

#message("Reading table...");
d = read.table(dm.file, header=TRUE, stringsAsFactors=FALSE, colClasses=cs, nrows=rows);

#message("Clustering...");
d.dist = as.dist(as.matrix(d));
h = hclust(d.dist);

#message("Extracting groups...");
seqs = h$labels[h$order];

groups = list();
groupIndex = 1;
groups[[groupIndex]] = c(seqs[1]);

cols = gsub("[:.]", "_", colnames(d));

for (i in 1:(length(seqs)-1)) {
    seq1 = gsub("[:.]", "_", seqs[i]);
    seq2 = gsub("[:.]", "_", seqs[i+1]);

    id1 = which(cols == seq1);
    id2 = which(cols == seq2);

    val = d[id1, id2];

    if (val == 1) { # new group
        groupIndex = groupIndex + 1;
        groups[[groupIndex]] = c(seq2);
    } else {
        groups[[groupIndex]] = c(groups[[groupIndex]], seq2);
    }
}

numGroups.unfiltered = length(groups);

# Filter gene lists
elementsToKeep = c();
for (i in 1:length(groups)) {
    if (length(groups[[i]]) > 3) {
        elementsToKeep = c(elementsToKeep, i);
    }
}

newGroups = groups[elementsToKeep];
numGroups.filtered = length(newGroups);

#message("Writing groups...");
write(paste("groupName", "groupMembers", sep="\t"), groups.out, append=FALSE, ncolumns=1000);
for (i in 1:length(newGroups)) {
    #groupName = names(sort(table(desc[newGroups[[i]],]$description), decreasing=TRUE))[1];
    groupName = i;
    groupMembers = paste(newGroups[[i]], collapse=",");

    write(paste(groupName, groupMembers, sep="\t"), groups.out, append=TRUE, ncolumns=1000);
}

#cat("Unfiltered groups: ", numGroups.unfiltered, "\n");
#cat("Filtered groups: ", numGroups.filtered, "\n");
