# synteny_2d_plot.R
# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# created by WRF 2019-04-01
# last modified 2020-03-01

#args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py
#all2Dfile = args[1]
all2Dfile = "~/genomes/hormiphora_californensis/hormiphora/synteny/hcalv1_v_ml2_2d_synteny_points_local_ml109.tab"

# read optional species names
#genome1_arg = args[2]
#genome2_arg = args[3]

genome1_arg = "Mnemiopsis leidyi"
genome2_arg = "Hormiphora californensis"

if (!is.na(genome1_arg)) {
xlab = paste( gsub("-", " ", genome1_arg),"(total Mb)")
} else {
xlab = "Genome 1 (total Mb)"
}

if (!is.na(genome2_arg)) {
ylab = paste( gsub("-", " ", genome2_arg),"(total Mb)")
} else {
ylab = "Genome 2 (total Mb)"
}

# read all data in a single file
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
#head(all2Ddata)

# breakdown categories of either two genomes s1 and s2, or gene hits g
categories = all2Ddata[,1]

is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
scafdata1

longestscaf1 = max(scafdata1[,6])
#longestscaf1
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0009)
#is_longscafs1

longscafs1 = c(0, scafdata1[,6][is_longscafs1] )

is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.00001)
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )

is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
#head(pointsdata)


genome_x = pointsdata[,7]
genome_y = pointsdata[,6]
xmax = tail( pretty(longscafs2), n=1)
ymax = tail( pretty(longscafs1), n=1)
longscafs_x = longscafs2
longscafs_y = longscafs1
nscafs_x = length(longscafs2)
nscafs_y = length(longscafs1)

xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)

# larger bitscores make larger points
pointsize = log10(as.numeric(pointsdata[,8])) / 4


matchcounts_by_scaf = table(pointsdata[,3], pointsdata[,5])
#matchcounts_by_scaf
match_by_ml = t(matchcounts_by_scaf)
max_by_ml = apply(match_by_ml,1,max)
sum_by_ml = apply(match_by_ml,1,sum)
#hist(max_by_ml/sum_by_ml, breaks=20, xlim=c(0,1), col="#2345a9")
pdf(file="~/genomes/hormiphora_californensis/hormiphora/synteny/hcalv1_v_ml2_2d_synteny_points.most_to_scaf.pdf", width=6, height=6)
par(mar=c(4.5,4.5,3,1))
plot(max_by_ml/sum_by_ml,sum_by_ml, type='p', pch=21, bg="#2345c9aa",xlab="Most to single Hcal scaffold/total", ylab="Num genes on Mlei scaffold", main="Hcal vs Mlei most hits to a single scaffold", cex.lab=1.4, cex.axis=1.4, cex.main=1.5, cex=1.5)
dev.off()

top_match_scaf = row.names(matchcounts_by_scaf)[apply(matchcounts_by_scaf, 2, which.max)]
#top_match_scaf
top_match_scaf_index = match(top_match_scaf, scafdata1[,2])
#top_match_scaf_index
sorted_top_match_index = sort(top_match_scaf_index, index.return=TRUE)
#sorted_top_match_index
scafname_sorted_by_match = colnames(matchcounts_by_scaf)[sorted_top_match_index$ix]
#scafname_sorted_by_match
scafsize_sorted_by_match = scafdata2[,4][match(scafname_sorted_by_match, scafdata2[,2])]
#scafsize_sorted_by_match
sorted_offsets1 = c(0, cumsum(scafdata1[,4]) )
sorted_offsets2 = c(0, cumsum(scafsize_sorted_by_match)[1:(length(scafsize_sorted_by_match)-1)])

genome_x = pointsdata[,7] + sorted_offsets2[match(pointsdata[,5],scafname_sorted_by_match)]
genome_y = pointsdata[,6] + sorted_offsets1[match(pointsdata[,3],scafdata1[,2])]



# dotcolor = "#18935188" # green
# dotcolor = "#c51b8a88" # pink
# dotcolor = "#1c909988" # teal
# dotcolor = "#2071d388" # blue
# dotcolor = "#88419d88" # purple


# make PDF
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",all2Dfile,perl=TRUE)
pdf(file=outputfile, width=8, height=11)
par( mar=c(4.5,4.5,1,1) )
plot(genome_x, genome_y, pch=16, xlim=c(0,109000000), ylim=c(0,110000000), col="#88419d88", cex=pointsize, main=basename(all2Dfile), xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.4)

tickpoints = pretty(c(0,100))
axis(1, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3)
barpos_x = rep(c( ymax*-0.01, ymax*-0.02, ymax*-0.03),round(nscafs_x)/3)
segments( sorted_offsets2[1:(nscafs_x-1)], barpos_x[1:(nscafs_x-1)], sorted_offsets2[2:nscafs_x], barpos_x[0:(nscafs_x-1)], lwd=2)
#segments(sorted_offsets2, 0, sorted_offsets2, longscafs_y[nscafs_y], lwd=0.1, col="#777777")

tickpoints = pretty(c(0,100))
axis(2, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3, line=0.5)
barpos_y = rep(c( xmax*-0.012, xmax*-0.027),round(nscafs_y)/2)
segments( barpos_y[1:(nscafs_y-1)], longscafs_y[1:(nscafs_y-1)], barpos_y[0:(nscafs_y-1)], longscafs_y[2:nscafs_y], lwd=3)
#segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.1, col="#777777")

# display numbers beside y-axis scaffold segments
textpos_y = rep(c( xmax*-0.027, xmax*-0.012 ),round(13)/2)
textmidbar = as.numeric(scafdata1[1:13,6]) - as.numeric(scafdata1[1:13,4])/2
text(textpos_y, textmidbar, 1:13, cex=0.5)

dev.off()




#