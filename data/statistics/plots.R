setwd("/Users/enricoseiler/git/minimizer_ibf/data/statistics/")

titles <- c("23/20", "24/20", "25/20", "40/20", "100/20", "1000/20", "23/20 xor", "24/20 xor", "25/20 xor", "40/20 xor", "100/20 xor", "1000/20 xor")

for (file in c("plot00", "plot01", "plot02", "plot03", "plot04", "plot05", "plot06", "plot07", "plot08", "plot09", "plot10", "plot11")) {

  table <- read.csv(paste(getwd(), '/', file, sep=''), header=FALSE, sep=',') # 23, 20
  table <- table[2:nrow(table), ]
  table <- table[min(which(table$V2 !=0)):max(which(table$V2 !=0)),]
  table <- transform(table, V2=V2/sum(V2))

  png(paste(getwd(), '/', file, '.png', sep=''), width=12*300, height=8*300, res=300, pointsize = 8)
  barplot(table$V2, names.arg = table$V1, main=titles[1])
  dev.off()
  titles <- titles[-1]

}
