pack <- 'results/summary/chatgpt_pro_review_pack'
f5 <- read.csv('figures/fig5_left.csv')
rows <- lapply(f5$cell_id,function(id) {
 raw <- readRDS(sprintf('results/raw/cell_%04d.rds',id));ms <- raw$ms[raw$ms$est=='ML',]
 stopifnot(nrow(ms)==1000,!anyDuplicated(ms$rep),all(ms$cell_id==id))
 cbind(data.frame(cell_id=id,weak=raw$cell$weak),ms[setdiff(names(ms),'cell_id')])
})
out <- do.call(rbind,rows)
stopifnot(nrow(out)==24000,!anyDuplicated(out[c('cell_id','est','rep')]))
write.csv(out,file.path(pack,'evidence/figure5_left_replications.csv'),row.names=FALSE,na='NA')
cat('Figure5 left:',nrow(out),'MS rows;',ncol(out),'fields\n')
