root <- '.'
pack <- file.path(root,'results/summary/chatgpt_pro_review_pack')
ids <- scan(file.path(root,'evidence/selected_cell_ids.txt'),quiet=TRUE)
f <- readRDS(file.path(root,'results/fits1f.rds'))
selected <- f[f$cell_id %in% ids,]
stopifnot(!anyDuplicated(selected[c('cell_id','est','rep')]),
          setequal(unique(selected$cell_id),ids),sum(selected$min_psi<=0,na.rm=TRUE)==7)
con <- gzfile(file.path(pack,'evidence/review_fits1f.csv.gz'),'wt')
write.csv(selected,con,row.names=FALSE,na='NA');close(con)
pop <- readRDS(file.path(root,'results/population.rds'))
write.csv(pop,file.path(pack,'evidence/population.csv'),row.names=FALSE,na='NA')
tab <- read.csv(file.path(root,'figures/table4.csv'))
pairs <- lapply(seq_len(nrow(tab)),function(k) do.call(rbind,lapply(c(500,2000),function(N) {
 id <- tab[[paste0('cell_id_',N)]][k]
 raw <- readRDS(file.path(root,sprintf('results/raw/cell_%04d.rds',id)))
 ms <- raw$ms[raw$ms$est=='ML',]
 d <- selected[selected$cell_id==id & selected$est=='ML',]
 stopifnot(nrow(ms)==1000,!anyDuplicated(ms$rep),setequal(ms$rep,d$rep))
 d <- d[match(ms$rep,d$rep),]
 cols <- c('rep','best_fx','n_distinct','n_admissible','any_heywood_min','def_conv',
           'def_nonglobal','def_on_weak','def_rev_nonglobal','def_rev_on_weak',
           'global_on_weak','spurious_weak','frac_starts_weak')
 out <- cbind(data.frame(table_row=k,cell_id=id,cond=tab$cond[k],phi=tab$phi[k],N=N),
              d[c('fx','chisq','df','cfi','tli','rmsea','converged','min_psi')],ms[cols])
 names(out)[names(out)=='fx'] <- 'default_fx'
 out
})))
pairs <- do.call(rbind,pairs)
stopifnot(nrow(pairs)==18000,!anyDuplicated(pairs[c('cell_id','rep')]))
write.csv(pairs,file.path(pack,'evidence/table4_replications.csv'),row.names=FALSE,na='NA')
cat('Selected fits:',nrow(selected),'rows;',ncol(selected),'fields;',length(ids),'cells\n')
cat('All7 negative uniqueness fits included. Table4:',nrow(pairs),'paired rows. Population:',nrow(pop),'rows\n')
