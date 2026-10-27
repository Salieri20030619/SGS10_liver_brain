root <- "/mnt/BioDisk/wangmingyuan/SGS10_analysis/TableS3_clean_20260907_v1"
warnings_seen <- character(); identical_values <- logical()
for(organ in c("brain","liver")) {
 d <- read.csv(file.path(root,paste0(organ,"_top30_per_sample.csv")),check.names=FALSE)
 groups <- split(seq_len(nrow(d)),paste(d$cluster,d$sample,sep="|"))
 for(ii in groups) {
  z <- d[ii,];s<-z$sample[1]
  rebuilt <- withCallingHandlers(data.frame(organ=z$organ[1],cluster=z$cluster[1],
    gene=z$gene,sample=s,n_cells_target=z$n_cells_target[1],
    n_cells_rest=setNames(z$n_cells_rest[1],s),
    evaluable_for_consistency=setNames(z$evaluable_for_consistency[1],s),
    mean_normalized_target=setNames(z$mean_normalized_target,z$gene),
    mean_normalized_rest=setNames(z$mean_normalized_rest,z$gene),
    pct_target=setNames(z$pct_target,z$gene),pct_rest=setNames(z$pct_rest,z$gene),
    log2FC=setNames(z$log2FC,z$gene),delta_pct=setNames(z$delta_pct,z$gene)),
    warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
  identical_values<-c(identical_values,isTRUE(all.equal(z,rebuilt,check.attributes=FALSE)))
 }
}
result <- c(paste("table_groups_checked:",length(identical_values)),
            paste("all_values_identical:",all(identical_values)),
            paste("reproduced_warning_count:",length(warnings_seen)),
            "reproduced_warning_messages:",unique(warnings_seen),
            "This reproduces the named-scalar row-name warning in table assembly; statistical values are unchanged.")
writeLines(result,file.path(root,"table_assembly_warning_audit.txt"))
cat(paste(result,collapse="\n"),"\n")
stopifnot(all(identical_values))
