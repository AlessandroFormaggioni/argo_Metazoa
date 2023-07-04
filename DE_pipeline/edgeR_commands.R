library(edgeR)
t1=read.table("ALL_readcount.tab",row.names=1,header=F)
group=c(1,1,1,2,2,2)
group=factor(group)
y=DGEList(counts=t1,group=group)
design <- model.matrix(~0+group, data=y$samples)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y=estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf<- glmQLFTest(fit, contrast=c(-1,1))
sign_table=qlf$table[qlf$table$PValue < 0.001, ]



