### Input arguments
arg = commandArgs(T); ref.dat = arg[1]; info.file = arg[2]
tissue = arg[3]; qtl = arg[4]; chr = arg[5]; out.fname = arg[6]

### Packages
library(LAVA010)
library(parallel)

### Fixed args
univ.thresh=.0001
high.mem = c("brain_cortex", "brain_cerebellar", "brain_cerebellum", "cells_lymphocytes", "pituitary", "thyroid") # these are tissues that require a bit more memory and have to use fewer cores
high.mem.chrs = c("brain_caudate", "brain_frontal_cortex_b9", "brain_nucleus_accumbens", "adrenal_gland") 

### Set the number of parallel analyses
## for a subset of tissues and chrs the sqtl analyses fail, so they need to use fewer cores
n.parallel = 10
if (qtl == "sqtl") {
	if (tissue %in% high.mem | ((chr==19 | chr==16 | chr==8 | chr==9) & tissue %in% high.mem.chrs )) { n.parallel = 5 }
}
print(paste('Using',n.parallel,'cores!'))

### Read in data
gwas.input = process.input(info.file, sample.overlap.file=NULL, ref.dat, phenos="depression")
process.eqtl.input(gwas.input, qtl, tissue=tissue, is.sqtl=ifelse(qtl=="sqtl", T, F), chromosomes=chr)
genes = gwas.input$current.genes

### Analyse
print(paste("Starting analysis for",length(genes),"genes"))
progress = ceiling(quantile(1:length(genes), seq(.05,1,.05)))

u=b=NULL

out = mclapply(1:length(genes), mc.cores=n.parallel, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))

	locus = process.eqtl.locus(genes[i], gwas.input)
	
	if (!is.null(locus)) {
		loc.info = data.frame(locus = locus$id, chr = locus$chr, n.snps = locus$n.snps, n.pcs = locus$K)
		
		# run univ and bivar
		ub = run.univ.bivar(locus, univ.thresh=univ.thresh)
                u = cbind(loc.info, ub$univ)

                if (!is.null(ub$bivar)) { b = cbind(loc.info, ub$bivar) }
                return(list(univ=u, bivar=b))
	}
})

write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(out.fname,".univ"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind, lapply(out,"[[","bivar")), paste0(out.fname,".bivar"), row.names=F,quote=F,col.names=T)
print(paste0("Done! results written to ",out.fname,"*"))
