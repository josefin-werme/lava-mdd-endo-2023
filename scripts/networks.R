arg = commandArgs(T); ref.dat = arg[1]; loc.file = arg[2]; info.file = arg[3]
phenos = unlist(strsplit(arg[4],";")); sample.overlap.file = arg[5]; out.fname = arg[6]

library(LAVA); library(parallel)

### Fixed args
univ.thresh=.0001

### Read in data
loci = read.loci(loc.file)
n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.dat, phenos)	# process sumstats etc

### Analyse
print(paste("Starting analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))

u=b=NULL

out = mclapply(1:n.loc, mc.cores=15, function(i) {
        if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))  # print progress

        # process locus
        locus = process.locus(loci[i,], input)

        # in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
        if (!is.null(locus)) {
                # extract general locus info for output
                loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

                # run univ & bivar analysis functions & store results
                ub = run.univ.bivar(locus, univ.thresh=univ.thresh)
                u = cbind(loc.info, ub$univ)
                if (!is.null(ub$bivar)) {
                        b = cbind(loc.info, ub$bivar)
                }
		return(list(univ=u, bivar=b))
        }
})

### WRITE OUTPUT ###
write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(out.fname,".univ"), row.names=F, quote=F)
write.table(do.call(rbind, lapply(out, "[[", "bivar")), paste0(out.fname,".bivar"), row.names=F, quote=F)
print(paste0("Done! results written to ",out.fname,"*"))
