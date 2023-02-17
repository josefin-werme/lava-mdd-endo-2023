### Input arguments
arg = commandArgs(T); ref.dat = arg[1]; loc.file = arg[2]; info.file = arg[3]
phenos = unlist(strsplit(arg[4],";")); sample.overlap.file = arg[5]; out.fname = arg[6]

# ref.dat="g1000_eur.maf005"; loc.file="gencode.v26.GRCh38.protein_coding.1Mb-cis.SNP-IDs.loci"; info.file="input.info.brain.txt"; phenos=c("accumbens_area_left","accumbens_area_right","depression"); sample.overlap.file="sample.overlap.dat"

### Packages
library(LAVA010)
library(parallel)

### Fixed args
univ.thresh=.0001
meta.r2.thresh = .9
brain.phenos = phenos[phenos!="depression"]
n.regions = length(brain.phenos)

### Read in data
loci = read.loci(loc.file)
n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.dat, phenos) 

### Analyse
print(paste("Starting analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))

out = mclapply(1:n.loc, mc.cores=10, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))

	locus = process.locus(loci[i,], input)
	
	# check that depression and at least one more region is present
	if (!is.null(locus) & "depression" %in% locus$phenos & length(locus$phenos) > 1) {
		loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
		
		# run univ	
		univ = run.univ(locus)
		u = cbind(loc.info, univ)

		# check which phenos have univ signal
		univ.sig = as.character(subset(univ, p < univ.thresh)$phen)
		
		# check that at least dep + one brain phenotype has signal
		b = c = NULL
		if ("depression" %in% univ.sig & length(univ.sig) > 1) {
			if (length(univ.sig) > 2) {
				# if region is lateralised and both pass univ thresh, evaluate how strongly correlated homologous regions are
				brain.bivar = run.bivar(locus, phenos=brain.phenos, p.values=F)
				c = cbind(loc.info, brain.bivar$rho, brain.bivar$r2, brain.bivar$r2.upper)
				
				if (brain.bivar$r2.upper == 1 | brain.bivar$r2 > meta.r2.thresh) {
					# if upper r2 CI == 1 or r2 greater than the specified threshold (.9), meta-analyse the homologous regions 
					meta.locus = meta.analyse.locus(locus, brain.phenos)
					b = cbind(loc.info, run.bivar(meta.locus))
				} else {
					# if hemispheres are not correlated enough, analyse separately
					b = cbind(loc.info, run.bivar(locus, target="depression")) 
				}
			} else {
				# if region is not lateralised or only one hemisphere reaches univ thresh
				b = cbind(loc.info, run.bivar(locus, phenos = univ.sig)) # this will only ever be depression + one phenotype so no need to set target
			}
		}
		return(list(univ=u, bivar=b, cor=c))
	}
})

write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(out.fname,".univ"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind, lapply(out,"[[","bivar")), paste0(out.fname,".bivar"), row.names=F,quote=F,col.names=T)
if (n.regions==2) write.table(do.call(rbind, lapply(out,"[[","cor")), paste0(out.fname,".cor"), row.names=F,quote=F,col.names=T)
print(paste0("Done! results written to ",out.fname,"*"))
