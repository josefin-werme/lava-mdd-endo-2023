library(data.table)
setwd("~/Documents/GitHub/mdd-endo-scripts-and-data2023") ### NOTE: THIS NEEDS TO BE ADJUSTED

######################################
######## READ IN ALL THE DATA ########
######################################
## Gene file + GTEx tissues
tissues = read.table("scripts/tissues.txt", header=F, stringsAsFactors = F)$V1
g = fread("data/analysis_files/gencode.v26.GRCh38.protein_coding.gtf",header=F, sep=" ") # protein coding genes

#### read in eQTLs
e = list()
for (i in tissues) {
	e[[i]] = fread(paste0("data/results/local_rg/eqtl/",i,".chrall.bivar"), data.table=F, header=T)
	e[[i]]$gene = g[match(e[[i]]$locus, g$V7),]$V9
	e[[i]]$tissue = i
	e[[i]]$start = g[match(e[[i]]$locus, g$V7),]$V4
	e[[i]]$stop = g[match(e[[i]]$locus, g$V7),]$V5
}
e = do.call(rbind,e)
e = subset(e,!is.na(gene)) # remove non-coding genes (these were not meant to be analysed, but were present in the lava processed gtex files)


#### read in sQTLs
s = list()
for (i in tissues) {
	s[[i]] = NULL
	for (j in 1:22) {
		# check that all files exist since try is used below (with silent=T since those errors are not informative and clogs the terminal)
		infile = paste0("data/results/local_rg/sqtl/",i,".chr",j,".bivar")
		if (!file.exists(infile)) { print(paste(infile,"does not exist")) }
		
		# using try here since some output files have zero rows due to lack of univariate signal
		tmp = try(fread(infile, data.table=F, header=T), silent=T)
		if (class(tmp) != "try-error") { s[[i]] = rbind(s[[i]], tmp) }
	}
	s[[i]]$ens = sapply(strsplit(s[[i]]$locus,'__'),'[[',1)
	s[[i]]$gene = g[match(s[[i]]$ens, g$V7),]$V9
	s[[i]]$tissue = i
	s[[i]]$start = g[match(s[[i]]$ens, g$V7),]$V4
	s[[i]]$stop = g[match(s[[i]]$ens, g$V7),]$V5
}
s = do.call(rbind,s)
s = subset(s,!is.na(gene)) # remove non-coding genes (these were never meant to be analysed, but were present in the lava processed gtex files)


#### read in MRI
regions = read.table("scripts/regions.txt",header=F, stringsAsFactors=F)$V1

mri.b = list()
for (i in regions) {
	infile = paste0("data/results/local_rg/regions/",i,".bivar")
	mri.b[[i]] = fread(infile, data.table=F)
	mri.b[[i]]$gene = g[match(mri.b[[i]]$locus, g$V7),]$V9
}
m.dat = do.call(rbind, mri.b)

#### read in networks
networks = read.table("scripts/networks.txt",header=F, stringsAsFactors=F)$V1

net.b  = list()
for (i in networks) {
	infile = paste0("data/results/local_rg/networks/",i,".bivar")
	net.b[[i]] = fread(infile, data.table=F)
	net.b[[i]]$gene = g[match(net.b[[i]]$locus, g$V7),]$V9
}
net = do.call(rbind, net.b)





#####################################
######## ENRICHMENT ANALYSIS ########
#####################################

### Read in gene sets
sets = read.table("data/analysis_files/magma_GS_v6.2.txt", sep="\t", stringsAsFactors=F)
sets = lapply(sets, strsplit, " ")[[1]]

### read in gencode ids
gencode = fread("data/analysis_files/gencode.v26.GRCh38.protein_coding.gtf",header=F, sep=" ", data.table=F)[,c(7,9)]
gencode[,1] = sapply(strsplit(gencode$V7,"\\."), "[[", 1)
colnames(gencode) = c("ensid","name")

#### Enrichment function
enrichment = function(set, genes, gencode) {
	# reg="thalamus_proper"; dat=m.dat; quant=.01
	# set = sets[[which(sapply(sets,"[[",1)=="GO_bp:go_cellular_response_to_zinc_ion")]]
	
	### Extract gene set info
	set.id = set[1]
	set.genes = set[-1]
	
	### Get overlap
	all.genes = gencode$ensid
	in.set = all.genes %in% set.genes
	in.list = all.genes %in% genes
	
	### Test
	ft = fisher.test(table(in.list, in.set), alternative='greater')
	or = round(as.numeric(unlist(ft)[4]),2)
	p = as.numeric(unlist(ft)[1])
	
	### Format
	ens.ids = intersect(set.genes, genes)
	gene.names = gencode$name[match(ens.ids, gencode$ensid)]
	
	data.frame(id = set.id, or, p, ens.ids = paste(ens.ids, collapse=";"), gene.names = paste(gene.names, collapse=";"))
}


#### Function to run enrichment separately per endophenotype using the top 1% of genes
endo.enrich = function(endo, dat, sets, quant, gencode) {
	# endo = "thalamus_proper"; dat = m.dat; quant = .01
	# endo = "brain_substantia_nigra"; dat = e; quant = .01
	
	### Extract top genes
	quant.p = quantile(subset(dat, endoph == endo)$p, quant)
	genes = sapply(strsplit(unique(subset(dat, endoph == endo & p < quant.p)$locus),"\\."), "[[", 1)
	# print(head(genes))
	
	### Subset pathways to those that contain at least one of the significant genes
	## doing this as it is no point evaluating enrichment for gene sets where there are no relevant genes to begin with
	set.sub = sets[sapply(sets, function(x) any(x %in% genes))]
	
	### Analyse and format
	results = lapply(set.sub, function(x) enrichment(x, genes, gencode))
	r = do.call(rbind, results)
	if (!is.null(r)) {
		r = r[order(r$p),] # ordering the results according to significance (with the most sig gene sets on top)
	}
	return(r)
}



#########################################################
#### Perform enrichment separately per endophenotype ####
#########################################################


#### EQTL ####
e$endoph = e$tissue
exp.enrich = lapply(unique(e$endoph), function(x) { print(x); endo.enrich(x, e, sets, .01, gencode) })
names(exp.enrich) = unique(e$endoph)
n.eqtl.enrich = sum(sapply(exp.enrich, nrow))

#### SQTL ####
s$endoph = s$tissue
splic.enrich = lapply(unique(s$endoph), function(x) { print(x); endo.enrich(x, s, sets, .01, gencode) })
names(splic.enrich) = unique(s$endoph)
n.sqtl.enrich = sum(sapply(splic.enrich, nrow))

#### MRI ####
m.dat$endoph = sapply(strsplit(rownames(m.dat),"\\."),"[",1)
mri.enrich = lapply(unique(m.dat$endoph), function(x) { print(x); endo.enrich(x, m.dat, sets, .01, gencode)})
names(mri.enrich) = unique(m.dat$endoph)
n.mri.enrich = sum(sapply(mri.enrich, nrow)) # get n tests performed

#### NET ####
net$endoph = net$phen1
net.enrich = lapply(unique(net$endoph), function(x) { print(x); endo.enrich(x, net, sets, .01, gencode) })
names(net.enrich) = unique(net$endoph)
n.net.enrich = sum(sapply(net.enrich, nrow))



#### Check which are significant #####
enrich.alpha = .05 / (n.mri.enrich + n.eqtl.enrich + n.sqtl.enrich + n.net.enrich)
do.call(rbind,lapply(mri.enrich, function(x) x[x$p < enrich.alpha,]))
do.call(rbind,lapply(exp.enrich, function(x) x[x$p < enrich.alpha,]))
do.call(rbind,lapply(splic.enrich, function(x) x[x$p < enrich.alpha,]))
do.call(rbind,lapply(net.enrich, function(x) x[x$p < enrich.alpha,]))

## Save a table with sig results (since only brain regions were sig, we can ignore the others)
m.enrich.out = do.call(rbind,lapply(mri.enrich, function(x) x[x$p < enrich.alpha,]))
m.enrich.out$endoph = sapply(strsplit(rownames(m.enrich.out), "\\."),"[[",1) # add region column
m.enrich.out = m.enrich.out[,-which(colnames(m.enrich.out)=="ens.ids")] # remove ens ids
m.enrich.out = m.enrich.out[,c("endoph",colnames(m.enrich.out)[-5])] # reorder cols
m.enrich.out = m.enrich.out[c(3:15,1,16,2),] # reorder rows
m.enrich.out$p = signif(m.enrich.out$p, 2) # format p-values
write.table(m.enrich.out, "enrichment.significant.txt", row.names=F, quote=F, sep="\t")



# Save all other results
for (i in tissues) {
	write.table(exp.enrich[[i]], paste0('data/results/enrichment/eqtl/',i,".txt"), row.names=F, quote=F, sep="\t")
	write.table(splic.enrich[[i]], paste0('data/results/enrichment/sqtl/',i,".txt"), row.names=F, quote=F, sep="\t")
}
for (i in regions) { write.table(mri.enrich[[i]], paste0('data/results/enrichment/regions/',i,".txt"), row.names=F, quote=F, sep="\t") }
for (i in networks) { write.table(net.enrich[[i]], paste0('data/results/enrichment/networks/',i,".txt"), row.names=F, quote=F, sep="\t") }
