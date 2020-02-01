# functions from MEDIPS R package
# Modified to use limma instead of edgeR


MEDIPS.meth = function(
		MSet1 = NULL, 
		MSet2 = NULL, 
		CSet = NULL, 
		ISet1 = NULL, 
		ISet2 = NULL, 
		chr = NULL, 
		p.adj="bonferroni", 
		diff.method="edgeR", 
		CNV=FALSE,
		MeDIP=FALSE,
		minRowSum=10,
		diffnorm="tmm",
		batch=NULL,
		detRate=FALSE
		)
{
	nMSets1 = length(MSet1)	
	nMSets2 = length(MSet2)	
	nISets1 = length(ISet1)	
	nISets2 = length(ISet2)	
	if(is.list(CSet)) CSet=CSet[[1]]
	if(!is.list(MSet1)) MSet1=c(MSet1)
	if(!is.list(MSet2)) MSet2=c(MSet2)
	if(!is.list(ISet1)) ISet1=c(ISet1)
	if(!is.list(ISet2)) ISet2=c(ISet2)
	
	##Proof of correctness
	#######################
	if (!is.null(batch) & diff.method != "limma")
	    stop('Batch correction only implemented in limma')	

	if(MeDIP){
		if(class(CSet)!="COUPLINGset"){stop("You have to state a COUPLINGset object!")}
	}

   if (detRate & diff.method != "limma")
	    stop('Detection rate adjustment only implemented in limma')	



	controlSet = MSet1[[1]]

	for(i in 1:nMSets1){
		if(class(MSet1[[i]])!="MEDIPSset" & class(MSet1[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
		if(length(genome_count(MSet1[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
	}
	if(!is.null(MSet2)){
    		for(i in 1:nMSets2){
			if(class(MSet2[[i]])!="MEDIPSset" & class(MSet2[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(MSet2[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
    		}
	}	
	if(!is.null(ISet1)){
		for(i in 1:nISets1){
			if(class(ISet1[[i]])!="MEDIPSset" & class(ISet1[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(ISet1[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
		}
	}
	if(!is.null(ISet2)){
		for(i in 1:nISets2){
			if(class(ISet2[[i]])!="MEDIPSset" & class(ISet2[[i]])!="MEDIPSroiSet"){stop("You have to state a MEDIPSset or MEDIPSroiSet object!")}
			if(length(genome_count(ISet2[[i]]))!=length(genome_count(controlSet)))stop("MEDIPSset/MEDIPSroiSet objects are of different length!")
		}
	}


	##Data preparation
	###################
	##Calculate genomic coordinates
	##################################
	if(class(controlSet)=="MEDIPSset"){
		window_size = window_size(controlSet)
		no_chr_windows = ceiling(chr_lengths(controlSet)/window_size(controlSet))
		supersize_chr = cumsum(no_chr_windows)	
		GRanges.genome = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chr_names(controlSet), chr_lengths(controlSet), window_size(controlSet))
	}else if(class(controlSet)=="MEDIPSroiSet"){
		GRanges.genome = rois(controlSet)
		window_size=as.numeric(width(GRanges.genome))
	} 
	
	##Create data frame for all genomic windows and MEDIPS SETs
	if(MeDIP){
		base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), CF=genome_CF(CSet), stringsAsFactors=F)
	}else{
		base = data.frame(chr=as.vector(seqnames(GRanges.genome)), start=start(GRanges.genome), stop=end(GRanges.genome), stringsAsFactors=F)
	}

	rm(controlSet)
	gc()

	counts.medip = NULL
	rpkm.medip = NULL
	rms = NULL	
	counts.input = NULL
	rpkm.input = NULL
	
	##Add counts
	if(!is.null(MSet1)){
		for(i in 1:nMSets1){
			cat(paste("Preprocessing MEDIPS SET ", i, " in MSet1...\n", sep=""))
			counts.medip = cbind(counts.medip, MSet1=genome_count(MSet1[[i]]))
			rpkm.medip = cbind(rpkm.medip, round(((genome_count(MSet1[[i]])*10^9)/(window_size*number_regions(MSet1[[i]]))), digits=2))
			if(MeDIP){
				ccObj = MEDIPS.calibrationCurve(MSet=MSet1[[i]], CSet=CSet, input=F)
				rms = cbind(rms, MEDIPS.rms(MSet1[[i]], CSet, ccObj=ccObj))
			}		
		}
	}	
	if(!is.null(MSet2)){
		 for(i in 1:nMSets2){
			cat(paste("Preprocessing MEDIPS SET ", i, " in MSet2...\n", sep=""))
			counts.medip = cbind(counts.medip, MSet2=genome_count(MSet2[[i]]))
			rpkm.medip = cbind(rpkm.medip, round(((genome_count(MSet2[[i]])*10^9)/(window_size*number_regions(MSet2[[i]]))), digits=2))
			if(MeDIP){
				ccObj = MEDIPS.calibrationCurve(MSet=MSet2[[i]], CSet=CSet, input=F)
				rms = cbind(rms, MEDIPS.rms(MSet2[[i]], CSet, ccObj=ccObj))
			}
	  	 }
	}	
	if(!is.null(ISet1)){		
		for(i in 1:nISets1){
			cat(paste("Preprocessing INPUT SET ", i, " in ISet1...\n", sep=""))
			counts.input = cbind(counts.input, ISet1=genome_count(ISet1[[i]]))
			rpkm.input = cbind(rpkm.input, round(((genome_count(ISet1[[i]])*10^9)/(window_size*number_regions(ISet1[[i]]))), digits=2))
		   }
	}	
	if(!is.null(ISet2)){
		for(i in 1:nISets2){
			cat(paste("Preprocessing INPUT SET ", i, " in ISet2...\n", sep=""))
			counts.input = cbind(counts.input, ISet2=genome_count(ISet2[[i]]))
			rpkm.input = cbind(rpkm.input, round(((genome_count(ISet2[[i]])*10^9)/(window_size*number_regions(ISet2[[i]]))), digits=2))
		 	}
	}		
	
	##Extract data for selected chromosome
	#######################################
	if(!is.null(chr)){
		fi=base[,1]%in%chr
		
		cat("Extracting data for", chr, "...\n", sep=" ")
		if(length(fi)==0){stop("Stated chromosome does not exist in the COUPLING SET.")}
		if(!is.null(counts.medip)){
			counts.medip = counts.medip[fi,]
			rpkm.medip = rpkm.medip[fi,]
			rms = rms[fi,]
		}
		if(!is.null(counts.input)){
			counts.input = counts.input[fi,]
			rpkm.input = rpkm.input[fi,]
		}
		base = base[fi,]
		cat(nrow(base), "windows on", chr, "\n",sep=" ")		
	}
	
	##Set colnames and transform to data.frames
	############################################
	col.names.count = NULL
	col.names.rpkm = NULL
	col.names.rms = NULL
	if(nMSets1!=0){
		for(i in 1:nMSets1){
			col.names.count = c(col.names.count, paste(sample_name(MSet1[[i]]), ".counts", sep=""))
			col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet1[[i]]), ".rpkm", sep=""))
			if(MeDIP){
				col.names.rms = c(col.names.rms, paste(sample_name(MSet1[[i]]), ".rms", sep=""))
			}
		}		
	}
	if(nMSets2!=0){
		for(i in 1:nMSets2){
			col.names.count = c(col.names.count, paste(sample_name(MSet2[[i]]), ".counts", sep=""))
			col.names.rpkm = c(col.names.rpkm, paste(sample_name(MSet2[[i]]), ".rpkm", sep=""))
			if(MeDIP){
				col.names.rms = c(col.names.rms, paste(sample_name(MSet2[[i]]), ".rms", sep=""))
			}
		}		
	}
	if(nMSets1!=0 | nMSets2!=0){
		
		counts.medip = data.frame(counts.medip)
		colnames(counts.medip) = col.names.count
		rpkm.medip =  data.frame(rpkm.medip)
		colnames(rpkm.medip) = 	col.names.rpkm
		if(MeDIP){
			rms = data.frame(rms)
			colnames(rms) = col.names.rms
		}
	}		
	
	col.names.count.input = NULL
	col.names.rpkm.input = NULL
	if(nISets1!=0){
		for(i in 1:nISets1){
			col.names.count.input = c(col.names.count.input, paste(sample_name(ISet1[[i]]), ".counts", sep=""))
			col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet1[[i]]), ".rpkm", sep=""))
		}		
	}
	if(nISets2!=0){
		for(i in 1:nISets2){
			col.names.count.input = c(col.names.count.input, paste(sample_name(ISet2[[i]]), ".counts", sep=""))
			col.names.rpkm.input = c(col.names.rpkm.input, paste(sample_name(ISet2[[i]]), ".rpkm", sep=""))
		}		
	}
	
	if(nISets1!=0 | nISets2!=0){		
		
		counts.input = data.frame(counts.input)
		colnames(counts.input) = 	col.names.count.input
		rpkm.input = data.frame(rpkm.input)	
		colnames(rpkm.input) = 		col.names.rpkm.input
	}
		
	
	##If two groups of MEDIPS SETs are given 
	##calculate differential coverage
	##################################
	if(!is.null(MSet1) & !is.null(MSet2)){
		cat(paste("Differential coverage analysis...\n", sep=" "))
		
		##Correct for test selection if necessary
		if((nMSets1<3 | nMSets2<3) & diff.method=="ttest"){
			stop("Method 'ttest' is not valid for less than 3 replicates per group. Method 'edgeR' can be applied in this case.")
		}		
		
		if(diff.method=="edgeR"|diff.method=="limma"){
			
			if(diffnorm=="rpkm" | diffnorm=="rms"){
				stop("rpkm and rms normalization are not available for the edgeR mode. Here, differential enrichment is performed on the count data. Please set diffnorm to tmm or quantile.\n")
			}
			
			##Extract number of reads per sample
			if(!is.null(MSet1)){
				n.r.M1 = NULL
				for(i in 1:nMSets1){
					n.r.M1 = c(n.r.M1, number_regions(MSet1[[i]]))		
				}
			}	
			if(!is.null(MSet2)){
		 		n.r.M2 = NULL
		 		for(i in 1:nMSets2){
					n.r.M2 = c(n.r.M2, number_regions(MSet2[[i]]))	
	  	 		}
			}	

			#Quantile normalization
			if(diffnorm=="quantile"){
				cat("Performing quantile normalization on sequencing counts. Please note, the returned counts - but not the returned rpkm values - will be quantile normalized.\n")
				counts.medip.coln = colnames(counts.medip)				
				counts.medip = preprocessCore::normalize.quantiles(as.matrix(counts.medip), copy=FALSE)	
				counts.medip <- round(counts.medip)
				colnames(counts.medip) = counts.medip.coln
			}
		
		    if(diff.method=="edgeR"){
			  diff.results.list = MEDIPS.diffMeth(base=base, values=counts.medip, diff.method="edgeR", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, n.r.M1=n.r.M1, n.r.M2=n.r.M2, MeDIP=MeDIP, minRowSum=minRowSum, diffnorm=diffnorm)
			}else{
			  diff.results.list = MEDIPS.diffMeth(base=base, values=counts.medip, diff.method="limma", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, n.r.M1=n.r.M1, n.r.M2=n.r.M2, MeDIP=MeDIP, minRowSum=minRowSum, diffnorm=diffnorm, batch=batch, detRate=detRate)
			}
		}else if(diff.method=="ttest"){

			if(diffnorm=="quantile" | diffnorm=="tmm"){
				stop("Quantile and TMM normalization are not available for the t-test mode. Please chose edgeR or set diffnorm to rpkm or rms.\n")
			}

			if(diffnorm=="rpkm"){
				diff.results.list = MEDIPS.diffMeth(base=base, values=rpkm.medip, diff.method="ttest", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, MeDIP=MeDIP, minRowSum=minRowSum)
			}
			else if(diffnorm=="rms"){
				if(MeDIP){
					diff.results.list = MEDIPS.diffMeth(base=base, values=rms, diff.method="ttest", nMSets1=nMSets1, nMSets2=nMSets2, p.adj=p.adj, MeDIP=MeDIP, minRowSum=minRowSum)
				}
				else{
					stop("Invalid specification for parameter diffnorm because parameter MeDIP is FALSE (no rms values have been calculated).")
				}
			}
			else{
				stop("Unknown specification for parameter diffnorm.")
			}
		}else{
		stop("Selected method for calculating differential coverage not supported")
		}
		
		cat("Please note, log2 ratios are reported as log2(MSet1/MSet2).\n")
		diff.results = diff.results.list$diff.results
		diff.index = diff.results.list$diff.index
				
		rm(diff.results.list)
		gc()
	
	}else{
		cat("No differential coverage will be calculated- only one group of MEDIPS SETs given.\n")
	}
	
	##If two groups of INPUT SETs are given 
	##calculate CNV on mean per set
	##################################
	if(CNV){
		if(!is.null(ISet1) & !is.null(ISet2)){
			cat(paste("CNV analysis...\n", sep=" "))		
			cnv.combined = MEDIPS.cnv(base=base, rpkm.input=rpkm.input, nISets1=nISets1, nISets2=nISets2)
		}
		else{
			cat("Cannot perform CNV analysis- please specify two groups of INPUT SETs!\n")
		}
	}
	
	##Create results table	
	##################################
	cat(paste("Creating results table...\n", sep=" "))
	if(!is.null(counts.medip)){
		if(MeDIP){
			results = data.frame(base, counts.medip, rpkm.medip, rms, stringsAsFactors=F)
		}
		else{
			results = data.frame(base, counts.medip, rpkm.medip, stringsAsFactors=F)
		}
	}
	if(!is.null(counts.input)){
		if(!is.null(counts.medip))
		{
			results = data.frame(results, counts.input, rpkm.input, stringsAsFactors=F)
		}
		else{
			results = data.frame(base, counts.input, rpkm.input, stringsAsFactors=F)
		}
	}
	
	##Add mean counts, rpkm, rms columns
	set1idx=1:(nMSets1)
	counts.mean.C=numeric(dim(counts.medip)[1])
	rpkm.mean.C=numeric(dim(rpkm.medip)[1])
	for (i in set1idx){
		counts.mean.C=counts.mean.C+counts.medip[,i]
		rpkm.mean.C=rpkm.mean.C+rpkm.medip[,i]
	}
	counts.mean.C=counts.mean.C/nMSets1
	rpkm.mean.C=rpkm.mean.C/nMSets1
	if(MeDIP){
		rms.mean.C =numeric(dim(rms)[1])
		for (i in set1idx){
			rms.mean.C=rms.mean.C+rms[,i]
		}
		rms.mean.C=rms.mean.C/nMSets1
		results = data.frame(results, MSets1.counts.mean=counts.mean.C, MSets1.rpkm.mean=rpkm.mean.C, MSets1.rms.mean=rms.mean.C, stringsAsFactors=F)
		rm(counts.mean.C,rpkm.mean.C,set1idx,rms.mean.C)
	}else{
		results = data.frame(results, MSets1.counts.mean=counts.mean.C, MSets1.rpkm.mean=rpkm.mean.C, stringsAsFactors=F)
		rm(counts.mean.C,rpkm.mean.C,set1idx)
	}

	if(nMSets2>0){
		set2idx=(nMSets1+1):(nMSets1+nMSets2)
		counts.mean.T=numeric(dim(counts.medip)[1])
		rpkm.mean.T=numeric(dim(rpkm.medip)[1])
		for (i in set2idx){
			counts.mean.T=counts.mean.T+counts.medip[,i]
			rpkm.mean.T=rpkm.mean.T+rpkm.medip[,i]
		}	
		counts.mean.T=counts.mean.T/nMSets2
		rpkm.mean.T=rpkm.mean.T/nMSets2

		if(MeDIP){
			rms.mean.T =numeric(dim(rms)[1])
			for (i in set2idx){
				rms.mean.T=rms.mean.T+rms[,i]
			}
			rms.mean.T=rms.mean.T/nMSets2
			results = data.frame(results, MSets2.counts.mean=counts.mean.T, MSets2.rpkm.mean=rpkm.mean.T, MSets2.rms.mean=rms.mean.T, stringsAsFactors=F)
			rm(counts.mean.T,rpkm.mean.T,set2idx,rms.mean.T)
		}
		else{
			results = data.frame(results, MSets2.counts.mean=counts.mean.T, MSets2.rpkm.mean=rpkm.mean.T, stringsAsFactors=F)
			rm(counts.mean.T,rpkm.mean.T,set2idx)
		}
	}

	if(nISets1>1){
		setI1idx=1:(nISets1)
		counts.input.mean.C = counts.input[,setI1idx[1]]
		rpkm.input.mean.C = rpkm.input[,setI1idx[1]]
		for (i in setI1idx[-1]){
		  counts.input.mean.C =counts.input.mean.C+counts.input[,i]
		  rpkm.input.mean.C =rpkm.input.mean.C+rpkm.input[,i]
		}
		counts.input.mean.C =counts.input.mean.C/nISets1
		rpkm.input.mean.C =rpkm.input.mean.C/nISets1
		results = data.frame(results, ISets1.counts.mean=counts.input.mean.C, ISets1.rpkm.mean=rpkm.input.mean.C, stringsAsFactors=F)
		rm(counts.input.mean.C,rpkm.input.mean.C,setI1idx)
	}
	if(nISets2>1){
		setI2idx=(nISets1+1):(nISets1+nISets2)
		counts.input.mean.T = counts.input[,setI2idx[1]]
		rpkm.input.mean.T = rpkm.input[,setI2idx[1]]
		for (i in setI2idx[-1]){
		  counts.input.mean.T =counts.input.mean.T+counts.input[,i]
		  rpkm.input.mean.T =rpkm.input.mean.T+rpkm.input[,i]
		}
		counts.input.mean.T =counts.input.mean.T/nISets2
		rpkm.input.mean.T =rpkm.input.mean.T/nISets2

		results = data.frame(results, ISets2.counts.mean=counts.input.mean.T, ISets2.rpkm.mean=rpkm.input.mean.T, stringsAsFactors=F)
		rm(counts.input.mean.T,rpkm.input.mean.T,setI2idx)
	}
	
	if(MeDIP){rm(base, counts.medip, rpkm.medip, rms, counts.input, rpkm.input)}else{rm(base, counts.medip, rpkm.medip, counts.input, rpkm.input)}
	gc()
	
	##Add diff.meth results
	if(nMSets1!=0 & nMSets2!=0){
		cat(paste("Adding differential coverage results...\n", sep=" "))
		dummy.results = matrix(ncol=ncol(diff.results), nrow=nrow(results))
		
		if(diff.method=="edgeR"|diff.method=="limma"){
			c.names = colnames(diff.results)
			diff.results <- matrix(unlist(diff.results), ncol=ncol(diff.results), byrow=FALSE)
			colnames(diff.results)=c.names
			rm(c.names)
		}
		
		dummy.results[diff.index,] = diff.results
		colnames(dummy.results)=colnames(diff.results)		
		results = data.frame(results, dummy.results, stringsAsFactors=F)
			
		rm(diff.results, dummy.results, diff.index)
		gc()
	}
	
	##Add CNV results
	if(!is.null(ISet1) & !is.null(ISet2)){
		if(CNV){
			cat(paste("Adding CNV results...\n", sep=" "))
			dummy.results = matrix(ncol=1, nrow=(nrow(results)))
			for(i in 1:nrow(cnv.combined)){
				dummy.results[cnv.combined[i,1]:cnv.combined[i,2]] = cnv.combined[i,3]		
			}
		colnames(dummy.results)="CNV.log2.ratio"
		results = data.frame(results, dummy.results, stringsAsFactors=F)
		
		rm(dummy.results)	
		gc()
		}
	}
	
	rownames(results) = seq(1, nrow(results))
    	gc()
	return(results)
		
}
















MEDIPS.diffMeth = function(base = NULL, values=NULL, diff.method="ttest", nMSets1=NULL, nMSets2=NULL, n.r.M1=n.r.M1, n.r.M2=n.r.M2, p.adj="bonferroni", MeDIP, minRowSum=10, diffnorm="tmm", batch=NULL,
	detRate = FALSE)
{
	##edgeR##
	#########
	if(diff.method=="edgeR"){		
			
		##Extract non-zero MeDIP count windows
		cat(paste("Extracting count windows with at least", minRowSum, "reads...\n", sep=" "))
		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor rows
		if(MeDIP){
			cat(paste("Extracting non-zero coupling factor windows...\n", sep=" "))
			filter=filter & base[,4]!=0
		}

		cat(paste("Execute edgeR for count data of", sum(filter), "windows...\n", sep=" "))
		cat("(Neglecting parameter 'type')\n")

		cat("Creating a DGEList object...\n")						
		edRObj.group=c(rep(1, nMSets1), rep(2, nMSets2))
		edRObj.length=c(n.r.M1, n.r.M2)
		#d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group, lib.size=edRObj.length)
		d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group)

		if(diffnorm=="tmm"){
			cat("Apply trimmed mean of M-values (TMM) for library sizes normalization...\n")	
			d=edgeR::calcNormFactors(d, refColumn = 1)
		}
		if(diffnorm=="quantile" | diffnorm=="none"){
			cat("Skipping trimmed mean of M-values (TMM) library size normalization...\n")	
			d=edgeR::calcNormFactors(d, method="none")
		}
	    if(diffnorm!="tmm" & diffnorm!="quantile" & diffnorm!="none"){
	    	stop("diffnorm method unknown.")
	    }

		if(nMSets1!=1 | nMSets2!=1){
			cat("Estimating common dispersion...\n")			
			d=edgeR::estimateCommonDisp(d)
			cat("Estimating tagwise dispersion...\n")	
			d=edgeR::estimateTagwiseDisp(d)
			cat("Calculating differential coverage...\n")
			de.com=edgeR::exactTest(d,pair=c("2","1"))
		}else{
			cat("There is no replication, setting dispersion to bcv^2 where bcv=0.01.\n")
			cat("Please consider http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf section 2.9.\n")
			bcv = 0.01
			cat("Calculating differential coverage...\n")
			de.com = suppressWarnings(edgeR::exactTest(d, dispersion=bcv^2, pair=c("2","1")))
		}
		
		##Adjusting p.values for multiple testing
		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))
		colnames(de.com$table) = c("edgeR.logFC", "edgeR.logCPM", "edgeR.p.value")		
		diff.results = cbind(de.com$table, edgeR.adj.p.value=p.adjust(de.com$table$edgeR.p.value, p.adj)) 
		
		rm(de.com, edRObj.group, edRObj.length, d)
		gc()
	}else if (diff.method=="limma"){

	    ##Extract non-zero MeDIP count windows
		cat(paste("Extracting count windows with at least", minRowSum, "reads...\n", sep=" "))
		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor rows
		if(MeDIP){
			cat(paste("Extracting non-zero coupling factor windows...\n", sep=" "))
			filter=filter & base[,4]!=0
		}

		cat(paste("Execute limma for count data of", sum(filter), "windows...\n", sep=" "))
		cat("(Neglecting parameter 'type')\n")

		cat("Creating a DGEList object...\n")						
		edRObj.group=c(rep(1, nMSets1), rep(2, nMSets2))
		edRObj.length=c(n.r.M1, n.r.M2)
		#d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group, lib.size=edRObj.length)

		d <- edgeR::DGEList(counts = values[filter,], group = edRObj.group)			


		if(diffnorm=="tmm"){
			cat("Apply trimmed mean of M-values (TMM) for library sizes normalization...\n")	
			d=edgeR::calcNormFactors(d, refColumn = 1)
		}
		if(diffnorm=="quantile" | diffnorm=="none"){
			cat("Skipping trimmed mean of M-values (TMM) library size normalization...\n")	
			d=edgeR::calcNormFactors(d, method="none")
		}
	    if(diffnorm!="tmm" & diffnorm!="quantile" & diffnorm!="none"){
	    	stop("diffnorm method unknown.")
	    }

	    if(detRate)
	    	pctZero <- colMeans2(values==0)

		if(nMSets1!=1 | nMSets2!=1){
			edRObj.group <- as.factor(edRObj.group)
			if (!is.null(batch)){
              batch <- as.factor(batch)
              if(detRate){
               design <- model.matrix(~edRObj.group + batch + pctZero)
              }else{
			   design <- model.matrix(~edRObj.group + batch)
			  }

			}else{
              if(detRate){
               design <- model.matrix(~edRObj.group + pctZero)
              }else{
			   design <- model.matrix(~edRObj.group)
			  }
			}

			#cat("Fitting limma-trend...\n")			
		  	#logCPM <- edgeR::cpm(d, log=TRUE, prior.count=1)
            #de.com <- limma::lmFit(logCPM, design)
            #de.com <- limma::eBayes(de.com, trend=TRUE)
           
            cat("Fitting limma-voom...\n")	
            v <- limma::voom(d, design)
            de.com <- limma::lmFit(v, design)
            de.com <- limma::eBayes(de.com)

            tab <- limma::topTable(de.com, number = nrow(d), 
           	  coef=2, sort.by = "none")[,c(1,2,4)]

		}else{
			cat("There is no replication, setting dispersion to bcv^2 where bcv=0.01.\n")
			cat("Please consider http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf section 2.9.\n")
			bcv = 0.01
			cat("Calculating differential coverage...\n")
			de.com = suppressWarnings(edgeR::exactTest(d, dispersion=bcv^2, pair=c("2","1")))
		}
		
		##Adjusting p.values for multiple testing
		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))
		colnames(tab) = c("limma.logFC", "limma.logCPM", "limma.p.value")		
		diff.results = cbind(limma::topTable(de.com, number = nrow(d), 
			coef=2, sort.by = "none"), 
		  limma.adj.p.value=p.adjust(tab$limma.p.value, p.adj)) 
		
		rm(de.com, edRObj.group, edRObj.length, d)
		gc()

	}else{			
		##Extract non-zero MeDIP rows
		cat(paste("Extracting count windows with at least",minRowSum," reads...\n", sep=" "))

		filter= rowSums(values)>=minRowSum

		##Extract non-zero coupling factor windows
		if(MeDIP){
			cat(paste("Extracting non-zero coupling factor windows...\n", sep=" "))
			filter=filter & base[,4]!=0
		}
	
		#filter NA
		#if there is a na in a row, this row is sorted out?!
		filter=filter & !is.na(rowSums(values))

			
		cat(paste("Calculating score for", sum(filter), "windows...\n", sep=" "))
		ms1=1:nMSets1
		ms2=(nMSets1+1):(nMSets1+nMSets2)
		##Calculate ratios##
		####################
		ratio=rowSums(values[filter,ms1]+0.1)/rowSums(values[filter,ms2,drop=F]+0.1)*(nMSets2/nMSets1)
	
		##Calculate p.values##
		######################
		##Check for constant entries
		## Filter out constant entries, as they have no sd

		const = apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms1,drop=F],MARGIN=1,FUN=max) 	+
			apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=min) - apply(X=values[filter,ms2,drop=F],MARGIN=1,FUN=max) 	== 0

		t.test.p.value = rep(NA, sum(filter))				

		t.test.p.value[!const]=matTtest(values[filter,][!const,],groups=c(rep(1,nMSets1),rep(2,nMSets2)))$p.value

		##Calculate the final score##
		#############################
		score = (-log10(t.test.p.value)*10)*log(ratio)
		
		##Adjusting p.values for multiple testing
		cat(paste("Adjusting p.values for multiple testing...\n", sep=" "))	
			
		diff.results = cbind(score.log2.ratio=log2(ratio), score.p.value=t.test.p.value, score.adj.p.value=p.adjust(t.test.p.value, p.adj), score=score)
			
		rm(const, ratio, t.test.p.value, score)
		
	}		
	return(list(diff.results=diff.results, diff.index=which(filter)))
}