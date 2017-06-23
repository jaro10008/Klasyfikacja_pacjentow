library('MALDIquant')
library('MALDIquantForeign')
library('tools')
#library('matrixStats')

findPeaks <- function(object, threshold = 0.0){
	intensity <- object@intensity
	mass <- object@mass

	indices <- which(intensity > threshold)

	mass <- mass[indices]
	intensity <- intensity[indices]

	return(createMassPeaks(mass, intensity,snr=rep.int(NA_real_, length(mass)),metaData=list()))
}

anova_test <- function(dirname,threshold = 0.01,pThresh = 0.05, controlMin = -1, controlMax = -1, minMZ = 0.0, maxMZ = 100000.0){
	files = list.files(dirname) #list of files
	
	sort(files)	#filenames must be in an alphabetic order to select groups

	i=1

	data_set = list()

	for (filename in files) {
		exten = file_ext(filename)
		if(all(exten == 'csv')) {
			#CSV file
			table = read.csv(paste(dirname,filename,sep='/'),sep=';')
			spektrum = createMassSpectrum(mass=table[,1],intensity=table[,2],metaData=list(name=filename))
		}
		else if(all(exten == 'mzML')) {
			#MZML file
			spektrum = importMzMl(paste(dirname,filename,sep='/'))[[1]]
		}
		else{
			#unkown file type
			print(paste('ERROR (',filename,'): Unkown file type!',sep=''))
			return()
		}
		
		#progowo = which(spektrum@intensity > 300.0)

		#spektrum@mass = spektrum@mass[progowo]
		#spektrum@intensity = spektrum@intensity[progowo]

		mz = spektrum@mass

		intensity = spektrum@intensity

		s = length(intensity)	

		#Normalization
		if((controlMin < 0) | (controlMax < 0)){
			#Normalize total ion current
			norm_value = sum(diff(mz)*(intensity[2:s]+intensity[1:(s-1)]))/2 #Integration, total ion current
			#norm_value = totalIonCurrent(spektrum)
		}
		else{
			#Normalize indicator intensity
			subtable = which(mz>=controlMin & mz<=controlMax);
			norm_value = max(intensity[subtable])
		}

		if(norm_value==0){
			print(paste('ERROR (',filename,'): Normalization value is equal to 0!',sep=''))
			return()
		}

		checkedIndices = which(mz>=minMZ & mz <= maxMZ)	
		spektrum@intensity=spektrum@intensity[checkedIndices]/norm_value
		spektrum@mass=spektrum@mass[checkedIndices]

		# noise <- estimateNoise(spektrum,method="SuperSmoother")

		# intensity <- spektrum@intensity - noise[,2]
		# intensity[intensity<0] <- 0
		# spektrum@intensity <- intensity

		#peaks <- detectPeaks(spektrum, method="SuperSmoother",halfWindowSize=halfWindowSize2, SNR=snr)

		peaks <- findPeaks(spektrum,threshold) #own function finding peaks	

		peaks <- monoisotopicPeaks(peaks,minCor=0.95,tolerance=1e-4,distance=1.002,size=3L:10L) #deizotoping
		
		data_set[[i]] = peaks #add found peaks to the list

		if(i>1){
			#if(i>2) break
			#dev.new()
		}


		#print(peaks)
		#plot(spektrum)
		#points(peaks, col="red", pch=4)
		
		i=i+1
	}


	#Alignment
	peaks <- binPeaks(data_set,tolerance=0.0002,method='relaxed')	

	#Compute matrix
	matrixPeaks = intensityMatrix(peaks) #matrixPeaks - pierwszy indeks - plik, drugi indeks - masa
	matrixPeaks[is.na(matrixPeaks)] <- 0

	#preparing data for anova test(grouping data)
	groups = list()
	groupsAvg = list()
	groupsSd = list()

	i = 1
	prev_char = ""
	group = 0

	bIn <- 0


	for(filename in files){
		char = substr(filename,1,1)
		if(char != prev_char){
			if(group>0){
				groupsAvg[[group]] = colMeans(matrixPeaks[bIn:(i-1),])
				groupsSd[[group]] = apply(matrixPeaks[bIn:(i-1),],2,sd)

			}
			group = group + 1
			bIn = i
		}
		groups[[i]] = group
		prev_char = char
		i = i+1
	}
	groupsAvg[[group]] = colMeans(matrixPeaks[bIn:(i-1),])
	groupsSd[[group]] = apply(matrixPeaks[bIn:(i-1),],2,sd)
	

	group = unlist(groups)
	chosen = list()	
	wsk = 1

	for(i in 1:length(matrixPeaks[1,])){
		intensity = as.vector(matrixPeaks[,i])
		groupSymbols = substr(files,1,1)
		df = data.frame(groupSymbols,intensity)
		res = aov(intensity~groupSymbols,data=df)
		pvalue = summary(res)[[1]][["Pr(>F)"]]
		if(pvalue<pThresh){
			chosen[[wsk]] = i
			wsk = wsk+1
		}
	}


	macierz = matrixPeaks[,unlist(chosen)]
	macierz = t(macierz)

	#macierz - pliki w kolumnach, masy w wierszach


	#HEATMAP
	#dev.new()
	pdf(sprintf("heatmap_%f.pdf",threshold), max(7,length(files)/5), max(5,length(chosen)/5))

	#heatmapMatrix = t(scale(t(macierz)))
	heatmapMatrix = t(apply(macierz, 1, function(x)(x-min(x))/(max(x)-min(x))))
	heatmap(heatmapMatrix, col = colorRampPalette(c("green","yellow","red"))(256), scale="none", margins=c(5,10), labCol=files)
	dev.off()


	#BARPLOT
	#dev.new()
	pdf(sprintf("barplot_%f.pdf",threshold), max(7,length(chosen)/2), 5)
	#groupsAvg = matrix(unlist(groupsAvg))
	groupsAvg=do.call(rbind,groupsAvg)
	groupsSd=do.call(rbind,groupsSd)
	groupNames = unique(substr(files,1,1))
	groupsAvg = groupsAvg[,unlist(chosen)]
	groupsSd = groupsSd[,unlist(chosen)]
	groupsAvgScaled = apply(t(groupsAvg), 1, function(x)(x)/(max(x)))

	barPlot = barplot(groupsAvgScaled, beside = TRUE, legend.text = groupNames, 
        	args.legend = list(x = "topleft", bty="n"),col=rainbow(length(groupNames)),names.arg=round(attributes(matrixPeaks)$mass[unlist(chosen)],digits=2),las=2)
	title("Grouped barplot of chosen masses");
	# arrows(barPlot, groupsAvg - groupsSd, barPlot,
       #	groupsAvg + groupsSd, lwd = 0.1, angle = 90,
       #	code = 3, length = 0.005)

	dev.off()


	#PCA ANALYSIS
	PCAMatrix = t(apply(macierz, 1, function(x)(x)/(max(x))))
	# print(PCAMatrix)
	PCA = prcomp(PCAMatrix,center=TRUE,scale.=TRUE)
	PCA = PCA$rotation[,1:2]
	#dev.new()
	pdf(sprintf("pca_%f.pdf",threshold), 7, 5)
	par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

	plot(PCA,col=group)

	legend("topright",groupNames,col=unique(group),pch=rep(c(16,18),each=4),inset=c(-0.2,0))
	title("Results of PCA analysis");
	dev.off()

}
args = commandArgs(trailingOnly=TRUE)
len = length(args);
#dirname,threshold,pThresh, controlMin, controlMax, minMZ, maxMZ

if(len<2){
	print("Too few arguments!");
} else if(len==2) {
	anova_test(args[1],as.numeric(args[2]));
} else if(len==3){
	anova_test(args[1],as.numeric(args[2]),as.numeric(args[3]));
} else if(len==5){
	anova_test(args[1],as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]));
} else if(len==7){
	anova_test(args[1],as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]));
} else{
	print('Bad number of arguments');
}

#anova_test('folder',0.0001,0.05,-1,-1,0,10000)