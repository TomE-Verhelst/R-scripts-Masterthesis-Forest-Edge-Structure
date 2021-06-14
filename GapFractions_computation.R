##loading in packages
library("lidR")
library("rgdal")

##resolution
resol <- 0.25 #ENTER prefered resolution (high increase in runtime below 0.25)

##loading in pointcloud files
#for the script to work, the loaded pointclouds need to be aligned with the
#axes of the coordinate system
setwd() #ENTER working directory
L1 <- 'BeLoT1_dev12.las'
L3 <- 'BeLoT3_dev12.las'
M1 <- 'BeMeT1_dev12.las'
M3 <- 'BeMeT3_dev12.las'
H1 <- 'BeHiT1_dev12.las'
H3 <- 'BeHiT3_dev12.las'
Fr1 <- 'NFT1_dev12.las'
Fr2 <- 'NFT2_dev12.las'

transects <- list(L1,L3,M1,M3,H1,H3,Fr1,Fr2)

##loading other data
setwd("C:/Users/vtom1/OneDrive/Documenten/DataverwerkingMT/LiDAR")
canophgrad <- read.table("canopy_height_gradients1205.txt",header = TRUE)
attach(canophgrad)
transectgrads <- list(BeLoT1z,BeLoT3z,BeMeT1z,BeMeT3z,BeHiT1z,
                    BeHiT3z,NFT1z,NFT2z)
canoph <- read.table("summary_hgrad_zmax.txt",header = TRUE)
attach(canoph)

##Defining gapfrac() function; calculates the proportion of non NA values in
#a input matrix m
gapfrac <- function(m){
	dimension <- dim(m)
	r <- dimension[1]
	c <- dimension[2]
	size <- length(m)
	bernou <- matrix(1,nrow=r,ncol=c)
	for (i in 1:r){
		for (j in 1:c){
			if (is.na(m[i,j])==TRUE){
				bernou[i,j] <- 0		#switch to one if you prefer
			}					#the proportion of NA values
		}
	}
	nNA <- sum(bernou)
	proportionNA <- nNA/size	
}

##preallocation of the gap fraction vectors
gapfraction_canopy <- rep(0,8)
gapfraction_edge <- rep(0,8)

##main loop
t <- 1
for (transect in transects){
	#reading in pointcloud files and applying the voxel_metrics function
	LASfile <- transect
	las <- readLAS(LASfile, select = "xyz")
      vm <- voxel_metrics(las, ~length(Z), res = c(resol), all_voxels = TRUE)
	#voxel_metrics creates 3D voxel grid
	attach(vm)	#four vectors: X,Y,Z and V1 (V1 contains number of data
			#points in the voxels at postion x,y,z
	
	#creating 3D matrix
	xdiff <- max(X)-min(X)
	ydiff <- max(Y)-min(Y)
	zdiff <- max(Z)-min(Z)
	nxp <- xdiff/resol+1
	nyp <- ydiff/resol+1
	nzp <- zdiff/resol+1
	transect_matrix <- array(rep(0,length(V1)), dim=c(nzp,nxp,nyp))
	for (i in 1:nxp){
		for (j in 1:nyp){
			zvalues <- rev(V1[(1+(j-1)*nzp+(i-1)*nyp*nzp):
					      (nzp+(j-1)*nzp+(i-1)*nyp*nzp)])
			transect_matrix[,i,j] <- zvalues
		}
	}	#transect_matrix is the 3D matrix of the transect (values are PCs)
	
	#canopy gap fraction
	hcanop <- canop_height[t]
	vec <- unlist(transectgrads[t])
	loc <- 1
	height <- vec[loc]
	while (height < (hcanop*0.95)){	#lever 1
		loc <- loc + 1
		height <- vec[loc]
	}
	ycutoff <- loc
	ypoint <- (ycutoff*0.05)/resol
	zcutoff <- max(Z)-(round(hcanop,digits=0))+10
	zrow <- zcutoff/resol
	inter_forest <- transect_matrix[1:zrow,,ypoint:nyp]
	gapprojection <- matrix(0,nrow=nxp,ncol=(nyp-ypoint+1))
	for (xloc in 1:nxp){
		for (yloc in 1:(nyp-ypoint+1)){
			gapprojection[xloc,yloc] <- sum(inter_forest[,xloc,yloc],na.rm=TRUE)
			if (gapprojection[xloc,yloc]==0){
				gapprojection[xloc,yloc] <- NA
			}
		}
	}
	gapfraction_canopy[t] <- 1-gapfrac(gapprojection)
	
	#edge gap fraction
	forest_edge <- transect_matrix[,,1:ypoint]
	edgeprojection <- matrix(0,nrow=nzp,ncol=nxp)
	for (xloc in 1:nxp){
		for (zloc in 1:nzp){
			edgeprojection[zloc,xloc] <- sum(forest_edge[zloc,xloc,],na.rm=TRUE)
			if (edgeprojection[zloc,xloc]==0){
				edgeprojection[zloc,xloc] <- NA
			}
		}
	}
	overheight <- max(Z)-(hcanop*0.75)		#lever 2
	ohp <- overheight/resol
	edgeprojection <- edgeprojection[-1:-ohp,]
	if (min(Z)<0){
		underground <- 0-min(Z)
		undergr_points <- underground/resol
		depth <- dim(edgeprojection)[1]
		edgeprojection <- edgeprojection[-depth:-(depth-undergr_points),]
	}
	gapfraction_edge[t] <- 1-gapfrac(edgeprojection)				

	#next transect
	t <- t + 1
}

