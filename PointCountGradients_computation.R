##loading in packages
library("lidR")
library("rgdal")

##resolution
resol <- 0.25 #ENTER prefered resolution (big increase in runtime below 0.25)

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

##Defining gapfrac() function; calculates the proportion of non NA values in a input matrix m
gapfrac <- function(m){
	dimension <- dim(m)
	r <- dimension[1]
	c <- dimension[2]
	size <- length(m)
	bernou <- matrix(1,nrow=r,ncol=c)
	for (i in 1:r){
		for (j in 1:c){
			if (is.na(m[i,j])==TRUE){
				bernou[i,j] <- 0		#switch to one if you prefer the proportion of NA values
			}					#in the input matrix m				
		}
	}
	nNA <- sum(bernou)
	proportionNA <- nNA/size	
}

##preallocation of the PC gradient matrices
#PC gradient matrices contain 3 rows (3 layers) for each transect
n <- 8				 #ENTER number of transects
max_length <- 150              #ENTER transect length (y-axis in meters)
max_nyp <- max_length/resol+1  #calculates number of points along y-axis
PCP_matrix <- matrix(NA,nrow=(3*n),ncol=max_nyp)	#point count proportion
TPPC_matrix <- matrix(NA,nrow=(3*n),ncol=max_nyp)  	#total plane point count
APC_matrix <- matrix(NA,nrow=(3*n),ncol=max_nyp)	#averaged point count

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
	}		#transect_matrix is the 3D matrix of the transect (values are PCs)
	
	#defining layer height + layers
	hls <- 0.5/resol	#height of the different vegetation layers
	hhs <- 3.5/resol	#low shrub (ls):0-0.5 m/ high shrub (hs):0.5-4 m
	hla <- 4/resol	#low arboreal (la):4-8 m
	lowshrub <- matrix(0,nrow=hls,ncol=nxp)
	highshrub <- matrix(0,nrow=hhs,ncol=nxp)
	lowarbor <- matrix(0,nrow=hla,ncol=nxp)

	#preallocation of the gradient vectors
	PCP_vector_ls <- rep(NA,max_nyp)
	PCP_vector_hs <- rep(NA,max_nyp)
	PCP_vector_la <- rep(NA,max_nyp)
	TPPC_ls <- rep(NA,max_nyp)
	TPPC_hs <- rep(NA,max_nyp)
	TPPC_la <- rep(NA,max_nyp)
	APC_ls <- rep(NA,max_nyp)
	APC_hs <- rep(NA,max_nyp)
	APC_la <- rep(NA,max_nyp)

	#gradient creating loop
	#loop cycles through each location p along the y-axis,
	#then through each location q along the x-axis	
	for (p in 1:nyp){
		intersect <- transect_matrix[,,p] #xz-slice at y location p
		for (q in 1:nxp){
			#defines ground height at each individual location q
			h <- dim(transect_matrix)[1]	#h = lowest z-axis point
			while (is.na(intersect[h,q])==TRUE && h > 1){
				h <- h - 1
			}
			if ((h-hls-hhs-hla+1) < 10){   #fail safe
				lowshrub[,q] <- NA
				highshrub[,q] <- NA
				lowarbor[,q] <- NA
			} else {
				#filling in the slices for each layer
				lowshrub[,q] <- intersect[h:(h-hls+1),q]
				highshrub[,q] <- intersect[(h-hls):(h-hls-hhs+1),q]
				lowarbor[,q] <- intersect[(h-hls-hhs):(h-hls-hhs-hla+1),q]
			}
		}

		#creating PCP value at location p via the gapfrac function
		PCP_vector_ls[p] <- gapfrac(lowshrub)
		PCP_vector_hs[p] <- gapfrac(highshrub)
		PCP_vector_la[p] <- gapfrac(lowarbor)

		#creating TPPC value at location p via the sum function
		TPPC_ls[p] <- sum(lowshrub,na.rm=TRUE)
		TPPC_hs[p] <- sum(highshrub,na.rm=TRUE)
		TPPC_la[p] <- sum(lowarbor,na.rm=TRUE)

		#creating APC value at location p
		lowshrub[is.na(lowshrub)] <- 0	#NA values need to be changed
		highshrub[is.na(highshrub)] <- 0	#to 0 values to impact the
		lowarbor[is.na(lowarbor)] <- 0	#average
		meanplane_ls <- rowMeans(lowshrub)	#averaging along the x-axis
		meanplane_hs <- rowMeans(highshrub)
		meanplane_la <- rowMeans(lowarbor)
		APC_ls[p] <- sum(meanplane_ls) #summation
		APC_hs[p] <- sum(meanplane_hs)
		APC_la[p] <- sum(meanplane_la)
	}

	#filling the gradient matrices
	PCP_matrix[1+3*(t-1),] <- PCP_vector_ls
	PCP_matrix[2+3*(t-1),] <- PCP_vector_hs
	PCP_matrix[3+3*(t-1),] <- PCP_vector_la
	TPPC_matrix[1+3*(t-1),] <- TPPC_ls
	TPPC_matrix[2+3*(t-1),] <- TPPC_hs
	TPPC_matrix[3+3*(t-1),] <- TPPC_la
	APC_matrix[[1+3*(t-1),] <- APC_ls
	APC_matrix[[2+3*(t-1),] <- APC_hs
	APC_matrix[[3+3*(t-1),] <- APC_la
	
	#next transect
	t <- t + 1
}

##save gradient matrices as txt files
write.table(PCP_matrix,file="ENTER DIRECTORY/prefered_name.txt",
            sep = "\t",row.names = FALSE, col.names = FALSE)
write.table(TPPC_matrix,file="ENTER DIRECTORY/prefered_name.txt",
            sep = "\t",row.names = FALSE, col.names = FALSE)
write.table(APC_matrix,file="ENTER DIRECTORY/prefered_name.txt",
            sep = "\t",row.names = FALSE, col.names = FALSE)	
		






