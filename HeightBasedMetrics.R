##loading in required packages
library("lidR")
library("rgdal")

##resolution
resol <- 0.05 #ENTER prefered resolution

##loading in pointcloud data
setwd()				#ENTER las file directory
L1 <- 'BeLoT1_dev12.las'	#ENTER correct las file names
L3 <- 'BeLoT3_dev12.las'
M1 <- 'BeMeT1_dev12.las'
M3 <- 'BeMeT3_dev12.las'
H1 <- 'BeHiT1_dev12.las'
H3 <- 'BeHiT3_dev12.las'
Fr1 <- 'NFT1_dev12.las'
Fr2 <- 'NFT2_dev12.las'

transects <- list(L1,L3,M1,M3,H1,H3,Fr1,Fr2)

##choose the desired analysis function

canoph <- function(z){
	height <- max(z)
	return(height)
}

vegh <- function(z){
	height <- max(z)-min(z)
	return(height)
}

terrainh <- function(z){
	height <- min(z)
	return(height)
}

##main script
t <- t + 1
height_gradient <- matrix(NA, nrow = 8, ncol = max_size)
for (transect in transects){
	LASfile <- transect
	las <- readLAS(LASfile, select = "xyz")
	vm <- grid_metrics(las, ~canoph(Z), res = resol)
	h <- as.matrix(vm)
		dimensions <- dim(grid)
	r <- dimensions[1]
	gradient <- rep(0,r) 
	for (n in 1:r){
		vector <- grid[n,]
		a <- median(vector,na.rm = TRUE)
      	gradient[n] <- a
	}
	height_gradient[t,1:r] <- rev(gradient)

	#next transect
	t <- t + 1
}

##saving height gradients in a text file
write.table(height_gradient, file = "C:/Users/vtom1/OneDrive/Documenten/DataverwerkingMT/LiDAR/vegetation_height_gradients1205.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE)
	
##defining mode function
getmode <- function(v){
	uniqv <- unique(v)
	len <- length(uniqv)
	n <- 1
	uniqv_edit <- 0
	for (value in 1:len){
		if (uniqv[value] > 5){
			uniqv_edit[n] <- uniqv[value]
			n <- n + 1
		}
	}
	print(uniqv_edit)
	uniqv_edit[which.max(tabulate(match(v, uniqv_edit)))]
}	



	
