##required packages
library(tuneR)
library(seewave)
library(matrixStats)
library(soundecology)

##set directory and load in sound files
#ENTER desired working directory below (dir)
dir <- ""						#directory to recordings
files <-list.files(path=dir, full.names = TRUE)		#list of recordings


##main loop
rec <- 1
n <- 168				#ENTER number of recordings in directory

#preallocation of the AI vectors
ACIvec <- rep(0, n)
ADIvec <- rep(0, n)
AEIvec <- rep(0, n)
BIvec <- rep(0, n)
Hvec <- rep(0, n)

#for loop that cycles through all recordings in the directory
for (file in files){
	w <- readWave(file)	#reads sound recording (tuneR package)
	#all AI functions are from the soundecology package
	ndi <- ndsi(w, fft_w = 512, anthro_min = 1000, anthro_max = 2000, 
	            bio_min = 2000, bio_max = 10000)
	NDSI <- ndi$ndsi_left		#selection of the wanted value
	if (NDSI < 0.7){
		ACIvec[rec] <- NA
		ADIvec[rec] <- NA
		AEIvec[rec] <- NA
		BIvec[rec] <- NA
		Hvec[rec] <- NA	
	} else{
		aci <- acoustic_complexity(w, max_freq = 10000, j = 5)
      	adi <- acoustic_diversity(w, max_freq = 10000, freq_step = 1000)
		aei <- acoustic_evenness(w, max_freq = 10000, freq_step = 1000)
		bi <- bioacoustic_index(w, min_freq = 2000, max_freq = 10000)
		ENT <- H(w, f = 48000)
		ACI <- aci$AciTotAll_left_bymin 	#selection of the wanted value
		ADI <- adi$adi_left	         	#selection of the wanted value
		AEI <- aei$aei_left		  	#selection of the wanted value
		BI <- bi$left_area		   	#selection of the wanted value
      	ACIvec[rec] <- ACI
		ADIvec[rec] <- ADI
		AEIvec[rec] <- AEI
		BIvec[rec] <- BI
		Hvec[rec] <- ENT
	}
	
	#next recording
	rec <- rec + 1
	print(rec)
}

##manually structuring AI values into a workable matrix
#this is based on the particular recording schedule in this thesis and needs
#to be adjusted accordingly (and to one's preferences)
Mm_ACI <- rbind(ACIvec[1:16],ACIvec[25:40],ACIvec[49:64],ACIvec[73:88],
                ACIvec[97:112],ACIvec[121:136],ACIvec[145:160])
Mm_ADI <- rbind(ADIvec[1:16],ADIvec[25:40],ADIvec[49:64],ADIvec[73:88],
                ADIvec[97:112],ADIvec[121:136],ADIvec[145:160])
Mm_AEI <- rbind(AEIvec[1:16],AEIvec[25:40],AEIvec[49:64],AEIvec[73:88],
                AEIvec[97:112],AEIvec[121:136],AEIvec[145:160])
Mm_BI <- rbind(BIvec[1:16],BIvec[25:40],BIvec[49:64],BIvec[73:88],
               BIvec[97:112],BIvec[121:136],BIvec[145:160])
Mm_H <- rbind(Hvec[1:16],Hvec[25:40],Hvec[49:64],Hvec[73:88],
              Hvec[97:112],Hvec[121:136],Hvec[145:160])

Me_ACI <- rbind(ACIvec[17:24],ACIvec[41:48],ACIvec[65:72],ACIvec[89:96],
                ACIvec[113:120],ACIvec[137:144],ACIvec[161:168])
Me_ADI <- rbind(ADIvec[17:24],ADIvec[41:48],ADIvec[65:72],ADIvec[89:96],
                ADIvec[113:120],ADIvec[137:144],ADIvec[161:168])
Me_AEI <- rbind(AEIvec[17:24],AEIvec[41:48],AEIvec[65:72],AEIvec[89:96],
                AEIvec[113:120],AEIvec[137:144],AEIvec[161:168])
Me_BI <- rbind(BIvec[17:24],BIvec[41:48],BIvec[65:72],BIvec[89:96],
               BIvec[113:120],BIvec[137:144],BIvec[161:168])
Me_H <- rbind(Hvec[17:24],Hvec[41:48],Hvec[65:72],Hvec[89:96],
              Hvec[113:120],Hvec[137:144],Hvec[161:168])

Mm <- rbind(Mm_ACI, Mm_ADI, Mm_AEI, Mm_BI, Mm_H)
Me <- rbind(Me_ACI, Me_ADI, Me_AEI, Me_BI, Me_H)

##saving AI values in a dawn and dusk text file
write.table(Mm, file = ".txt", sep = "\t",		#ENTER desired text file name (name.txt)
            row.names = FALSE, col.names = TRUE)
write.table(Me, file = ".txt", sep = "\t",		#ENTER desired text file name (name.txt)
            row.names = FALSE, col.names = TRUE)
