# TODO: Add comment
# 
# Author: tiffnk
###############################################################################


files <- list.files("/Users/tiffnk/Google Drive/Kaweah/data/", full.names=TRUE, recursive=FALSE, pattern=".txt")

kaweahraw <- vector("list", length(files))
for(i in 1:length(kaweahraw)){
	kaweahraw[[i]] <- read.table(files[[i]], sep="\t", header=TRUE)
	names(kaweahraw[[i]])[[4]] <- c("discharge_cfs")
	names(kaweahraw)[[i]] <- kaweahraw[[i]]$site_no[[1]]
}

cdecfiles <- list.files("/Users/tiffnk/Google Drive/Kaweah/data/cdec", full.names=TRUE )
cdecnames <- list.files("/Users/tiffnk/Google Drive/Kaweah/data/cdec", full.names=FALSE )
cdecnames <- c(strsplit(cdecnames,"_")[[1]][[1]],strsplit(cdecnames,"_")[[2]][[1]])
for(i in 1:length(cdecfiles)){
	kaweahraw[[length(files)+i]] <- read.csv(cdecfiles[[i]])
	names(kaweahraw)[[length(files)+i]] <- cdecnames[[i]]
	names(kaweahraw[[length(files)+i]]) <- c("datetime","PST","discharge_cfs")
	for(j in 1:length(kaweahraw[[length(files)+i]]$datetime)){
		kaweahraw[[length(files)+i]]$datetime[[j]] <- paste(paste0(strsplit(as.character(kaweahraw[[length(files)+i]]$datetime[[j]]),"")[[1]][1:4],collapse=""),paste0(strsplit(as.character(kaweahraw[[length(files)+i]]$datetime[[j]]),"")[[1]][5:6],collapse=""),paste0(strsplit(as.character(kaweahraw[[length(files)+i]]$datetime[[j]]),"")[[1]][7:8],collapse=""),sep="-")
	}
	kaweahraw[[length(files)+i]]$PST <- NULL
	kaweahraw[[length(files)+i]]$site_no <- cdecnames[[i]]
	kaweahraw[[length(files)+i]]$agency_cd <- "CDEC"
	kaweahraw[[length(files)+i]]$discharge_cfs[which(kaweahraw[[length(files)+i]]$discharge_cfs=="m")] <- NA
	kaweahraw[[length(files)+i]]$discharge_cfs <- as.numeric(as.character((kaweahraw[[length(files)+i]]$discharge_cfs)))
	kaweahraw[[length(files)+i]]$discharge_cfs[which(kaweahraw[[length(files)+i]]$discharge_cfs<0)] <- 0
}

for(i in 1:length(kaweahraw)){
	kaweahraw[[i]]$datetime <- as.Date(kaweahraw[[i]]$datetime,"%Y-%m-%d")
}

kaweahraw$'11203220c' <- rbind.data.frame(kaweahraw$'11203200'[,c(1:4)], kaweahraw$'11203220'[,c(1:4)])

kall <- data.frame(date=seq.Date(from=as.Date("1900-01-01",format="%Y-%m-%d"),to=as.Date("2016-09-26",format="%Y-%m-%d"),by="day"))
for(i in 1:length(kaweahraw)){
	kall <- merge(kall, kaweahraw[[i]][,c("datetime","discharge_cfs")],by.x="date",by.y="datetime", all.x=TRUE)
	names(kall)[which(names(kall)=="discharge_cfs")] <- names(kaweahraw)[[i]]
}

kall$kw_final <- NA
kall$kw_final[1:which(kall$date==as.Date("1958-09-30",format="%Y-%m-%d"))] <- kall$'11210500'[1:which(kall$date==as.Date("1958-09-30",format="%Y-%m-%d"))] 
kall$kw_final[which(kall$date==as.Date("1958-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1990-09-30",format="%Y-%m-%d"))] <- kall$'11209900'[which(kall$date==as.Date("1958-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1990-09-30",format="%Y-%m-%d"))] + kall$'11210100'[which(kall$date==as.Date("1958-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1990-09-30",format="%Y-%m-%d"))] 
kall$kw_final[which(kall$date==as.Date("1990-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("2016-09-26",format="%Y-%m-%d"))] <- kall$trm[which(kall$date==as.Date("1990-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("2016-09-26",format="%Y-%m-%d"))]

kall$tul_final <- NA
kall$tul_final[1:which(kall$date==as.Date("1957-09-30",format="%Y-%m-%d"))] <- kall$'11204500'[1:which(kall$date==as.Date("1957-09-30",format="%Y-%m-%d"))] 
kall$tul_final[which(kall$date==as.Date("1957-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1991-09-30",format="%Y-%m-%d"))] <- kall$'11203220c'[which(kall$date==as.Date("1957-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1991-09-30",format="%Y-%m-%d"))] + kall$'11204500'[which(kall$date==as.Date("1957-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("1991-09-30",format="%Y-%m-%d"))] 
kall$tul_final[which(kall$date==as.Date("1991-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("2016-09-26",format="%Y-%m-%d"))] <- kall$scc[which(kall$date==as.Date("1991-10-1",format="%Y-%m-%d")):which(kall$date==as.Date("2016-09-26",format="%Y-%m-%d"))]



perc90 <- NA
for(i in 2:length(kall)){
	perc90[[i-1]] <- quantile(kall[,i],0.9, na.rm=TRUE)
}

final_data_cfs <- kall[,c(1,8,9,13,14)]
write.csv(final_data_cfs,"/Users/tiffnk/Google Drive/kaweah/final_data_cfs.csv")

finalperc90 <- data.frame(gauge=names(final_data_cfs)[2:5],perc90=NA)
finalperc90$perc90[[1]] <- quantile(final_data_cfs[2],0.9, na.rm=TRUE)
finalperc90$perc90[[2]] <- quantile(final_data_cfs[3],0.9, na.rm=TRUE)
finalperc90$perc90[[3]] <- quantile(kall$trm,0.9, na.rm=TRUE)
finalperc90$perc90[[4]] <- quantile(kall$scc,0.9, na.rm=TRUE)

write.csv(finalperc90,"/Users/tiffnk/Google Drive/kaweah/finalperc90.csv")
