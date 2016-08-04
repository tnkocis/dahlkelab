# TODO: Add comment
# 
# Author: tiffnk
###############################################################################
updateplots(new_file_path="/Users/tiffnk/Desktop/DNAoutput/T12_8_1_2016.csv",
		most_recent_master_dataset="/Users/tiffnk/Desktop/DNAoutput/MasterData_created_2016-08-03 18-09-45.RData", 
		plot_output_folderpath="/Users/tiffnk/Desktop/DNAoutput/",
		master_dataset_output_folderpath="/Users/tiffnk/Desktop/DNAoutput/",
		width=16,
		height=8)


### change height and width 
### possibly make all DNA separate plots
###check faceting variable
updateplots <- function(new_file_path, most_recent_master_dataset, plot_output_folderpath, master_dataset_output_folderpath, width, height){
	library(ggplot2)
	load(most_recent_master_dataset)
	filelocation <- new_file_path
	output <- read.csv(filelocation)
	filename <- basename(filelocation)
	
	if(filename%in%data_list$all$original_filename){
		print("file already in database!")
	}else{
		filedate <- as.Date(paste(strsplit(strsplit(filename, split="[.]")[[1]][1],"_")[[1]][2:4],collapse="-"), format="%m-%d-%Y")
		DNAtype <- paste(strsplit(strsplit(filename, split="[.]")[[1]][1],"_")[[1]][1])
		##define as fn argument
		start_quant <- 1000
		###
		output$Sample <- as.character(output$Sample)
		output$Starting.Quantity..SQ <- start_quant
		output$DNAtype <- DNAtype
		output$original_filename <- filename
		
		for(i in 1:length(output$Sample)){
			output$Sample_date[[i]] <- format(as.POSIXct(paste(strsplit(output$Sample[[i]],"_")[[1]][1]),format="%m/%d/%Y %H:%M"),format="%m/%d/%Y %H:%M")
			output$Sample_type[[i]] <- strsplit(output$Sample[[i]],"_")[[1]][2]
		}
		
		output$duplicate <- NA
		for(i in 1:length(output$Sample)){
			output$duplicate[[i]] <- strsplit(output$Sample[[i]],"_")[[1]][3]
			output$Sample[[i]] <- paste(strsplit(output$Sample[[i]],"_")[[1]][1:2], collapse="_")
			if(is.na(output$Cq[[i]])){
				output$Cq.edit[[i]] <- 40
			}else{
				output$Cq.edit[[i]] <- output$Cq[[i]]
			}
		}
		
		for(i in 1:length(output$Sample)){
			output$Cq.edit.Dup.Mean[[i]] <- mean(output$Cq.edit[which(output$Sample==output$Sample[[i]])])
			output$Cq.edit.Dup.Std.Dev[[i]] <- sd(output$Cq.edit[which(output$Sample==output$Sample[[i]])])
		}
		
		for(i in 1:length(output$Sample_type)){
			if(strsplit(output$Sample_type[[i]]," ")[[1]][1]=="LYS"){
				output$LYS_num[[i]] <- strsplit(output$Sample_type[[i]]," ")[[1]][2]
				output$Sample_type[[i]] <- "LYS"
			} else {
				output$LYS_num[[i]] <- NA
			}
		}
		
		calc_copies <- function(DNAtype, Cq){
			if(DNAtype=="T3"){
				copies <- 7e10*exp(-0.534*Cq)
				detection_limit <- 7e10*exp(-0.534*40)
			}else if (DNAtype=="T4"){
				copies <- 4e10*exp(-0.621*Cq)
				detection_limit <- 4e10*exp(-0.621*40)
			}else if (DNAtype=="T10"){
				copies <- 1e11*exp(-0.573*Cq)
				detection_limit <- 1e11*exp(-0.573*40)
			}else if (DNAtype=="T11"){
				copies <- 4e11*exp(-0.648*Cq)
				detection_limit <- 4e11*exp(-0.648*40)
			}else if (DNAtype=="T12"){
				copies <- 2e11*exp(-0.598*Cq)
				detection_limit <- 2e11*exp(-0.598*40)
			}
			return(copies)
		}
		output$copies <- NA
		for(i in 1:length(output$Sample)){
			output$copies[[i]] <- calc_copies(DNAtype=DNAtype, Cq=output$Cq.edit[[i]])
		}
		
		if(DNAtype=="T3"){
			detection_limit <- 7e10*exp(-0.534*40)
		}else if (DNAtype=="T4"){
			detection_limit <- 4e10*exp(-0.621*40)
		}else if (DNAtype=="T10"){
			detection_limit <- 1e11*exp(-0.573*40)
		}else if (DNAtype=="T11"){
			detection_limit <- 4e11*exp(-0.648*40)
		}else if (DNAtype=="T12"){
			detection_limit <- 2e11*exp(-0.598*40)
		}
		
		
		for(i in 1:length(output$Sample)){
			if(is.na(output$copies[[i]])){
				output$copies_abv_detection[[i]] <- 0
			}else if(output$copies[[i]]>detection_limit){
				output$copies_abv_detection[[i]] <- output$copies[[i]]
			}else{
				output$copies_abv_detection[[i]] <- 0
			}
		}
		
		outputSNCPC<- output[which(output$Sample_type!="MS"&output$Sample_type!="AS"&output$Sample_type!="LYS"),]
		outputSNPC_filt <- outputSNCPC[,c("Well","Sample","DNAtype","Cq","Sample_type","duplicate","copies","copies_abv_detection")]
		for(i in 1:length(output$Sample)){
			output$copies.Dup.Mean[[i]] <- mean(output$copies[which(output$Sample==output$Sample[[i]])])
			output$copies.Dup.Std.Dev[[i]] <- sd(output$copies[which(output$Sample==output$Sample[[i]])])
		}
		
			repeat{
				n <- readline(prompt="To acknowledge that you are paying attention\n and are about to review the controls for error\n  TYPE: Y")
			if(n=="Y"){
				break
				}
			}
			
			print(outputSNPC_filt)
			
			repeat{
				n <- readline(prompt="If there are no known errors and you wish to continue TYPE: Y\n  if there is an error and you wish to stop TYPE: STOP\n WARNING: IF THERE IS AN ERROR, AND YOU CONTINUE, \nYOU WILL RUIN THE DATABASE AND IT WILL BE ANNOYING TO FIX!")
				if(n=="Y"){
					break
				} else if(n=="STOP"){
					break
				}
			}

		if(n=="STOP"){
			print("function stopped by user")
		} else if(n=="Y"){
			
			output$copies_to_SQ <- output$copies/output$Starting.Quantity..SQ
			output$copies_abv_detection_to_SQ <- output$copies_abv_detection/output$Starting.Quantity..SQ
			
			outputMSAS <- output[which(output$Sample_type=="MS"|output$Sample_type=="AS"),]
			outputMSAS$Sample_date <- as.POSIXct(outputMSAS$Sample_date, format="%m/%d/%Y %H:%M")
			outputLYS <- output[which(output$Sample_type=="LYS"),]
			outputLYS$LYS_num <- factor(outputLYS$LYS_num,levels=c(11,6,1,12,7,2,13,8,3,14,9,4,15,10,5))
			outputLYS$Sample_date <- as.POSIXct(outputLYS$Sample_date, format="%m/%d/%Y %H:%M")
			
			#duplicate plots
			MSASd1copies <- outputMSAS[which(outputMSAS$duplicate=="d1"),]
			MSASd2copies <- outputMSAS[which(outputMSAS$duplicate=="d2"),]
			MSASdup <-merge(MSASd1copies,MSASd2copies,by="Sample")
			
			LYSd1copies <- outputLYS[which(outputLYS$duplicate=="d1"),]
			LYSd2copies <- outputLYS[which(outputLYS$duplicate=="d2"),]
			LYSdup <-merge(LYSd1copies,LYSd2copies,by="Sample")
			
			data_list$MSAS <- rbind.data.frame(data_list$MSAS,outputMSAS)
			data_list$LYS <- rbind.data.frame(data_list$LYS,outputLYS)
			data_list$all <- rbind.data.frame(data_list$all,output)
			data_list$SNCPC <- rbind.data.frame(data_list$SNCPC,outputSNCPC)
			data_list$MSASduplicates <- rbind.data.frame(data_list$MSASduplicates, MSASdup)
			data_list$LYSduplicates <- rbind.data.frame(data_list$LYSduplicates,LYSdup)
			
			for(i in 1:length(data_list)){
				if(length(data_list[[i]][[1]])==0){
					data_list[i] <- list(NULL)
				}
			}
			
			plots["MSAS_copies_date_plot"] <- tryCatch({list(ggplot(data=data_list$MSAS,aes(x=Sample_date, y=copies_abv_detection, color=duplicate, group=duplicate)) +
					geom_line()+geom_point(aes(shape=Sample_type), size=4)+scale_x_datetime(date_labels=c("%m/%d %H:%M"), date_breaks="12 hours",
							date_minor_breaks="2 hours")+
					xlab("\nDatetime") + ylab("Copies\n") + ggtitle("DNA") + facet_wrap(~DNAtype, ncol=1 ))},
				error = function(e) list(NULL))
			plots["MSAS_copies_SQ_date_plot"] <- tryCatch({list(ggplot(data=data_list$MSAS,aes(x=Sample_date, y=copies_abv_detection_to_SQ, color=duplicate, group=duplicate)) +
					geom_line()+geom_point(aes(shape=Sample_type), size=4)+scale_x_datetime(date_labels=c("%m/%d %H:%M"), date_breaks="12 hours",
							date_minor_breaks="2 hours")+
					xlab("\nDatetime") + ylab("Copies/SQ\n") + ggtitle("DNA") + facet_wrap(~DNAtype, ncol=1 ))},
				error = function(e) list(NULL))
			
			# GET EQUATION AND R-SQUARED AS STRING
			# SOURCE: http://goo.gl/K4yh
			
			lm_eqn <- function(df){
				m <- lm(copies_abv_detection.y ~ copies_abv_detection.x, df);
				eq <- substitute(italic(r)^2~"="~r2, 
						list(r2 = format(summary(m)$r.squared, digits = 3)))
				as.character(as.expression(eq));                 
			}
			
			MSAST3 <- data_list$MSASduplicates[which(data_list$MSASduplicates$DNAtype.x=="T3"),]
			MSAST4 <- data_list$MSASduplicates[which(data_list$MSASduplicates$DNAtype.x=="T4"),]
			MSAST10 <- data_list$MSASduplicates[which(data_list$MSASduplicates$DNAtype.x=="T10"),]
			MSAST11 <- data_list$MSASduplicates[which(data_list$MSASduplicates$DNAtype.x=="T11"),]
			MSAST12 <- data_list$MSASduplicates[which(data_list$MSASduplicates$DNAtype.x=="T12"),]
			
			LYST3 <- data_list$LYSduplicates[which(data_list$LYSduplicates$DNAtype.x=="T3"),]
			LYST4 <- data_list$LYSduplicates[which(data_list$LYSduplicates$DNAtype.x=="T4"),]
			LYST10 <- data_list$LYSduplicates[which(data_list$LYSduplicates$DNAtype.x=="T10"),]
			LYST11 <- data_list$LYSduplicates[which(data_list$LYSduplicates$DNAtype.x=="T11"),]
			LYST12 <- data_list$LYSduplicates[which(data_list$LYSduplicates$DNAtype.x=="T12"),]
			
			
			plots["MSASdupdupplotT3"] <-  tryCatch({list(ggplot(data=MSAST3, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
								geom_smooth(method="lm", show.legend=TRUE) +
								geom_text(x = 1500, y = 3500, label = lm_eqn(MSAST3), parse = TRUE)+
								xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype))},
					error = function(e) list(NULL))
			plots["MSASdupdupplotT4"] <-  tryCatch({list(ggplot(data=MSAST4, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
								geom_smooth(method="lm", show.legend=TRUE) +
								geom_text(x = 1500, y = 3500, label = lm_eqn(MSAST4), parse = TRUE)+
								xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype))},
					error = function(e) list(NULL))
			plots["MSASdupdupplotT10"] <-  tryCatch({list(ggplot(data=MSAST10, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
								geom_smooth(method="lm", show.legend=TRUE) +
								geom_text(x = 1500, y = 3500, label = lm_eqn(MSAST10), parse = TRUE)+
								xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype))},
					error = function(e) list(NULL))
			plots["MSASdupdupplotT11"] <-  tryCatch({list(ggplot(data=MSAST11, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
								geom_smooth(method="lm", show.legend=TRUE) +
								geom_text(x = 1500, y = 3500, label = lm_eqn(MSAST11), parse = TRUE)+
								xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype))},
					error = function(e) list(NULL))
			plots["MSASdupdupplotT12"] <-  tryCatch({list(ggplot(data=MSAST12, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
								geom_smooth(method="lm", show.legend=TRUE) +
								geom_text(x = 1500, y = 3500, label = lm_eqn(MSAST12), parse = TRUE)+
								xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype))},
					error = function(e) list(NULL))
			
			
			plots["LYS_copies_date_plot"] <- tryCatch({list(ggplot(data=data_list$LYS,aes(x=Sample_date, y=copies_abv_detection, color=duplicate, group=duplicate)) +
					 geom_line()+geom_point(aes(shape=Sample_type), size=6) +facet_wrap(~LYS_num, ncol=3)+
					 scale_x_datetime(date_labels=c("%m/%d %H:%M")))},
		 		error = function(e) list(NULL))
			 
			 plots["LYS_copies_SQ_date_plot"] <-tryCatch({list(ggplot(data=data_list$LYS,aes(x=Sample_date, y=copies_abv_detection_to_SQ, color=duplicate, group=duplicate)) +
					 geom_line()+geom_point(aes(shape=Sample_type), size=6) +facet_wrap(~LYS_num, ncol=3)+
					 scale_x_datetime(date_labels=c("%m/%d %H:%M")))},
				 error = function(e) list(NULL))
			
		plots["LYSdupdupplotT3"] <-  tryCatch({list(ggplot(data=LYST3, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
							geom_smooth(method="lm", show.legend=TRUE) +
					 		geom_text(x = 1500, y = 3500, label = lm_eqn(LYST3), parse = TRUE)+
					 		xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype)+facet_wrap(~LYS_num, ncol=3))},
		 		error = function(e) list(NULL))
		plots["LYSdupdupplotT4"]  <-  tryCatch({list(ggplot(data=LYST4, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
							geom_smooth(method="lm", show.legend=TRUE) +
							geom_text(x = 1500, y = 3500, label = lm_eqn(LYST4), parse = TRUE)+
							xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype)+facet_wrap(~LYS_num, ncol=3))},
				error = function(e) list(NULL))
		plots["LYSdupdupplotT10"]  <-  tryCatch({list(ggplot(data=LYST10, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
							geom_smooth(method="lm", show.legend=TRUE) +
							geom_text(x = 1500, y = 3500, label = lm_eqn(LYST10), parse = TRUE)+
							xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype)+facet_wrap(~LYS_num, ncol=3))},
				error = function(e) list(NULL))
		plots["LYSdupdupplotT11"]  <-  tryCatch({list(ggplot(data=LYST11, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
							geom_smooth(method="lm", show.legend=TRUE) +
							geom_text(x = 1500, y = 3500, label = lm_eqn(LYST11), parse = TRUE)+
							xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype)+facet_wrap(~LYS_num, ncol=3))},
				error = function(e) list(NULL))
		plots["LYSdupdupplotT12"]  <-  tryCatch({list(ggplot(data=LYST12, aes(x=copies_abv_detection.x, y=copies_abv_detection.y)) +geom_point()+
							geom_smooth(method="lm", show.legend=TRUE) +
							geom_text(x = 1500, y = 3500, label = lm_eqn(LYST12), parse = TRUE)+
							xlab("\nCopies (Duplicate 2)")+ylab("Copies (Duplicate 1)\n")+ ggtitle(DNAtype)+facet_wrap(~LYS_num, ncol=3))},
				error = function(e) list(NULL))
			
			 currentdatetime <- gsub(":","-", as.character(Sys.time()))
			 
			 	
		save("data_list","plots",file=paste(master_dataset_output_folderpath,"/MasterData_created_",currentdatetime,".RData", sep=""))
		
			for(i in 1:length(plots)){
				if(length(plots[[i]][[1]])>1){
					plotname <- names(plots)[[i]]
					ggsave(filename=paste(plot_output_folderpath,"/",plotname,".png",sep=""),plots[[i]],height=height,width=width)
				}else{
					plotname <-names(plots)[[i]]
					print(paste(plotname," plot not created because of empty dataset", sep=""))
				}
			}
			print("completed")
		 }
	 }
}

 #########
 data_list <- vector("list",6)
 names(data_list) <- c("MSAS","LYS","SNCPC","all","MSASduplicates","LYSduplicates")
 data_list$MSAS <- outputMSAS
 data_list$LYS <- outputLYS
 data_list$all <- output
 data_list$SNCPC <- output$SNCPC
 data_list$MSASduplicates <- MSASdup
 data_list$LYSduplicates <- LYSdup
 
 
 plots <- vector("list", 6)
 names(plots) <- c("MSAS_copies_date_plot","MSAS_copies_SQ_date_plot","MSASdupdupplot", "LYS_copies_date_plot",
		 "LYS_copies_SQ_date_plot","LYSdupdupplot")
 plots$MSAS_copies_date_plot <- MSAS_copies_date_plot
 plots$MSAS_copies_SQ_date_plot <- MSAS_copies_SQ_date_plot
 plots$MSASdupdupplotT3 <- MSASdupdupplotT3
 plots$MSASdupdupplotT4 <- MSASdupdupplotT4
 plots$MSASdupdupplotT10 <- MSASdupdupplotT10
 plots$MSASdupdupplotT11 <- MSASdupdupplotT11
 plots$MSASdupdupplotT12 <- MSASdupdupplotT12
 plots$LYS_copies_date_plot <- LYS_copies_date_plot
 plots$LYS_copies_SQ_date_plot <- LYS_copies_SQ_date_plot
 plots$LYSdupdupplotT3 <- LYSdupdupplotT3
 plots$LYSdupdupplotT4 <- LYSdupdupplotT4
 plots$LYSdupdupplotT10 <- LYSdupdupplotT10
 plots$LYSdupdupplotT11 <- LYSdupdupplotT11
 plots$LYSdupdupplotT12 <- LYSdupdupplotT12
 
 
 
 data_list <- vector("list",6)
 names(data_list) <- c("MSAS","LYS","SNCPC","all","MSASduplicates","LYSduplicates")
 data_list["MSAS"] <- list(NULL)
 data_list["LYS"] <- list(NULL)
 data_list["all"] <- list(NULL)
 data_list["SNCPC"] <- list(NULL)
 data_list["MSASduplicates"] <- list(NULL)
 data_list["LYSduplicates"] <- list(NULL)
 
 
 plots <- vector("list", 4)
 names(plots) <- c("MSAS_copies_date_plot","MSAS_copies_SQ_date_plot", "LYS_copies_date_plot",
		 "LYS_copies_SQ_date_plot")
 plots["MSAS_copies_date_plot"] <- list(NULL)
 plots["MSAS_copies_SQ_date_plot"] <- list(NULL)
 plots["MSASdupdupplotT3"] <-list(NULL)
 plots["MSASdupdupplotT4"] <- list(NULL)
 plots["MSASdupdupplotT10"] <- list(NULL)
 plots["MSASdupdupplotT11"] <- list(NULL)
 plots["MSASdupdupplotT12"] <- list(NULL)
 plots["LYS_copies_date_plot"] <- list(NULL)
 plots["LYS_copies_SQ_date_plot"] <- list(NULL)
 plots["LYSdupdupplotT3"] <- list(NULL)
 plots["LYSdupdupplotT4"] <- list(NULL)
 plots["LYSdupdupplotT10"] <- list(NULL)
 plots["LYSdupdupplotT11"] <- list(NULL)
 plots["LYSdupdupplotT12"] <- list(NULL)