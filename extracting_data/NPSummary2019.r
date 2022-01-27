setwd("~/Documents/Jill/postdoc/RMBL/2019/final_rosette")
getwd()

library(writexl)

##Data
NP <-read.delim("NPCensus2019.txt", header=T)
sapply(NP,class)


#NP_input<-read.delim("NPCensus2019.txt", header=T)
#sapply(NP_input,class)
##Sort by plantid
#NP <- NP_input[order(NP_input$PlantID),]
#row.names(NP) <- 1:nrow(NP)

head(NP,20)
tail(NP)

NP$Master_notes<-as.character(NP$Master_notes)

Reproductive<-subset(NP, State=="Reproductive")
head(Reproductive,20)


##Empty data set that we will populate with summary data
Empty <-read.csv("NPEmpty.csv", header=T)

#Empty_input<-read.csv("NPEmpty.csv", header=T)
#Empty <- Empty_input[order(Empty_input$PlantID),]
#row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes19<-as.character(Empty $Master_notes19)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(NP $PlantID))) {

        ##Subset for each plant
target_ID<-unique(NP $PlantID)[i]
sub<-subset(NP, PlantID==target_ID)

        ##If the subset has less than 8 rows, tell me the plant ID number
if((nrow(sub)==8)==F) 
{ print(sub$PlantID[1]) }

}



##The loop begins - 
for (i in 1:length(unique(NP$PlantID))) {
	
##Subset the data for plantID i
	target_ID<-unique(NP$PlantID)[i]
	sub<-subset(NP, PlantID==target_ID)
	
	
##Identify the row number in Empty that corresponds to this ID number 
	target_row<-as.numeric(rownames(Empty[which(Empty$PlantID==target_ID),]))
	

#Overwinter_survival - Census 5 is the leaf damage census, alive=1 dead=0
	ifelse((sub$Status[1]=="Alive" | sub$Status[2]=="Alive" | sub$Status[3]=="Alive" |sub$Status[4]=="Alive" | sub$Status[5]=="Alive"), Empty$Overwinter_survival[target_row]<-1, Empty$Overwinter_survival[target_row]<-0)

	
#Season_survival - Census 8 was the last full census, alive=1 dead=0
	ifelse((sub$Status[8]=="Alive"), Empty$Season_survival[target_row]<-1, Empty$Season_survival[target_row]<-0)
		
	
##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census or if they had leaves counted during the damage census (note that we can only ask if there are >= 1 leaf if the element is not NA)
	if(Empty$Season_survival[target_row]==1 | (Empty$Leaf_Number[target_row]>=1 & is.na(Empty$Leaf_Number[target_row])==F)) {
		Empty$Overwinter_survival[target_row]<-1
	}

##I measured leaf damage at census 5	

Empty$Leaf_Number[target_row]<-sub$Leaf_Number[5] 
Empty$Damaged_Leaf_Number[target_row]<-sub$Damaged_Leaf_Number[5] 
Empty$Percent_Leaf_Damage[target_row]<-sub$Percent_Leaf_Damage[5] 
Empty$LAR[target_row]<-sub$LAR[5] 

## Percent vegetation at census 8
Empty$Percent_Veg_Cover[target_row]<-sub$Percent_Veg_Cover[8] 

## Master notes
Empty$Master_notes19[target_row]<-sub$Master_notes[8] 
	
#Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	

	if(Empty$Season_survival[target_row]==0) { #Outer if statement
	
	#Subset data where Status==Dead	
	First_death<-subset(sub, Status=="Dead")
	
	#Extract row name of last census where status==Alive
	Last_alive<-rownames(tail(sub[which(sub$Status=="Alive"),], n=1))
	
	##Adjust First_death data set to include only data 	after the last Alive status was observed
	if(length(Last_alive)>=1) { #Inner if statement
	
	##Adjust subset of First_death
	First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
		
	} #Closes inner if statement
	
	##Assign the date of first death
	Empty$Date_first_death_observed[target_row]<-First_death$Ordinal_date[1] 
	
} #Closes outer if statement

	
	
##Gopher damage or Censor - was coded as GopherDamage or Censor in the Exclude column - sometimes a plant marked as being exposed to gopher damage doesn't die - for that reason, this has to meet conditions that it was censused as gopher damaged at least once and that we didn't observe it as alive in the last census
	if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive")) {
	Empty$Exclude[target_row]<-"GopherDamage"
	}

	if(any(sub$Exclude=="Censor") & (sub$Status[8]!="Alive")) {
	Empty$Exclude[target_row]<-"Censor"
	}
	
	if(is.na(Empty$Exclude[target_row])==T) {
	Empty$Exclude[target_row]<-"Include"
	}

	
##Date_first_gopher, again only noted if the plant was not marked as alive in the last census
##This code is a little crazy because we are also calculating the date a plant was last alive prior to gopher damage and whether the living plant died prior to gopher damage
if(any(sub$Exclude=="GopherDamage") & (sub$Status[8]!="Alive" )) {
		
	Empty$Date_first_gopher[target_row]<-sub$Ordinal_date[which(sub$Exclude=="GopherDamage")][1] 
	
	##Last date of being alive prior to gopher damage
	Last_alive<-tail(sub$Ordinal_date[which(sub$Status=="Alive" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)
		
	#Last date of being dead prior to gopher damage
	Last_dead<-tail(sub$Ordinal_date[which(sub$Status=="Dead" & sub$Ordinal_date<Empty$Date_first_gopher[target_row])], n=1)	
	
	#If the plant has both a Last_alive and Last_dead date
	if(length(Last_alive)>=1 & length(Last_dead)>=1) {
	
	#If the plant was alive before it was dead and if it died before gopher was first observed - have to separate this if statement from the last because you can't ask if something is less than something else if never occured
	if(Last_alive<Last_dead & Last_dead<Empty$Date_first_gopher[target_row]) {
		
		#Include plant
		Empty$Exclude[target_row]<-"Include"
			
	} } #closes inner if statements
	
	##Lastly, if plant died prior to gopher damage, list it as a natural death
	if(Empty$Date_first_death_observed[target_row]<Empty$Date_first_gopher[target_row]) {
		
		Empty$Exclude[target_row]<-"Include"
				
	} #Closes inner if statement
	
} #Closes outer if statement

	
	

	
#Bolted and Flowered, yes=1 no=0
	#to exclude plants that were erroneously marked as bolting, I added a codition that a plant needed to be censused as reproductive at least two times
	ifelse((sum(sub$State=="Reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted[target_row]<-1, Empty$Bolted[target_row]<-0)
	
	#First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
	ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)




####Subset reproductive data
Repro<-subset(sub, State=="Reproductive" & (Flower_number>=1 | Silique_number>=1))


#Date_flowering (really this is the date that a plant is first observed to be reproductive)
	if(nrow(Repro)>=1) {
	Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
	

#Flower_color - note the nested ifelse statements - this is needed to assign a 1 for plants with pigmented flowers, a 0 to plants with white flowers, and an NA to plants that never flowered (or for which we never observed flowers) 
	if(any(Repro$Flower_color=="Pigmented")) {
	Empty$Flower_color[target_row]<-1 }
	
	if(all(Repro$Flower_color!="Pigmented") & sum(Repro$Flower_number, na.rm=T)>=1) {
	Empty$Flower_color[target_row]<-0 }


#Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
	if(nrow(Repro)>=1) {
	Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }

	if(nrow(Repro)>=1) {
	Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }

	if(nrow(Repro)>=1) {
	Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
	
	if(nrow(Repro)>=1) {
	Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] }
	

#Height1 through Height3 at first flower
	if(nrow(Repro)>=1) {
	Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
			
	if(nrow(Repro)>=1) {
	Empty$Height2_flowering[target_row]<-Repro$Height2[1] }

	if(nrow(Repro)>=1) {
	Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
	
		if(nrow(Repro)>=1) {
	Empty$Height4_flowering[target_row]<-Repro$Height4[1] }
			
	if(nrow(Repro)>=1) {
	Empty$Height5_flowering[target_row]<-Repro$Height5[1] }

	if(nrow(Repro)>=1) {
	Empty$Height6_flowering[target_row]<-Repro$Height6[1] }


########Note - this if statement contains summary commands for all peak-flowering related data

if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below

#Date_peak_flowering		
	max_fl<-max(Repro$Flower_number, na.rm=T)
	
	Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
			
			
#Flower_number_peak and Silique_number_peak and Stem_number_peak
	Empty$Flower_number_peak[target_row]<-max_fl
	
	Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
	
	Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]


#Height1_peak through Height5_peak
	Empty$Height1_peak[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
		
	Empty$Height2_peak[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]

	Empty$Height3_peak[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]

	Empty$Height4_peak[target_row]<-Repro$Height4[which(Repro$Flower_number==max_fl)][1]
		
	Empty$Height5_peak[target_row]<-Repro$Height5[which(Repro$Flower_number==max_fl)][1]

	Empty$Height6_peak[target_row]<-Repro$Height6[which(Repro$Flower_number==max_fl)][1]
	

		
} #Closes if statement




########Data for date of first SO - Note that this is the date of first SO observed after the last flower was observed

if(sum(Repro$Silique_number, na.rm=T)>=1) { #outer if statement

	#Extract the rowname for the last day a flower was observed, and then the rownumber of that rowname within sub (note that the rownames in sub will come from its name within the master, long-format data file)
	last_Fl<-rownames(tail(Repro[which(Repro$Flower_number>=1),], n=1))
	
	##If there is a date of last flower, adjust the repro subset to reflect only data after the date of last flower
	if(length(last_Fl)>0) { #inner if statement

	##Row number of date of last flower in sub
	last_Fl<-which(grepl(last_Fl, rownames(Repro)))

	##Subset data to exclude censuses on or before the last flower was observed
	
	Repro<-Repro[-c(1:last_Fl),]
		
} #Closes inner if statement


#Subset for censuses where a silique was observed (we have to do this second subset to exlude plants that became b or b, tco after flowering)
	Repro<-subset(Repro, Silique_number>=1)

	
##Now estimate the date of first silique
	Empty$Date_silique[target_row]<-Repro$Ordinal_date[1]


#Stem_number_silique, Silique_number_silique, and Silique_length_silique
	Empty$Stem_number_silique[target_row]<-Repro$Stem_number[1] 
	
	Empty$Silique_number_silique[target_row]<-Repro$Silique_number[1]
	
	Empty$Silique_length_silique[target_row]<-Repro$Silique_length[1]


#Height1_silique
	Empty$Height1_silique[target_row]<-Repro$Height1[1]
	
	Empty$Height2_silique[target_row]<-Repro$Height2[1]

	Empty$Height3_silique[target_row]<-Repro$Height3[1]

	Empty$Height4_silique[target_row]<-Repro$Height4[1]
	
	Empty$Height5_silique[target_row]<-Repro$Height5[1]

	Empty$Height6_silique[target_row]<-Repro$Height6[1]

	
} #Closing bracket from outer if statement





#Mature_silique_number
	Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)


#Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
	if (sum(sub$Silique_number, na.rm=T)>=1) { 
		
	Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)

	##If the number of failed siliques is negative, set to zero
	if(Empty$Failed_silique_number[target_row]<0) {
		
		Empty$Failed_silique_number[target_row]<-0
		
	} #closes inner if statement
} #closes outer if statement


#Mature_length_siliques - note that the column numbers have to be updated
	Empty$Mature_length_siliques[target_row]<-sum(sub[,c(34,36,38,40,42,44,46,48,50,52,54,56,58,60,62)], na.rm=T)


#Total_flower_number
	Empty$Total_flower_number[target_row]<-sum(sub$Flower_number, na.rm=T)
	
#Mature_entire_silique_length (when %PD is measured) - first vector has column numbers of silique lengths, second vector has column numbers of PD% - columns must be updated before calculating
	Empty$Mature_PD[target_row]<-(sum(sub[,c(35,37,39,41,43,45,47,49,51,53,55,57,59,61,63)], na.rm=T))/sum(sub$Number_collected_siliques, na.rm=T)
	
	
} ##Closes entire loop


N3684<-subset(Empty, PlantID=="N3684")
N8385<-subset(Empty, PlantID=="N8385")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "NPSummary2019.xlsx", sheetName="Summary2019", col.names=TRUE, row.names=TRUE, keepNA=TRUE) 
