setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/greenhouse/data/census_files")

getwd()

library(writexl)


###### season 1 ######


##Data
season1 <-read.csv("2gen_season_1_prevern_census.csv", header=T)
sapply(season1,class)

head(season1,20)
tail(season1)

season1$Master_notes<-as.character(season1$Master_notes)
season1$Flower_number <as.integer(season1$Flower_number)

reproductive<-subset(season1, state=="reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("2gen_season_1_prevern_census_empty.csv", header=T)
Empty <- Empty_input[order(Empty_input$exp_ID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(season1 $exp_ID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(season1 $exp_ID)[i]
  sub<-subset(season1, exp_ID==target_ID)
  
  #        ##If the subset has less than 2 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==2)==F)
  { print(sub$exp_ID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(season1$exp_ID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(season1$exp_ID)[i]
  sub<-subset(season1, exp_ID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$exp_ID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0, not in this experiment, is a greenhouse study
  ifelse((sub$status[1]=="alive"), Empty$Oververn_survival[target_row]<-1, Empty$Oververn_survival[target_row]<-0)
  
  
  #season_survival - Census 2 was the last full census, alive=1 dead=0
  ifelse((sub$status[2]=="alive"), Empty$Season_survival[target_row]<-1, Empty$Season_survival[target_row]<-0)
  
  #season_Include - Census 2 was the last full census, alive=1 dead=0
  #ifelse((sub$Include[2]=="Include"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude[target_row]<-sub$Include[2] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
 # if(Empty$Season_survival_2022[target_row]==1) {
 #   Empty$Overwinter_survival_2022[target_row]<-1
 # }
  
  
  # leaf damage not in file
  
  #Empty$num_total3_2022[target_row]<-sub$Leaf_number[3] 
  #Empty$num_damaged3_2022[target_row]<-sub$Damaged_leaves[3] 
  #Empty$total_damage3_2022[target_row]<-sub$total_damage[3]
  #Empty$avg_damage3_2022[target_row]<-sub$average_damage[3]
  #Empty$LAR_3_2022[target_row]<-sub$LAR[3] 
  
  #Empty$num_total7_2022[target_row]<-sub$Leaf_number[7] 
  #Empty$num_damaged7_2022[target_row]<-sub$Damaged_leaves[7] 
  #Empty$total_damage7_2022[target_row]<-sub$total_damage[7]
  #Empty$avg_damage7_2022[target_row]<-sub$average_damage[7]
  #Empty$LAR_7_2022[target_row]<-sub$LAR[7] 
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes[2] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival[target_row]==0) { #Outer if statement
    
    #Subset data where status==Dead	
    First_death<-subset(sub, status=="dead")
    
    #Extract row name of last census where status==alive
    Last_alive<-rownames(tail(sub[which(sub$status=="alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
 
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$state=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted[target_row]<-1, Empty$Bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, state=="reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak[target_row]<-max_fl
    
    Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
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
    
    
    #Height3_silique
    Empty$Height1_silique[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques[target_row]<-sum(sub[,c(45,
                                                       47,
                                                       49,
                                                       51,
                                                       53,
                                                       55,
                                                       57,
                                                       59,
                                                       61,
                                                       63,
                                                       65,
                                                       67,
                                                       69,
                                                       71,
                                                       73,
                                                       75,
                                                       77,
                                                       79,
                                                       81,
                                                       83,
                                                       85,
                                                       87)], na.rm=T)
  
} #closes inner if statement



#check1 <-subset(Empty, PlantID=="G16375")
#check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2gen_season1_summary.xlsx", sheetName="SummarySeason1", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 


###### season 2 ######


##Data
season2 <-read.csv("2gen_season_2_census.csv", header=T)
sapply(season2,class)

head(season2,20)
tail(season2)

season2$Master_notes<-as.character(season2$Master_notes)
season2$Flower_number <as.integer(season2$Flower_number)

reproductive<-subset(season2, state=="reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("2gen_season_2_census_empty.csv", header=T)
Empty <- Empty_input[order(Empty_input$exp_ID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(season2 $exp_ID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(season2 $exp_ID)[i]
  sub<-subset(season2, exp_ID==target_ID)
  
  #        ##If the subset has less than 12 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==12)==F)
  { print(sub$exp_ID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(season2$exp_ID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(season2$exp_ID)[i]
  sub<-subset(season2, exp_ID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$exp_ID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0, not in this experiment, is a greenhouse study
  ifelse((sub$status[1]=="alive"), Empty$Oververn_survival[target_row]<-1, Empty$Oververn_survival[target_row]<-0)
  
  
  #season_survival - Census 2 was the last full census, alive=1 dead=0
  ifelse((sub$status[2]=="alive"), Empty$Season_survival[target_row]<-1, Empty$Season_survival[target_row]<-0)
  
  #season_Include - Census 2 was the last full census, alive=1 dead=0
  #ifelse((sub$Include[2]=="Include"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude[target_row]<-sub$Include[2] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
  # if(Empty$Season_survival_2022[target_row]==1) {
  #   Empty$Overwinter_survival_2022[target_row]<-1
  # }
  
  
  # leaf damage not in file
  
  #Empty$num_total3_2022[target_row]<-sub$Leaf_number[3] 
  #Empty$num_damaged3_2022[target_row]<-sub$Damaged_leaves[3] 
  #Empty$total_damage3_2022[target_row]<-sub$total_damage[3]
  #Empty$avg_damage3_2022[target_row]<-sub$average_damage[3]
  #Empty$LAR_3_2022[target_row]<-sub$LAR[3] 
  
  #Empty$num_total7_2022[target_row]<-sub$Leaf_number[7] 
  #Empty$num_damaged7_2022[target_row]<-sub$Damaged_leaves[7] 
  #Empty$total_damage7_2022[target_row]<-sub$total_damage[7]
  #Empty$avg_damage7_2022[target_row]<-sub$average_damage[7]
  #Empty$LAR_7_2022[target_row]<-sub$LAR[7] 
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes[2] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival[target_row]==0) { #Outer if statement
    
    #Subset data where status==Dead	
    First_death<-subset(sub, status=="dead")
    
    #Extract row name of last census where status==alive
    Last_alive<-rownames(tail(sub[which(sub$status=="alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
  
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$state=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted[target_row]<-1, Empty$Bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, state=="reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak[target_row]<-max_fl
    
    Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
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
    
    
    #Height3_silique
    Empty$Height1_silique[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques[target_row]<-sum(sub[,c(45,
                                                       47,
                                                       49,
                                                       51,
                                                       53,
                                                       55,
                                                       57,
                                                       59,
                                                       61,
                                                       63,
                                                       65,
                                                       67,
                                                       69,
                                                       71,
                                                       73,
                                                       75,
                                                       77,
                                                       79,
                                                       81,
                                                       83,
                                                       85,
                                                       87)], na.rm=T)
  
} #closes inner if statement



#check1 <-subset(Empty, PlantID=="G16375")
#check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2gen_season2_summary.xlsx", sheetName="SummarySeason2", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 

###### season 3 ######


##Data
season3 <-read.csv("2gen_season_3_census.csv", header=T)
sapply(season3,class)

head(season3,20)
tail(season3)

season3$Master_notes<-as.character(season3$Master_notes)
season3$Flower_number <as.integer(season3$Flower_number)

reproductive<-subset(season3, state=="reproductive")
head(reproductive,20)
tail(reproductive,20)




##Empty data set that we will populate with summary data


Empty_input<-read.csv("2gen_season_3_census_empty.csv", header=T)
Empty <- Empty_input[order(Empty_input$exp_ID),]
row.names(Empty) <- 1:nrow(Empty)


sapply(Empty,class)
head(Empty,20)
Empty $Master_notes<-as.character(Empty $Master_notes)




##Do any plants have an incomplete number of censuses listed in the dataset?

for (i in 1:length(unique(season3 $exp_ID))) {
  
  #        ##Subset for each plant
  target_ID<-unique(season3 $exp_ID)[i]
  sub<-subset(season3, exp_ID==target_ID)
  
  #        ##If the subset has less than 12 rows (number of censuses), tell me the plant ID number
  if((nrow(sub)==4)==F)
  { print(sub$exp_ID[1]) }
  
}



##The loop begins - 
for (i in 1:length(unique(season3$exp_ID))){ 
  
  ##Subset the data for plantID i
  target_ID<-unique(season3$exp_ID)[i]
  sub<-subset(season3, exp_ID==target_ID)
  
  
  ##Identify the row number in Empty that corresponds to this ID number 
  target_row<-as.numeric(rownames(Empty[which(Empty$exp_ID==target_ID),]))
  
  #overwinter_survival - Census 1 was the first full census, alive=1 dead=0, not in this experiment, is a greenhouse study
  ifelse((sub$status[1]=="alive"), Empty$Oververn_survival[target_row]<-1, Empty$Oververn_survival[target_row]<-0)
  
  
  #season_survival - Census 2 was the last full census, alive=1 dead=0
  ifelse((sub$status[2]=="alive"), Empty$Season_survival[target_row]<-1, Empty$Season_survival[target_row]<-0)
  
  #season_Include - Census 2 was the last full census, alive=1 dead=0
  #ifelse((sub$Include[2]=="Include"), Empty$Season_survival_2022[target_row]<-1, Empty$Season_survival_2022[target_row]<-0)
  Empty$Exclude[target_row]<-sub$Include[2] 
  
  #error here
  ##In some cases, plants initially marked as dead are later found alive - we need to adjust Overwinter_survival to equal 1 if plants were alive in the final census
  
  # if(Empty$Season_survival_2022[target_row]==1) {
  #   Empty$Overwinter_survival_2022[target_row]<-1
  # }
  
  
  # leaf damage not in file
  
  #Empty$num_total3_2022[target_row]<-sub$Leaf_number[3] 
  #Empty$num_damaged3_2022[target_row]<-sub$Damaged_leaves[3] 
  #Empty$total_damage3_2022[target_row]<-sub$total_damage[3]
  #Empty$avg_damage3_2022[target_row]<-sub$average_damage[3]
  #Empty$LAR_3_2022[target_row]<-sub$LAR[3] 
  
  #Empty$num_total7_2022[target_row]<-sub$Leaf_number[7] 
  #Empty$num_damaged7_2022[target_row]<-sub$Damaged_leaves[7] 
  #Empty$total_damage7_2022[target_row]<-sub$total_damage[7]
  #Empty$avg_damage7_2022[target_row]<-sub$average_damage[7]
  #Empty$LAR_7_2022[target_row]<-sub$LAR[7] 
  
  
  ## Master notes
  Empty$Master_notes[target_row]<-sub$Master_notes[2] 
  
  #Date_first_death_observed - only assigns a date of first death for individuals who are not alive at the end of the season - This code does account for the individuals who were listed as dead, then were found as alive, and later died	
  
  if(Empty$Season_survival[target_row]==0) { #Outer if statement
    
    #Subset data where status==Dead	
    First_death<-subset(sub, status=="dead")
    
    #Extract row name of last census where status==alive
    Last_alive<-rownames(tail(sub[which(sub$status=="alive"),], n=1))
    
    ##Adjust First_death data set to include only data 	after the last alive status was observed
    if(length(Last_alive)>=1) { #Inner if statement
      
      ##Adjust subset of First_death
      First_death<-subset(First_death, as.numeric(rownames(First_death))>as.numeric(Last_alive))
      
    } #Closes inner if statement
    
    ##Assign the date of first death
    Empty$Date_first_death_observed[target_row]<-First_death$Ordinal_date[1] 
    
  } #Closes outer if statement
  
  
  #Bolted and Flowered, yes=1 no=0
  #to exclude plants that were erroneously marked as bolting, Jill added a condition that a plant needed to be censused as reproductive at least two times
  ifelse((sum(sub$state=="reproductive", na.rm=T)>=2 | sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Bolted[target_row]<-1, Empty$Bolted[target_row]<-0)
  
  #First, assign a 0 if it never had flowers or siliques, otherwise assign a 1
  ifelse((sum(sub$Flower_number, na.rm=T)>=1 | sum(sub$Silique_number, na.rm=T)>=1), Empty$Flowered[target_row]<-1, Empty$Flowered[target_row]<-0)
  
  ####Subset reproductive data
  Repro<-subset(sub, state=="reproductive" & (Flower_number>=1 | Silique_number>=1)) #this filters out the rows where the plant was reproductive, but is missing flower or silique number
  
  
  #Date_flowering (really this is the date that a plant is first observed to be reproductive)
  if(nrow(Repro)>=1) {
    Empty$Date_flowering[target_row]<-Repro$Ordinal_date[1] }
  
  
  
  #Stem_number_flowering, Flower_number_flowering, Silique_number_flowering, and Silique_length_flowering
  if(nrow(Repro)>=1) {
    Empty$Stem_number_flowering[target_row]<-Repro$Stem_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Flower_number_flowering[target_row]<-Repro$Flower_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_number_flowering[target_row]<-Repro$Silique_number[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Silique_length_flowering[target_row]<-Repro$Silique_length[1] 
    #<-Repro$Silique_length[1] 
  }
  
  
  
  
  #Height1 through Height3 at first flower
  if(nrow(Repro)>=1) {
    Empty$Height1_flowering[target_row]<-Repro$Height1[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height2_flowering[target_row]<-Repro$Height2[1] }
  
  if(nrow(Repro)>=1) {
    Empty$Height3_flowering[target_row]<-Repro$Height3[1] }
  
  ########Note - this if statement contains summary commands for all peak-flowering related data
  
  if(sum(Repro$Flower_number, na.rm=T)>=1) { #closes below
    
    #Date_peak_flowering		
    max_fl<-max(Repro$Flower_number, na.rm=T)
    
    Empty$Date_peak_flowering[target_row]<-Repro$Ordinal_date[which(Repro$Flower_number==max_fl)][1]
    
    
    #Flower_number_peak and Silique_number_peak and Stem_number_peak
    Empty$Flower_number_peak[target_row]<-max_fl
    
    Empty$Silique_number_peak[target_row]<-Repro$Silique_number[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Stem_number_peak[target_row]<-Repro$Stem_number[which(Repro$Flower_number==max_fl)][1]
    
    
    #Height1_peak through Height3_peak
    Empty$Height1_peak[target_row]<-Repro$Height1[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height2_peak[target_row]<-Repro$Height2[which(Repro$Flower_number==max_fl)][1]
    
    Empty$Height3_peak[target_row]<-Repro$Height3[which(Repro$Flower_number==max_fl)][1]
    
    
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
    
    
    #Height3_silique
    Empty$Height1_silique[target_row]<-Repro$Height1[1]
    
    Empty$Height2_silique[target_row]<-Repro$Height2[1]
    
    Empty$Height3_silique[target_row]<-Repro$Height3[1]
    
    
    
  } #Closing bracket from outer if statement
  
  
  
  
  
  
  
  #Mature_silique_number
  Empty$Mature_silique_number[target_row]<-sum(sub$Number_collected_siliques, na.rm=T)
  
  
  #Failed_silique_number 
  Empty$Failed_silique_number[target_row]<-sum(sub$Failed_silique_number, na.rm=T)
  
  ##Failed_silique_number - if statement needed because you can't take the maximum of a vector of NAs (when vegetative)
  #if (sum(sub$Silique_number, na.rm=T)>=1) { 
  
  # Empty$Failed_silique_number[target_row]<-max(sub$Silique_number, na.rm=T)-sum(sub$Number_collected_siliques, na.rm=T)
  
  ##If the number of failed siliques is negative, set to zero
  #if(Empty$Failed_silique_number[target_row]<0) {
  
  # Empty$Failed_silique_number[target_row]<-0
  
  #Mature_length_siliques - note that the column numbers have to be updated
  Empty$Mature_length_siliques[target_row]<-sum(sub[,c(45,
                                                       47,
                                                       49,
                                                       51,
                                                       53,
                                                       55,
                                                       57,
                                                       59,
                                                       61,
                                                       63,
                                                       65,
                                                       67,
                                                       69,
                                                       71,
                                                       73,
                                                       75,
                                                       77,
                                                       79,
                                                       81,
                                                       83,
                                                       85,
                                                       87)], na.rm=T)
  
} #closes inner if statement



#check1 <-subset(Empty, PlantID=="G16375")
#check2 <-subset(Empty, PlantID=="G16316")

library(openxlsx)

#Empty <- Empty[order(Empty $Datafile_position),]
#row.names(Empty) <- 1:nrow(Empty)

write.xlsx(Empty, "2gen_season3_summary.xlsx", sheetName="SummarySeason3", colnames=TRUE, rownames=TRUE, keepNA=TRUE) 






