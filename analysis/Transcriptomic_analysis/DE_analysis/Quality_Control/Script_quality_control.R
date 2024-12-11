## RNASEQ Quality Control
## Coded by Marion Verdenaud and Gabrièle Adam
## Run using functions developped in https://forgemia.inra.fr/GNet/dicoexpress/

## Set working directory
Snakemake_Work_Directory<-getwd()

setwd("Template_scripts")
source("Sources/Load_Functions.R")

Load_Functions()

## Inputs
Data_Directory<-paste0(Snakemake_Work_Directory,"/input")
Results_Directory<-paste0(Snakemake_Work_Directory,"/Quality_Control")

## Creation of output repository
dir.create(Results_Directory)
# print (Data_Directory)

## Load Data Files
Project_Name <- "Cts"

##Filter Dataset by condition
Filter=NULL
Filter_condition=''
# print (paste("filterdondition",Filter_condition))

Sep=","
# print (paste("Filter",Filter))

Data_Files <- Load_Data_Files(Data_Directory, Project_Name,Filter,Sep)

Project_Name <- Data_Files$Project_Name
Target <- Data_Files$Target
Raw_Counts <- Data_Files$Raw_Counts
Annotation <- Data_Files$Annotation
Reference_Enrichment<-Data_Files$Reference_Enrichment

## Quality Control
Filter_Strategy="NbReplicates"
CPM_Cutoff=5
Normalization_Method="TMM"


## Model
Replicate=TRUE
Interaction=TRUE

Model <- GLM_Contrasts(Results_Directory, Project_Name,
Target, Replicate, Interaction)

GLM_Model <- Model$GLM_Model
Contrasts <- Model$Contrasts

#Save contrast matrix to use with DEseq2
saveRDS(Contrasts,paste0(Results_Directory,"/",Project_Name,"/DiffAnalysis/","Contrasts.rds"))

Quality_Control(Data_Directory, Results_Directory, Project_Name,
                Target, Raw_Counts, Filter_Strategy,
                Color_Group=NULL,CPM_Cutoff,
                Normalization_Method)

sessioninfopath<-paste0(Results_Directory,"/",Project_Name,"/","SessionInfo.txt")
