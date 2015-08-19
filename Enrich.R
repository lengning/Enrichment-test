
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # score file, rows are genes
Local=options[2] # local list; rows are lists
Lowersetsize=options[3]
if(length(options)<2)Local=NULL
if(length(options)<3)Lowersetsize=10



# csv or txt
tmp=strsplit(File, split="\\.")[[1]]
FileType=tmp[length(tmp)]

if(FileType=="csv"){
	cat("\n Read in csv file \n")
	prefix=strsplit(File,split="\\.csv")[[1]][1]
	In=read.csv(File,stringsAsFactors=F,row.names=1)
}
if(FileType!="csv"){
	cat("\n Read in tab delimited file \n")
	prefix=strsplit(File,split=paste0("\\.",FileType))[[1]][1]
	In=read.table(File,stringsAsFactors=F,row.names=1)
	if(!is.numeric(In[[1]]))In=read.table(File,stringsAsFactors=F,row.names=1,header=T)
}



if(!is.null(Local)){
	tmp2=strsplit(Local, split="\\.")[[1]]
	FileType2=tmp2[length(tmp2)]

	if(FileType2=="csv"){
		cat("\n Read in csv file (gene list)\n")
		ListIn=read.csv(Local,stringsAsFactors=F,row.names=1)
	}
	if(FileType2!="csv"){
		cat("\n Read in tab delimited file (gene list)\n")
		ListIn=read.table(Local,stringsAsFactors=F,row.names=1)
	}

	List=sapply(1:nrow(ListIn),function(i)setdiff(as.vector(ListIn[i,]),c(""," ")))
	names(List)=rownames(ListIn)

} else List=NULL

library(allez)
Score=In[[1]]
names(Score)=rownames(In)
Out=allez(score=Score,lib="org.Hs.eg",idtype="SYMBOL",locallist=List)

Mat=Out$setscores[,c("Term","set.mean","set.sd","set.size","z.score")]
Mat$p.value <- pnorm(-abs(Mat$z.score))
Mat$p.adj <- p.adjust(Mat$p.value, method="BH")
Mat <- Mat[which(Mat$set.size>Lowersetsize),]

MatOut=Mat[order(Mat$p.value),c("Term","p.value","p.adj","z.score","set.size","set.mean","set.sd")]

message("sets with size < ",Lowersetsize, " are not considered" )
LocalOut=MatOut[which(is.na(MatOut[,"Term"])),]
write.table(MatOut,file=paste0(prefix,"_enrichment_allsets.txt"),sep="\t")
write.table(LocalOut,file=paste0(prefix,"_enrichment_localsets.txt"), sep="\t")



