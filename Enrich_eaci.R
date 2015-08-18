
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # score file, rows are genes
Local=options[2] # local list; rows are lists
if(length(options)<2)Local=NULL


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

library(EACI)
Score=In[[1]]
names(Score)=rownames(In)
Out=eacitest(score=Score,lib="org.Hs.eg",idtype="SYMBOL",locallist=List,iter=10)

Mat=Out[[1]][,c("Term","set.mean","set.sd","set.size","pval")]
Mat$p.adj <- p.adjust(Mat$pval, method="BH")
MatOut=Mat[order(Mat$pval),c("Term","pval","p.adj","set.size","set.mean","set.sd")]

LocalOut=MatOut[which(is.na(MatOut[,"Term"])),]
write.table(MatOut,file=paste0(prefix,"_EACIenrichment_allsets.txt"),sep="\t")
write.table(LocalOut,file=paste0(prefix,"_EACIenrichment_localsets.txt"), sep="\t")



