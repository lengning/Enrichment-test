
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # score file, rows are genes
Local=options[2] # local list; rows are lists
Lowersetsize=as.numeric(options[3])
Uppersetsize=as.numeric(options[4])
db.v <- options[5]
if(Local=="NULL" | length(options)<2 ) Local <- NULL
if(length(options)<3)Lowersetsize=10
if(length(options)<4)Uppersetsize=800
if(length(options)<5)db.v <- "human" # annotation
message(c("annotation: ", db.v))
if(db.v=="human")lib.v <- "org.Hs.eg"
if(db.v=="mouse")lib.v <- "org.Mm.eg"

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

  if(nrow(ListIn)>1){
  List=sapply(1:nrow(ListIn),function(i)setdiff(as.vector(ListIn[i,]),c(""," ")))
  }
  if(nrow(ListIn)==1) {
    List=vector("list",1)
    List[[1]]=setdiff(unlist(ListIn),c(""," "))}
	names(List)=rownames(ListIn)

} else List=NULL

library(EACI)
Score=In[[1]]
names(Score)=rownames(In)
Out=eacitest(score=Score,lib=lib.v,idtype="SYMBOL",locallist=List,iter=10, minsetsize=Lowersetsize)

Mat=Out[[1]][,c("Term","set.mean","set.sd","set.size","pval")]
Mat$p.adj <- p.adjust(Mat$pval, method="BH")
Mat <- Mat[which(Mat$set.size>Lowersetsize),]
Mat <- Mat[which(Mat$set.size<Uppersetsize),]
MatOut=Mat[order(Mat$pval),c("Term","pval","p.adj","set.size","set.mean","set.sd")]
message("sets with size < ",Lowersetsize, " or > ", Uppersetsize, " are not considered" )


LocalOut=MatOut[which(is.na(MatOut[,"Term"])),]

MatOut2 <-  cbind(rownames(MatOut), MatOut)
LocalOut2 <- cbind(rownames(LocalOut), LocalOut)
colnames(MatOut2)[1] = colnames(LocalOut2)[1] = "GO_ID"
write.table(MatOut2,file=paste0(prefix,"_EACIenrichment_allsets.txt"),sep="\t", row.names=F)
write.table(LocalOut2,file=paste0(prefix,"_EACIenrichment_localsets.txt"), sep="\t", row.names=F)


