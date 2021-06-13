# This RNA interaction program
# Version 2021.6.13
# including Simith matrix<chen_smis(seq1,seq2)>; <SWChen(x,y,m,t,large,message)>;
# including Volcano trace<Sefmap(chen_smis)>; Scatter plot martrix <MatChen(x,y,noise,...)>
# 核心函数：执行代码SWChen(目标序列.csv,数据库序列.csv,默认值=RNA互作评分rv【还含有s\smi\mis\hair四个参数】)

library(ggplot2) # need package and core calculation theroy is not necessary

# LncRNA miRNA mRNA靶向预测及RNA区域二级结构程序，提供rv\s\smi\mis\hair五种预测形式
# default <rv>: interaction calculate
# <s>: similar align
# <smi>: seed region of miRNA sequence to match sequence database
# <mis>: RNA sequence to match miRNA database sequence
# <hair>: predict hair structure

# public functions
# sequence convert function
chen_sc=function(x,n=1){
T=as.character(x)
z=c()
for (i in 1:length(T)){
change=c()
change=unlist(strsplit(T[i],NULL))
if(n==3){
z[i]=paste(rev(change),collapse="")
}else if(n==2){
change[which(change=="A")]="T@"
change[which(change=="T")]="A@"
change[which(change=="C")]="G@"
change[which(change=="G")]="C@"
change[which(change=="A@")]="A"
change[which(change=="T@")]="T"
change[which(change=="C@")]="C"
change[which(change=="G@")]="G"
z[i]=paste(change,collapse="")
}else{
change[which(change=="A")]="T@"
change[which(change=="T")]="A@"
change[which(change=="C")]="G@"
change[which(change=="G")]="C@"
change[which(change=="A@")]="A"
change[which(change=="T@")]="T"
change[which(change=="C@")]="C"
change[which(change=="G@")]="G"
z[i]=paste(rev(change),collapse="")
}}
z}

# match score: this is very inportant for your prediction
S=function(x,y){
SM=matrix(c(2,-7,-5,-7,-7,2,-7,-5,-5,-7,2,-7,-7,-5,-7,2),4,4)
I=4;J=4
if(x=="A"){I=1}else if(x=="C"){I=2}else if(x=="G"){I=3}else{I=4}
if(y=="A"){J=1}else if(y=="C"){J=2}else if(y=="G"){J=3}else{J=4}
SorM=0;SorM=SM[I,J]
SorM}

# F matrix fill values: this is very inportant for your prediction
chen_smis=function(x,y){
if(nchar(as.character(x))>nchar(as.character(y))){A=as.character(y);B=as.character(x)
}else{A=as.character(x);B=as.character(y)}
F=matrix(0,nrow=nchar(A)+1,ncol=nchar(B)+1)
Match=0;Delete=0;Insert=0;d=-5
for(i in 1:nchar(A)){
for(j in 1:nchar(B)){
SD=S(substring(A,i,i),substring(B,j,j))
Match = F[i,j] + SD
Delete = F[i, j+1] + d
Insert = F[i+1, j] + d
F[i+1,j+1] = max(Match, Insert, Delete,0)
}};F}

# F matrix for re-socre: this is very inportant for your prediction
chen_score=function(F){
pscore=0;pscore=max(F)
num=0;num=length(which(F==pscore))
rg=which(F==pscore)%%(nrow(F))
rg[which(rg==0)]=nrow(F)
cg=ceiling(which(F==pscore)/nrow(F))
i=rg[length(rg)];j=cg[length(cg)]
n=1;miss=0;stat=1
score=0;ogap=0;egap=0
misscore=0
repeat{M=max(F[i-1,j-1],F[i-1,j],F[i,j-1],0)
pstat=stat
if(i==2|j==2|M==0){break}else if(M==F[i-1,j-1]){
i=i-1;j=j-1;stat=1
if(F[i,j]-F[i-1,j-1]>0){n=n+1}else{miss=miss+1}
}else if(M==F[i-1,j]){i=i-1;stat=0
}else{j=j-1;stat=0}
if(stat==0){egap=egap+1} 
if(pstat!=stat){ogap=ogap+1}
misscore=5*ogap+2*egap+miss*5
score=(n+1)*10-misscore
if(misscore>=50|n+1>=nrow(F)){break}} # Back count and score again
srg="";scg=""
srg=paste("第",as.character(rg),"位",sep="",collapse=",")
scg=paste("第",as.character(cg),"位",sep="",collapse=",")
best=matrix();best=matrix(c(score,num,srg,scg),ncol=4)
best}

# core map function
SWChenM=function(x,t=1,Classer="single",passname="SWChenM",RNA="RNA",Larg=300,message="Drafted by Chen"){
address=getwd()
best=x
cat("开始绘制作用位点图... ...\n")
# table making
Tmar=data.frame()
Tmar=data.frame(name=as.character(best[,1]),HighestScore=best$HighestScore,SCG=best$SCG_Rverse_longer)
data_mapM=data.frame()
for(mj in 1:nrow(Tmar)){
if(grepl(",",as.character(Tmar$SCG[mj]))==FALSE){
mapM=data.frame()
short1=as.integer(substr(as.character(Tmar$SCG[mj]),2,nchar(as.character(Tmar$SCG[mj]))-1))
mapM=data.frame(name=Tmar$name[mj],xp=short1,Yscore=as.integer(Tmar$HighestScore[mj]))
}else{
xsp=c()
xsp=unlist(strsplit(as.character(Tmar$SCG[mj]),split=","))
xsp=substr(xsp,2,nchar(xsp)-1)
mapM=data.frame()
mapM=data.frame(name=Tmar$name[mj],xp=as.integer(xsp),Yscore=as.integer(Tmar$HighestScore[mj]))}
data_mapM=rbind(data_mapM,mapM)}
#made table
data_mapM[,2]=as.integer(data_mapM[,2])
cat("绘制数据处理完毕\n")
graph=ggplot()+
geom_point(aes(x=data_mapM$xp,y=data_mapM$Yscore,color=data_mapM$Yscore,shape="★"))+
scale_colour_gradient(low = "black",high = "red")
if(t==1){
mns=which(data_mapM[,3]==max(data_mapM[,3]))
mtext=""; mx=0; my=0
mtext=as.character(data_mapM[mns,1])
mx=data_mapM[mns,2]
my=data_mapM[mns,3]
graph=graph+geom_text(aes(x=mx,y=my),label=paste("↙",substr(mtext[loop],1,15)),hjust=0,vjust=0,check_overlap = T)
}else{cat("不标注最大信息\n")}
egraph=graph+labs(title=paste("RNA INTERACTION"," [",Classer,"]     ",message,sep=""), 
x="Seqence Length", y="Interaction Score")+
theme(plot.title = element_text(hjust = 0.5))
cat("绘制完毕\n")
if (file.exists("./GraphRNAi")==TRUE){cat("阁下目标文件夹 GraphRNAi 已存在\n")}else{
dir.create("./GraphRNAi", recursive=TRUE)
cat("目标文件夹 GraphRNAi 已为阁下创建\n")}
setwd("./GraphRNAi")
gNAME=paste("阁下",passname,RNA,"已绘制完成",gsub(":","_",Sys.time()),".png")
ggsave(filename=gNAME,egraph,dpi=Larg,width=24,height=8)
cat("阁下",passname,RNA,"已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)
egraph}

# Core function
SWChen=function(x,y,Classer="rv",m=1,t=0,Larg=600,message="★SWChen★"){
address=getwd()
cat("Welcome to SWChen RNA interaction predictor\nCurrent Version： 2021.6.13 Version\n")
cat("Copyright reserved\nupdate:\nhttp://ylchen-swchen.lofter.com/ \nhttps://github.com/YLCHEN1992\nStart working\n")
TYPE=Classer
cat(paste("You selected ",Classer,"type\n",sep="")) 
X=read.csv(deparse(substitute(x)))
Y=read.csv(deparse(substitute(y)))
RNA=""; X12=""; Y02=c()
if(TYPE=="s"){
X12<-as.character(X[1,2])
Y02<-as.character(Y[,2])
RNA="RNA序列局部相似分析"
}else if(TYPE=="smi"){
X12<-chen_sc(substr(as.character(X[1,2]),2,8))
Y02<-as.character(Y[,2])
RNA="目标miRNA短种子"
}else if(TYPE=="mis"){
X12<-as.character(X[1,2])
Y02<-chen_sc(substring(as.character(Y[,2]),2,8))
RNA="miRNA数据库短种子"
}else if(TYPE=="hair"){
X12<-as.character(X[1,2])
Y02<-chen_sc(as.character(Y[,2]))
RNA="内部发夹结构预测"
}else{
X12<-chen_sc(as.character(X[1,2]))
Y02<-as.character(Y[,2])
RNA=" RNA相互作用评分 "}
N=nrow(Y)
Mg=matrix(0,nrow=N,ncol=4)
for(i in 1:N){
T=Sys.time()
Mg[i,]=chen_score(chen_smis(as.character(X12),Y02[i]))
cat("序列【",as.character(Y[i,1]),"】已完成\n",
as.character(round(as.numeric((i/N)*100),3)),"% 运算完成！\n")
TP=Sys.time()
DT=as.numeric(TP-T)
cat("%\n 当前速度为",as.character(round(as.numeric(DT),3)),"秒/个\n",
"大约还需要",as.character(round(as.numeric(DT),3)*(N-i)),"秒\n")}
bst=data.frame()
bst=data.frame(Mg)
colnames(bst)=c("HighestScore","NumberBinding","SRG","SCG_Rverse_longer")
bst=cbind(Y,bst)
short2=chartr('.',"_",as.character(X[1,1]))
NAME=paste("阁下",short2,RNA,"已预测完成",gsub(":","_",Sys.time()),".csv")
if (file.exists("./RNAi")==TRUE){cat("阁下目标文件夹 RNAi 已存在\n")}else{
dir.create("./RNAi", recursive=TRUE)
cat("目标文件夹 RNAi 已为阁下创建\n")}
setwd("./RNAi")
write.csv(bst,NAME,row.names=FALSE)
cat("阁下",short2,RNA,"已预测完成,文件保存在",as.character(getwd()),"目录下\n")
best0=read.csv(NAME)
setwd(address)
if(m==1){SWChenM(x=best0,t=t,Larg=Larg,Classer=Classer,passname=short2,RNA=RNA,message=message)
}else{cat("选择命令m=1不画图")}}

# affiliate function for simple compute 1 need binding [chen_smis] function
Sefmap=function(F,name="Interaction Region",message="Interaction Region",Larg=600){
address=getwd()
x=c(1:ncol(F)); y=c()
for(i in 1:ncol(F)){y=c(y,sum(F[,i]))}
map=data.frame(Sites=x,Score=y)
graphmap=ggplot(map)+
geom_point(aes(x=Sites,y=Score,color=Score))+
scale_colour_gradient(low = "green3",high = "black")+
labs(title=message)+theme(plot.title = element_text(hjust = 0.5))
if (file.exists("./GraphRNAi")==TRUE){cat("阁下目标文件夹 GraphRNAi 已存在\n")}else{
dir.create("./GraphRNAi", recursive=TRUE)
cat("目标文件夹 GraphRNAi 已为阁下创建\n")}
setwd("./GraphRNAi")
gNAME=paste("阁下",name,"已绘制完成",gsub(":","_",Sys.time()),".png")
ggsave(filename=gNAME,graphmap,dpi=Larg,width=24,height=4)
cat("阁下",name,"已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)
graphmap}

# affiliate function for simple compute 2
MatChen=function(x,y,noise=3,name="Interaction Point",message="Interaction Point",Larg=600){
address=getwd()
if(nchar(as.character(x))>nchar(as.character(y))){
A=as.character(y)
B=as.character(x)}else{
A=as.character(x)
B=as.character(y)}
F=matrix(0,nrow=nchar(A)+1,ncol=nchar(B)+1)
for(i in seq(1,nchar(A),noise)){
for(j in seq(1,nchar(B),noise)){
if(substring(A,i,i+noise-1)==substring(B,j,j+noise-1)){F[i+1,j+1]="MARK"}}}
X=c();Y=c()
for(i in 1:ncol(F)){
J=c()
J=nrow(F)-which(F[,i]=="MARK")
X=c(X,rep(i,length(J)))
Y=c(Y,J)}
map=data.frame(A=X,B=Y)
graphmap=ggplot(map)+
geom_point(aes(x=A,y=B), shape=17)+
labs(title=message)+theme(plot.title = element_text(hjust = 0.5))
if (file.exists("./GraphRNAi")==TRUE){cat("阁下目标文件夹 GraphRNAi 已存在\n")}else{
dir.create("./GraphRNAi", recursive=TRUE)
cat("目标文件夹 GraphRNAi 已为阁下创建\n")}
setwd("./GraphRNAi")
gNAME=paste("阁下",name,"已绘制完成",gsub(":","_",Sys.time()),".png")
ggsave(filename=gNAME,graphmap,dpi=Larg,width=24,height=4)
cat("阁下",name,"已绘制完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)
graphmap}


