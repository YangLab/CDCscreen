x<-read.table("04_join_FC_negP_FPBcirc.txt")
CDCscreen=scale(-log10(x$V3))+scale(abs(log2(x$V2)))
data=cbind(x,CDCscreen)
write.table(data,file="05_CDCscreen_score.txt",sep ="\t",row.names =FALSE,col.names=FALSE, quote =FALSE)
