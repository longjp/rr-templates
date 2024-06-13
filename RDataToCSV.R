## this code outputs elements of templates
## into .csv files which can be used without R

rm(list=ls())


load("template_des.RData")
unlink("template_des/",recursive=TRUE)
dir.create("template_des/")
write.table(data.frame(names(tem_des$dust),unname(tem_des$dust)),
            file="template_des/dust.csv",row.names=FALSE,col.names=FALSE,sep=",")
write.table(tem_des$templates,file="template_des/templates.csv",
            col.names=FALSE,sep=",")
write.csv(tem_des$betas,file="template_des/betas.csv")




load("template_sdss.RData")
unlink("template_sdss/",recursive=TRUE)
dir.create("template_sdss/")
write.table(data.frame(names(tem_sdss$dust),unname(tem_sdss$dust)),
            file="template_sdss/dust.csv",row.names=FALSE,col.names=FALSE,sep=",")
write.table(tem_sdss$templates,file="template_sdss/templates.csv",
            col.names=FALSE,sep=",")
write.csv(tem_sdss$betas,file="template_sdss/betas.csv")