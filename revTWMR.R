#library(OpenMx)

cmd_args=commandArgs()
file<-cmd_args[3]
trait<-cmd_args[4]
NeQTLsfile<-cmd_args[5]

#######function SE
SE<-function(beta,gamma,Ngwas,NeQTLs){

D<-diag(length(beta))
D_inv <- solve(D)
GCG_inv <- t(beta)%*%solve(D, beta)
GCG_inv<-solve(GCG_inv)

     df_dg <- GCG_inv %*% t(beta) %*% D_inv
df_dG <- (GCG_inv %x% (t(gamma) %*% D_inv %*% (-(beta %*% GCG_inv %*% t(beta)) %*% D_inv + diag(nrow(beta))))) + ((-t(gamma) %*% D_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% D_inv))

alp = GCG_inv %*% (t(beta)%*% D_inv %*% gamma)
Vgam = pmax(1/NeQTLs,apply((gamma - beta %*% alp),2,var)%x%rep(1,length(gamma[,1])))
SEs<-c(1/sqrt(Ngwas),sqrt(Vgam))

     J <- cbind(df_dG, df_dg)

 
 R<-diag(length(beta[1,])+1)
 Sigma <- (SEs %*% t(SEs)) * (D %x% R)   
     V <- J %*% Sigma %*% t(J)
se<- sqrt(V[1,1])
return(se)
}



NeQTLs<-read.table(NeQTLsfile,header=F,sep=" ",dec=".")

out<-c("gene","alpha_original","SE_original","P_original","N_original","Phet_original","alpha","SE","P","Nend","Phet")


filecluster<-read.table(file,header=T,sep=" ",dec=".")
for (i in 2:(length(filecluster[1,])-3)) {

gamma<-as.matrix(filecluster[,i])
beta<-as.matrix(filecluster[,(length(filecluster[1,])-2)])
rsname<-filecluster[,1]
Ngwas<-as.matrix(filecluster[,length(filecluster[1,])])

nozero<-which(!is.na(gamma[,1]) & abs(beta[,1])>abs(gamma[,1]))
gamma<-as.matrix(gamma[nozero,])
beta<-as.matrix(beta[nozero,])
Ngwas<-as.matrix(Ngwas[nozero,])

rsname<-rsname[nozero]
D<-diag(length(beta))

S <- t(beta)%*%solve(D, beta) 
alpha <- solve(S, t(beta) %*% solve(D, gamma))
alpha<-as.vector(alpha)

se<-SE(beta,gamma,Ngwas,NeQTLs[i,2])

Z<-alpha[1]/se
pval<-2*pnorm(abs(Z),lower.tail=FALSE)

alphaORIGINAL<-alpha[1]
seORIGINAL<-se
pvalORIGINAL<-pval

Nstart=length(gamma[,1])
d<-gamma - alpha[1]*beta
var_d_vec<- (1/NeQTLs[i,2])+se*se*beta*beta+alpha[1]*alpha[1]*(1/Ngwas)+(1/Ngwas)*se*se
var_d<-diag(as.vector(var_d_vec))


z<-t(as.matrix(d))%*%solve(as.matrix(var_d))%*%as.matrix(d)
N<-length(d)
p_hetORIGINAL<- 1-pchisq(z,N-1)

d<-abs(d)
outlier<-which(d==max(d))[1]
p_het<-p_hetORIGINAL

while (p_het<0.05 & N>3) {
rsname<-rsname[-outlier]
gamma<-as.matrix(gamma[-outlier,])
beta<-as.matrix(beta[-outlier,])
D<-D[-outlier,-outlier]
Ngwas<-as.matrix(Ngwas[-outlier,])

S <- t(beta)%*%solve(D, beta) 
alpha <- solve(S, t(beta) %*% solve(D, gamma))

alpha<-as.vector(alpha)
se<-SE(beta,gamma,Ngwas,NeQTLs[i,2])

Z<-alpha[1]/se
pval<-2*pnorm(abs(Z),lower.tail=FALSE)

d<-gamma - alpha[1]*beta
var_d_vec<- (1/NeQTLs[i,2])+se*se*beta*beta+alpha[1]*alpha[1]*(1/Ngwas)+(1/Ngwas)*se*se
var_d<-diag(as.vector(var_d_vec))


z<-t(as.matrix(d))%*%solve(as.matrix(var_d))%*%as.matrix(d)
N<-length(d)
p_het<- 1-pchisq(z,N-1)
d<-abs(d)
outlier<-which(d==max(d))[1]

}  #end while


line<-c(colnames(filecluster)[i],round(alphaORIGINAL,digits=5),round(seORIGINAL,digits=5),signif(pvalORIGINAL,digits=3),Nstart,signif(p_hetORIGINAL,digits=3),round(alpha[1],digits=5),round(se,digits=5),signif(pval,digits=3),N,signif(p_het,digits=3))
out<-rbind(out,line)
print(i)
}



write.table(out,file=paste(trait,"_revTWMR_results.txt",sep=""),quote=F,col.names=F,row.names=F)

