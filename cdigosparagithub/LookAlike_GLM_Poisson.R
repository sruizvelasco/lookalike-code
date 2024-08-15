##################################################
### Look Alike 
### Generalize Linear Models (GLM)
### Simulations 
### Poisson
##################################################

##################################################

##################################################
### Metodo 1: Estimator de Rao-Blackwell
##################################################

### Inverse of the cdf 
F1_poisson <- function(zh,tjh,j){
	A = sum(zh[1:(j-1)])
	w = rbinom(1, size=tjh, prob=1/(A+1))
	w
}

LA.RaoBlackwell.Poisson <- function(Yord,Z, Cord, n,k){
### Yla = y^* = look alike sample
### Tla = t^* = sufficient statistic from look alike sample
Yla = Yord
Tno = t(Z)%*%Yord
Tla = matrix(NA,n,k+1)
Tla[n,] = Tno

### Generate Yla_n from F^{-1}_n, recalculate Tla_{n-1} from Tno_n and Yla_n
h = Cord[n]
### F1_poisson(zh,tjh,j)
Yla[n] = F1_poisson( Z[,h], Tla[n,h], n) 
Tla[n-1,] = Tla[n,]
Tla[n-1,h] = Tla[n,h]-Yla[n]

### Generate Yla_{n-1} from F^{-1}_{n-1}, recalculate Tla_{n-2} from Tno_{n-1} and Yla_{n-1}
### Keeps going until generate Yla_{k+1} from F^{-1}_{k+1}, recalculate Tla_k from Tno_{k+1} and Yla_{k+1}
for(i in (n-1):(k+2)){
	h = Cord[i]
	### F1_poisson(zh,tjh,j)
	Yla[i] = F1_poisson( Z[,h], Tla[i,h], i) 
	Tla[i-1,] = Tla[i,]
	Tla[i-1,h] = Tla[i,h]-Yla[i]
}

### Calculate Yla_1,...,Yla_k from the unique order solution of Tno
Yla[1:(k+1)] = Tla[k+1,]

return(Yla)
}

##################################################

DevPoisson <- function(Y,Yfit){
	dd1 = Y*log(Y/Yfit) 
	dd2 = (Y-Yfit) 
	p3 = 1/sqrt(Yfit)
	dev = 2*sum(dd1[is.finite(dd1)]) - 2*sum(dd2[is.finite(dd2)]) - sum(p3[is.finite(p3)])/6
	dev
}

X2Poisson <- function(Y,Yfit){
	xx = ((Y-Yfit)^2) /Yfit
	sum(xx[is.finite(xx)])
}

##################################################

##################################################
### To generate data 

nn = c(15, 30, 50, 100, 200)
kk = c(5, 10, 15, 25, 50)
i2k = c(1,2,3,4,5)

BB = array(NA,dim=c(5,50+1))
BB[1,1:6] = c((5:1)/5, 1)
BB[2,1:11] = c((5:1)/5, (5:1)/5, 1)
BB[3,1:16] = c((5:1)/5, (5:1)/5, (5:1)/5, 1)
BB[4,1:26] = c((5:1)/5, (5:1)/5, (5:1)/5, rep(0.3,10), 1)
BB[5,1:51] = c((5:1)/5, (5:1)/5, (5:1)/5, rep(0.3,10), rep(0.2,25),1)

Xdata = array(NA,dim=c(length(nn),length(kk),  max(nn),max(kk)))
Cdata = array(NA,dim=c(length(nn),length(kk), max(nn)))

Zdata = array(NA,dim=c(length(nn),length(kk), max(nn), max(kk)+1))
Tdata = array(NA,dim=c(length(nn),length(kk), max(nn), max(kk)+1))

set.seed(12345)
for(i1 in 1:length(nn)){
for(i2 in 1:i2k[i1]){
	n = nn[i1]
	k = kk[i2]

	Xc = rmultinom(n-k-1, size=1, prob=rep(1/(k+1),(k+1)))
	X = rbind(diag(1,k), 0, t(Xc[1:k,]))
	Xdata[i1,i2,1:n,1:k] = X 

	cat = X%*%as.matrix((1:k))
	cat[cat==0] = k+1
	Cdata[i1,i2,1:n] = cat

	Z = matrix(0,n,k+1)
	for(i in 1:n){
	for(h in 1:k){
		Z[i,h] = X[i,h]*prod(1-X[i,-c(h,k+1)])
	}
	Z[i,k+1] = prod(1-X[i,-c(k+1)])
	}		
	Zdata[i1,i2,1:n,1:(k+1)] = Z 
	
}	}

##################################################
### change address
dir = "~/Documents/Articulos/Look Alike/Code GitHub/Poisson/LApois"

##################################################
### Simulate Rao-Blackwell

par(mfrow=c(3,1))
 
SIM = 100

DEV = array(NA,dim=c(SIM, length(nn),length(kk)))
DEVADJ = array(NA,dim=c(SIM, length(nn),length(kk)))
XX2 = array(NA,dim=c(SIM, length(nn),length(kk)))

Pdev = array(NA,dim=c(SIM, length(nn),length(kk)))
Pdevadj = array(NA,dim=c(SIM, length(nn),length(kk)))
Pxx2 = array(NA,dim=c(SIM, length(nn),length(kk)))

DEV.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
DEVADJ.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
XX2.LA = array(NA,dim=c(SIM, length(nn),length(kk)))

DEVsd.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
DEVADJsd.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
XX2sd.LA = array(NA,dim=c(SIM, length(nn),length(kk)))

Pdev.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
Pdevadj.LA = array(NA,dim=c(SIM, length(nn),length(kk)))
Pxx2.LA = array(NA,dim=c(SIM, length(nn),length(kk)))

M = 1000

#set.seed(12345)

for(i1 in 1:length(nn)){ 
for(i2 in 1:i2k[i1]){

n = nn[i1]
k = kk[i2]	
X = cbind(Xdata[i1,i2,1:n,1:k],1)  		
B = as.vector(BB[i2,1:(k+1)])
cat = Cdata[i1,i2,1:n]  		
Z = Zdata[i1,i2,1:n,1:(k+1)] 

Ymu = exp(X%*%B)

for(sim in 1:SIM){ 
### Data

Y = matrix(NA,n)
for(i in 1:n){
	Y[i] = rpois(1,lambda=Ymu[i])
}

### Fit glm model

mod <- glm(Y ~ 0+factor(cat), family=poisson(link="log"))
#summary(mod)

devY = deviance(mod)   
devadjY = DevPoisson(Y,mod$fit)    
xx2Y = X2Poisson(Y,mod$fit) 

DEV[sim,i1,i2] = devY 
DEVADJ[sim,i1,i2] = devadjY 
XX2[sim,i1,i2] = xx2Y 
 
Pdev[sim,i1,i2] = pchisq( devY, df=n-k-1, ncp=0, lower.tail=FALSE)
Pdevadj[sim,i1,i2] = pchisq( devadjY, df=n-k-1, ncp=0, lower.tail=FALSE)
Pxx2[sim,i1,i2] = pchisq( xx2Y, df=n-k-1, ncp=0, lower.tail=FALSE)

### Generate Look Alike
### Rao-Blackwell
devla = devadjla = xx2la = matrix(NA,M)
YLA = array(NA,dim=c(M,n))

for(mla in 1:M){
	
Yla = LA.RaoBlackwell.Poisson(Y,Z, cat, n,k)

modla <- glm(Yla ~ factor(cat), family=poisson(link="log"))

devla[mla] = deviance(modla) 
devadjla[mla] = DevPoisson(Yla,modla$fit)  
xx2la[mla] = X2Poisson(Yla,modla$fit) 

YLA[mla,] = Yla
}


DEV.LA[sim,i1,i2] = mean(devla)
DEVADJ.LA[sim,i1,i2] = mean(devadjla)
XX2.LA[sim,i1,i2] = mean(xx2la)

DEVsd.LA[sim,i1,i2] = var(devla)
DEVADJsd.LA[sim,i1,i2] = var(devadjla)
XX2sd.LA[sim,i1,i2] = var(xx2la)

Pdev.LA[sim,i1,i2] = sum(ifelse(devla>=devY, 1,0))/M
Pdevadj.LA[sim,i1,i2] = sum(ifelse(devadjla>=devadjY, 1,0))/M
Pxx2.LA[sim,i1,i2] = sum(ifelse(xx2la>=xx2Y, 1,0))/M

}	

write.table( DEV ,  paste(dir," dev.csv",sep=""))
write.table( DEVADJ ,  paste(dir," devadj.csv",sep=""))
write.table( XX2 ,  paste(dir," xx2.csv",sep=""))

write.table( Pdev ,  paste(dir," Pdev.csv",sep=""))
write.table( Pdevadj ,  paste(dir," Pdevadj.csv",sep=""))
write.table( Pxx2 ,  paste(dir," Pxx2.csv",sep=""))

write.table( DEV.LA ,  paste(dir," dev.LA.csv",sep=""))
write.table( DEVADJ.LA ,  paste(dir," devadj.LA.csv",sep=""))
write.table( XX2.LA ,  paste(dir," xx2.LA.csv",sep=""))

write.table( DEVsd.LA ,  paste(dir," devsd.LA.csv",sep=""))
write.table( DEVADJsd.LA ,  paste(dir," devadjsd.LA.csv",sep=""))
write.table( XX2sd.LA ,  paste(dir," xx2sd.LA.csv",sep=""))

write.table( Pdev.LA ,  paste(dir," Pdev.LA.csv",sep=""))
write.table( Pdevadj.LA ,  paste(dir," Pdevadj.LA.csv",sep=""))
write.table( Pxx2.LA ,  paste(dir," Pxx2.LA.csv",sep=""))

}	}

##################################################

hist(devla,probability=TRUE,nclass=50)
lines((1:3000)/10,dchisq((1:3000)/10,n-k-1),col="blue")
lines((1:3000)/10,dnorm((1:3000)/10, n-k-1,sqrt(2*(n-k-1))),col="red")
points(devY,0,col="red",pch=19)

hist(devadjla,probability=TRUE,nclass=50)
lines((1:3000)/10,dchisq((1:3000)/10,n-k-1),col="blue")
lines((1:3000)/10,dnorm((1:3000)/10, n-k-1,sqrt(2*(n-k-1))),col="red")
points(devadjY,0,col="red",pch=19)

hist(xx2la,probability=TRUE,nclass=50)
lines((1:3000)/10,dchisq((1:3000)/10,n-k-1),col="blue")
lines((1:3000)/10,dnorm((1:3000)/10, n-k-1,sqrt(2*(n-k-1))),col="red")
points(xx2Y,0,col="red",pch=19)

##################################################

##################################################
