## Calculating PID

## simplest case - creating a matrix of all combinations
df2 = matrix(, nrow = 8, ncol = 3)
n=1
firstvec = c(1,2)
secondvec = c(1,2)
thirdvec = c(1,2)

for(i in 1:2){
	for(j in 1:2){
		for(k in 1:2){
			df2[n,] = c(firstvec[i],secondvec[j],thirdvec[k])
			n = n + 1 	
		}
	}
}

## sample data set, hexaploid, 7 alleles

rw1freq1 = c(4,3,5,4,8,3,7)
rw1freq2 = rw1freq1/34 
names(rw1freq2) = c("rw1.104","rw1.108","rw1.112","rw1.116","rw1.120","rw1.124","rw1.128")
alnamevec = c("rw1.104","rw1.108","rw1.112","rw1.116","rw1.120","rw1.124","rw1.128")

df1 <- matrix(, nrow = 117649, ncol = 12)
n=1
for(i in 1:length(rw1freq2)){
	for(j in 1:length(rw1freq2)){
		for(k in 1:length(rw1freq2)){
			for(l in 1:length(rw1freq2)){
				for(m in 1:length(rw1freq2)){
					for(o in 1:length(rw1freq2)){
						df1[n,] = c(rw1freq2[i],rw1freq2[j],rw1freq2[k],rw1freq2[l],rw1freq2[m],rw1freq2[o],alnamevec[i],alnamevec[j],alnamevec[k],alnamevec[l],alnamevec[m],alnamevec[o])
			n = n + 1
					} 	
				}
			}
		}
	}
}

unique(df1[117648,7:12]) ## list alleles 

unique(df1[117648,7:12])[order(unique(df1[117648,7:12]))] ##reorder alleles

phenotype = paste(unique(df1[117648,7:12])[order(unique(df1[117648,7:12]))],collapse=",") ## create allele name
i
phenotype = rep(2,117649)
for (i in 1:length(phenotype)){
	phenotype[i] = paste(unique(df1[i,7:12])[order(unique(df1[i,7:12]))],collapse=",")
}

##checking
tail(df1)
tail(phenotype)
df1[10246,]
phenotype[10246]

df1 = as.data.frame(df1)

## adding columns to df1

##converting frequencies from factor to numeric

df1[,1] = as.numeric(as.character(df1[,1]))
df1[,2] = as.numeric(as.character(df1[,2]))
df1[,3] = as.numeric(as.character(df1[,3]))
df1[,4] = as.numeric(as.character(df1[,4]))
df1[,5] = as.numeric(as.character(df1[,5]))
df1[,6] = as.numeric(as.character(df1[,6]))

## finding the probability of each genotype
prod = apply((df1[,1:6]),1,prod)
head(prod)
length(prod)
df1$prob = prod
df1$pheno = phenotype
head(df1)

## summing the probability of each genotype in a new table

phenolist = unique(phenotype)
sumprobs = rep(2,length(phenolist))

for (i in 1:length(sumprobs)){
	sumprobs[i] = sum(df1$prob[which(df1$pheno==phenolist[i])])
}

sumprobs2d = sumprobs^2
sum(sumprobs2d)

# abbreviated method for calculating pid:  using example set, sample with replacement 200,000 values (100,000 sets of 2 trees), then see how many are identical -----

rw1freq2
exsamplevec = sample(alnamevec,600000,replace=TRUE,prob=rw1freq2)

## Code to count matches of two draws

dfex2 = matrix(exsamplevec,nrow=50000,ncol=12)
head(dfex2)

phenotype1 = rep(2,50000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfex2[i,1:6])[order(unique(dfex2[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,50000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfex2[i,7:12])[order(unique(dfex2[i,7:12]))],collapse=",")
}

match = rep(2,50000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/50000

## simple case to check code: triploid with 2 alleles

freq = c(0.6,0.4)
names(freq) = c("a","b")
alnamevec2 = c("a","b")

df3 <- matrix(, nrow = 8, ncol = 6)
n=1
for(i in 1:length(freq)){
	for(j in 1:length(freq)){
		for(k in 1:length(freq)){
						df3[n,] = c(freq[i],freq[j],freq[k],alnamevec2[i],alnamevec2[j],alnamevec2[k])
			n = n + 1
		}
	}
}

phenotype = rep(2,8)
for (i in 1:length(phenotype)){
	phenotype[i] = paste(unique(df3[i,4:6])[order(unique(df3[i,4:6]))],collapse=",")
}

df3 = as.data.frame(df3)

## adding columns to df1

##converting frequencies from factor to numeric

df3[,1] = as.numeric(as.character(df3[,1]))
df3[,2] = as.numeric(as.character(df3[,2]))
df3[,3] = as.numeric(as.character(df3[,3]))

## finding the probability of each genotype
prod = apply((df3[,1:3]),1,prod)
head(prod)
length(prod)
df3$prob = prod
df3$pheno = phenotype
0.4^3

## summing the probability of each genotype in a new table

phenolist = unique(phenotype)
sumprobs = rep(2,3)

for (i in 1:length(sumprobs)){
	sumprobs[i] = sum(df3$prob[which(df3$pheno==phenolist[i])])
}

sumprobs2d = sumprobs^2
sum(sumprobs2d)

##seq18d73

seq18d73freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq18d73.txt",sep = "\t",header = TRUE)
seq18d73freq = as.numeric(seq18d73freqtab)
seq18d73namevec = as.character(names(seq18d73freqtab))
length(seq18d73freq)

df4 <- matrix(, nrow = 4826809, ncol = 12)
freqvec = seq18d73freq
namevec = seq18d73namevec
n=1
for(i in 1:length(freqvec)){
	for(j in 1:length(freqvec)){
		for(k in 1:length(freqvec)){
			for(l in 1:length(freqvec)){
				for(m in 1:length(freqvec)){
					for(o in 1:length(freqvec)){
						df4[n,] = c(freqvec[i],freqvec[j],freqvec[k],freqvec[l],freqvec[m],freqvec[o],namevec[i],namevec[j],namevec[k],namevec[l],namevec[m],namevec[o])
			n = n + 1
					} 	
				}
			}
		}
	}
}

phenotype = rep(2,4826809)
for (i in 1:length(phenotype)){
	phenotype[i] = paste(unique(df4[i,7:12])[order(unique(df4[i,7:12]))],collapse=",")
}

##checking
tail(df4)
tail(phenotype)
df4[10246,]
phenotype[10246]

df4 = as.data.frame(df4)

## adding columns to df4

##converting frequencies from factor to numeric

df4[,1] = as.numeric(as.character(df4[,1]))
df4[,2] = as.numeric(as.character(df4[,2]))
df4[,3] = as.numeric(as.character(df4[,3]))
df4[,4] = as.numeric(as.character(df4[,4]))
df4[,5] = as.numeric(as.character(df4[,5]))
df4[,6] = as.numeric(as.character(df4[,6]))

## finding the probability of each genotype
prod = apply((df4[,1:6]),1,prod)
head(prod)
length(prod)
df4$prob = prod
df4$pheno = phenotype
tail(df4)

## summing the probability of each genotype in a new table

phenolist = unique(phenotype)
sumprobs = rep(2,length(phenolist))

for (i in 1:length(sumprobs)){
	sumprobs[i] = sum(df4$prob[which(df4$pheno==phenolist[i])])
}

sumprobs2d = sumprobs^2
sum(sumprobs2d)

## trying new method on seq18d73--works!

### abbreviated method for calculating pid:  using example set, sample with replacement 200,000 values (100,000 sets of 2 trees), then see how many are identical

seq18d73freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq18d73.txt",sep = "\t",header = TRUE)

seq18d73freq = as.numeric(seq18d73freqtab)
seq18d73namevec = as.character(names(seq18d73freqtab))

seq18d73samplevec = sample(seq18d73namevec,1200000,replace=TRUE,prob=seq18d73freq)

## Code to count matches of two draws

dfseq18d73 = matrix(seq18d73samplevec,nrow=100000,ncol=12)
tail(dfseq18d73)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfseq18d73[i,1:6])[order(unique(dfseq18d73[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfseq18d73[i,7:12])[order(unique(dfseq18d73[i,7:12]))],collapse=",")
}

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000

## rw39

rw39freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw39.txt",sep = "\t",header = TRUE)

rw39freq = as.numeric(rw39freqtab)
rw39namevec = as.character(names(rw39freqtab))

rw39samplevec = sample(rw39namevec,1200000,replace=TRUE,prob=rw39freq)

## Code to count matches of two draws

dfrw39 = matrix(rw39samplevec,nrow=100000,ncol=12)
tail(dfrw39)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfrw39[i,1:6])[order(unique(dfrw39[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfrw39[i,7:12])[order(unique(dfrw39[i,7:12]))],collapse=",")
}

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000

## rw28

rw28freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw28.txt",sep = "\t",header = TRUE)

rw28freq = as.numeric(rw28freqtab)
rw28namevec = as.character(names(rw28freqtab))

rw28samplevec = sample(rw28namevec,1200000,replace=TRUE,prob=rw28freq)

## Code to count matches of two draws

dfrw28 = matrix(rw28samplevec,nrow=100000,ncol=12)
tail(dfrw28)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfrw28[i,1:6])[order(unique(dfrw28[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfrw28[i,7:12])[order(unique(dfrw28[i,7:12]))],collapse=",")
}

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000

## seq8e8

seq8e8freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq8e8.txt",sep = "\t",header = TRUE)

seq8e8freq = as.numeric(seq8e8freqtab)
seq8e8namevec = as.character(names(seq8e8freqtab))

seq8e8samplevec = sample(seq8e8namevec,1200000,replace=TRUE,prob=seq8e8freq)

## Code to count matches of two draws

dfseq8e8 = matrix(seq8e8samplevec,nrow=100000,ncol=12)
tail(dfseq8e8)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfseq8e8[i,1:6])[order(unique(dfseq8e8[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfseq8e8[i,7:12])[order(unique(dfseq8e8[i,7:12]))],collapse=",")
}

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000

## rw56

rw56freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw56.txt",sep = "\t",header = TRUE)

rw56freq = as.numeric(rw56freqtab)
rw56namevec = as.character(names(rw56freqtab))

rw56samplevec = sample(rw56namevec,1200000,replace=TRUE,prob=rw56freq)

## Code to count matches of two draws

dfrw56 = matrix(rw56samplevec,nrow=100000,ncol=12)
tail(dfrw56)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfrw56[i,1:6])[order(unique(dfrw56[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfrw56[i,7:12])[order(unique(dfrw56[i,7:12]))],collapse=",")
}

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000

## rwdi11

rwdi11freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrwdi11.txt",sep = "\t",header = TRUE)

rwdi11freq = as.numeric(rwdi11freqtab)
rwdi11namevec = as.character(names(rwdi11freqtab))

rwdi11samplevec = sample(rwdi11namevec,1200000,replace=TRUE,prob=rwdi11freq)

## Code to count matches of two draws

dfrwdi11 = matrix(rwdi11samplevec,nrow=100000,ncol=12)

phenotype1 = rep(2,100000)  

for (i in 1:length(phenotype1)){
	phenotype1[i] = paste(unique(dfrwdi11[i,1:6])[order(unique(dfrwdi11[i,1:6]))],collapse=",")
}

phenotype2 = rep(2,100000)

for (i in 1:length(phenotype2)){
	phenotype2[i] = paste(unique(dfrwdi11[i,7:12])[order(unique(dfrwdi11[i,7:12]))],collapse=",")
}

tail(phenotype2)

match = rep(2,100000)
for (i in 1:length(match)){
	match[i] = ifelse(phenotype1[i]==phenotype2[i],1,0)
}

sum(match)/100000