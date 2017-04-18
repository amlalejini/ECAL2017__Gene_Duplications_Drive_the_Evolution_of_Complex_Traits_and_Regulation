setwd("~/Dropbox/ChangingEnvironments/")

# This analysis script includes information about some sets of runs we did not include in this
# submission: ones where we did not enforce a minimum genome size.  These are all referred to as
# no.min in their data sets.  While this would seem to indicate that we should adjust the 
# Bonferroni corrections of the Kruskal-Wallis tests to include that we only have 3 environments,
# not 6, the reality is that this correction for 3 would make the p values smaller, and they are
# already so small that we report them as << 0.0001, rather than their exact values.

#data.raw <- read.table(file="SlipMutations_FinalDoms.csv", header=T, sep=",")
data.raw <- read.table(file="SlipMutations_FinalDoms__DUPvSCRAM.csv", header=T, sep=",")

Q1.data <- subset(data.raw, question=="Q1")
Q2.data <- subset(data.raw, question=="Q2")
Q3.data <- subset(data.raw, question=="Q3")

Q1.data.no.min <- subset(Q1.data, min_genome_size=="0")
Q1.data.min <- subset(Q1.data, min_genome_size=="100")

Q2.data.no.min <- subset(Q2.data, min_genome_size=="0")
Q2.data.min <- subset(Q2.data, min_genome_size=="100")

Q3.data.no.min <- subset(Q3.data, min_genome_size=="0")
Q3.data.min <- subset(Q3.data, min_genome_size=="100")

############################# This section just for the additional data set DUPvSCRAM:

Q1.data.F <- subset(Q1.data, fancy_name=="Slip-scramble")
Q1.data.G <- subset(Q1.data, fancy_name=="Slip-duplicate")

Q2.data.F <- subset(Q2.data, fancy_name=="Slip-scramble")
Q2.data.G <- subset(Q2.data, fancy_name=="Slip-duplicate")

Q3.data.F <- subset(Q3.data, fancy_name=="Slip-scramble")
Q3.data.G <- subset(Q3.data, fancy_name=="Slip-duplicate")

p.Q1.min.7 <- wilcox.test(x=Q1.data.F$fdom_phenotype_score, y=Q1.data.G$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]
p.Q2.min.7 <- wilcox.test(x=Q2.data.F$fdom_phenotype_score, y=Q2.data.G$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]
p.Q3.min.7 <- wilcox.test(x=Q3.data.F$fdom_phenotype_score, y=Q3.data.G$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]

p.SCRAM.DUP <- p.adjust(c(p.Q1.min.7, p.Q2.min.7, p.Q3.min.7), method="bonferroni", n=3)

#################


Q1.min.A <- subset(Q1.data.min, fancy_name=="Baseline")
Q1.min.B <- subset(Q1.data.min, fancy_name=="High mutation rate")
Q1.min.C <- subset(Q1.data.min, fancy_name=="Slip-scatter")
Q1.min.D <- subset(Q1.data.min, fancy_name=="Slip-NOP")
Q1.min.E <- subset(Q1.data.min, fancy_name=="Slip-random")
Q1.min.F <- subset(Q1.data.min, fancy_name=="Slip-scramble")
Q1.min.G <- subset(Q1.data.min, fancy_name=="Slip-duplicate")

Q1.no.min.A <- subset(Q1.data.no.min, fancy_name=="Baseline")
Q1.no.min.B <- subset(Q1.data.no.min, fancy_name=="High mutation rate")
Q1.no.min.C <- subset(Q1.data.no.min, fancy_name=="Slip-scatter")
Q1.no.min.D <- subset(Q1.data.no.min, fancy_name=="Slip-NOP")
Q1.no.min.E <- subset(Q1.data.no.min, fancy_name=="Slip-random")
Q1.no.min.F <- subset(Q1.data.no.min, fancy_name=="Slip-scramble")
Q1.no.min.G <- subset(Q1.data.no.min, fancy_name=="Slip-duplicate")

Q2.min.A <- subset(Q2.data.min, fancy_name=="Baseline")
Q2.min.B <- subset(Q2.data.min, fancy_name=="High mutation rate")
Q2.min.C <- subset(Q2.data.min, fancy_name=="Slip-scatter")
Q2.min.D <- subset(Q2.data.min, fancy_name=="Slip-NOP")
Q2.min.E <- subset(Q2.data.min, fancy_name=="Slip-random")
Q2.min.F <- subset(Q2.data.min, fancy_name=="Slip-scramble")
Q2.min.G <- subset(Q2.data.min, fancy_name=="Slip-duplicate")

Q2.no.min.A <- subset(Q2.data.no.min, fancy_name=="Baseline")
Q2.no.min.B <- subset(Q2.data.no.min, fancy_name=="High mutation rate")
Q2.no.min.C <- subset(Q2.data.no.min, fancy_name=="Slip-scatter")
Q2.no.min.D <- subset(Q2.data.no.min, fancy_name=="Slip-NOP")
Q2.no.min.E <- subset(Q2.data.no.min, fancy_name=="Slip-random")
Q2.no.min.F <- subset(Q2.data.no.min, fancy_name=="Slip-scramble")
Q2.no.min.G <- subset(Q2.data.no.min, fancy_name=="Slip-duplicate")

Q3.min.A <- subset(Q3.data.min, fancy_name=="Baseline")
Q3.min.B <- subset(Q3.data.min, fancy_name=="High mutation rate")
Q3.min.C <- subset(Q3.data.min, fancy_name=="Slip-scatter")
Q3.min.D <- subset(Q3.data.min, fancy_name=="Slip-NOP")
Q3.min.E <- subset(Q3.data.min, fancy_name=="Slip-random")
Q3.min.F <- subset(Q3.data.min, fancy_name=="Slip-scramble")
Q3.min.G <- subset(Q3.data.min, fancy_name=="Slip-duplicate")

Q3.no.min.A <- subset(Q3.data.no.min, fancy_name=="Baseline")
Q3.no.min.B <- subset(Q3.data.no.min, fancy_name=="High mutation rate")
Q3.no.min.C <- subset(Q3.data.no.min, fancy_name=="Slip-scatter")
Q3.no.min.D <- subset(Q3.data.no.min, fancy_name=="Slip-NOP")
Q3.no.min.E <- subset(Q3.data.no.min, fancy_name=="Slip-random")
Q3.no.min.F <- subset(Q3.data.no.min, fancy_name=="Slip-scramble")
Q3.no.min.G <- subset(Q3.data.no.min, fancy_name=="Slip-duplicate")



kruskal.test(fdom_phenotype_score ~ treatment, data=Q1.data.min)
p.1.1 <- kruskal.test(fdom_phenotype_score ~ treatment, data=Q1.data.min)[[3]]

p.Q1.min.1 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.B$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.min.2 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.C$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.min.3 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.D$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.min.4 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.E$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.min.5 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.F$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.min.6 <- wilcox.test(x=Q1.min.A$fdom_phenotype_score, y=Q1.min.G$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 


p.Q1.min.raw <- c(p.Q1.min.1, p.Q1.min.2, p.Q1.min.3, p.Q1.min.4, p.Q1.min.5, p.Q1.min.6)

p.Q1.min.adj <- p.adjust(p.Q1.min.raw, method="bonferroni", n=6)

########### Q1, minimum genome, exactly as we expect.  There is an overall effect
# and the only ones different from the baseline are the Slip_DUP and Slip_SCRAM

kruskal.test(fdom_phenotype_score ~ treatment, data=Q1.data.no.min)
p.1.2 <- kruskal.test(fdom_phenotype_score ~ treatment, data=Q1.data.no.min)[[3]]

p.Q1.no.min.1 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.B$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.2 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.C$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.3 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.D$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.4 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.E$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.5 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.F$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.6 <- wilcox.test(x=Q1.no.min.A$fdom_phenotype_score, y=Q1.no.min.G$fdom_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q1.no.min.raw <- c(p.Q1.no.min.1, p.Q1.no.min.2, p.Q1.no.min.3, p.Q1.no.min.4, p.Q1.no.min.5, p.Q1.no.min.6)

p.Q1.no.min.adj <- p.adjust(p.Q1.no.min.raw, method="bonferroni", n=6)

########## Things are different in Q1 with no minimum genome size.  Here,
# Slip_NOP. Slip_RAND, and Slip_SCRAM all do *worse* than baseline; Slip_DUP
# still does better.
# No, in the new data set, Slip-scramble does better.  This is different.



kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q2.data.min)
p.1.3 <- kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q2.data.min)[[3]]

p.Q2.min.1 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.B$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.2 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.C$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.3 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.D$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.4 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.E$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.5 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.F$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.6 <- wilcox.test(x=Q2.min.A$end_min_1000_phenotype_score, y=Q2.min.G$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.min.raw <- c(p.Q2.min.1, p.Q2.min.2, p.Q2.min.3, p.Q2.min.4, p.Q2.min.5, p.Q2.min.6)

p.Q2.min.adj <- p.adjust(p.Q2.min.raw, method="bonferroni", n=6)

####### For Q2, minimum genome size, Slip_DUP and Slip_SCRAM are better 
# than baseline.  HIGHMUT is *worse* than baseline.


kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q2.data.no.min)
p.1.4 <- kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q2.data.no.min)[[3]]

p.Q2.no.min.1 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.B$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.2 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.C$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.3 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.D$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.4 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.E$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.5 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.F$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.6 <- wilcox.test(x=Q2.no.min.A$end_min_1000_phenotype_score, y=Q2.no.min.G$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q2.no.min.raw <- c(p.Q2.no.min.1, p.Q2.no.min.2, p.Q2.no.min.3, p.Q2.no.min.4, p.Q2.no.min.5, p.Q2.no.min.6)

p.Q2.no.min.adj <- p.adjust(p.Q2.no.min.6, method="bonferroni", n=6)


# For Q2 with no minimum, anything that's different than baseline is *worse*
# That includes HIGHMUT, Slip_NOP, Slip_RAND, Slip_SCRAM, and Slip_DUP


############## Q3: Changing environments

kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q3.data.min)
p.1.5 <-- kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q3.data.min)[[3]]

p.Q3.min.1 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.B$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.min.2 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.C$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]  

p.Q3.min.3 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.D$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]  

p.Q3.min.4 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.E$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]  

p.Q3.min.5 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.F$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]  

p.Q3.min.6 <- wilcox.test(x=Q3.min.A$end_min_1000_phenotype_score, y=Q3.min.G$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]]  

p.Q3.min.raw <- c(p.Q3.min.1, p.Q3.min.2, p.Q3.min.3, p.Q3.min.4, p.Q3.min.5, p.Q3.min.6)

p.Q3.min.adj <- p.adjust(p.Q3.min.raw, method="bonferroni", n=6)

####### For Q3, minimum genome size, significance after adjustment for multipl
####### comparisons is: YES for High Mutation Rate (worse), 
####### YES for Slip Scramble, and Slip Duplicate (better)
####### NO for all others.


kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q3.data.no.min)
p.1.6 <- kruskal.test(end_min_1000_phenotype_score ~ treatment, data=Q3.data.no.min)[[3]]

p.Q3.no.min.1 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.B$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.2 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.C$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.3 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.D$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.4 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.E$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.5 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.F$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.6 <- wilcox.test(x=Q3.no.min.A$end_min_1000_phenotype_score, y=Q3.no.min.G$end_min_1000_phenotype_score, alternative=c("two.sided"), conf.int=TRUE)[[3]] 

p.Q3.no.min.raw <- c(p.Q3.no.min.1, p.Q3.no.min.2, p.Q3.no.min.3, p.Q3.no.min.4, p.Q3.no.min.5, p.Q3.no.min.6)

p.Q3.no.min.adj <- p.adjust(p.Q3.no.min.raw, method="bonferroni", n=6)

##### For Q3, no minimum genome size, adjusted for multiple comparisons:
##### NOT significant for Slip-scatter or Slip-duplicate
##### YES for High mutation rate, Slip-NOP, Slip-random, and Slip-scramble (worse)


adjusted.p.1 <- p.adjust(c(p.1.1, p.1.2, p.1.3, p.1.4, p.1.5, p.1.6), method="bonferroni", n=6)


graphics.off()
png(file="Test_for_Alex.png")
plot(x=jitter(Q1.no.min.A$fdom_phenotype_score), y=jitter(Q1.no.min.F$fdom_phenotype_score), xlab="Baseline", ylab="Slip-scramble", xlim=c(0, 9), ylim=c(0, 9))
abline(a=0, b=1, lty=1, lwd=2, col="black")
dev.off()
