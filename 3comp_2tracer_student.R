install.packages(c("mvtnorm","loo","coda"), repos="https://cloud.r-project.org/",dependencies=TRUE)
options(repos=c(getOption('repos'), rethinking='http://xcelab.net/R'))
install.packages(c("coda","mvtnorm","devtools","loo"))
library(devtools)
devtools::install_github("rmcelreath/rethinking")
install.packages('ggplot2')
install.packages('bayesplot')
install.packages('openxlsx')
install.packages('reshape2')
install.packages('rstan')
library(rethinking)
library(rstan)
library(rstantools)
library(openxlsx)
library(reshape2)
library(cowplot)
library(ggplot2)
library(tidyr)
library(bayesplot)


U=read.xlsx('/Users/katie/Desktop/PCA_result.xlsx', sheet=2, colNames = TRUE)
U=as.data.frame(U)

###the first chunk is just creating some proportions and data, and then making a database for the average ###proportions times averages, including the averages (which are distributions)
N <- 22
N_weeks <- 22

## prior for PCA projected data (2H - solutes)
#event_2H <- rnorm(N,mean=2.78, sd=1.68) 
#gw_2H <- rnorm(N, mean=-1.74, sd=0.12) 
#soil_2H <- rnorm(N,  mean=-0.5, sd=0.1) 
#output_2H <- U$2H

## prior for PCA projected data (Si - isotopes)
#event_Si <- rnorm(N,mean=1.68, sd=.96) 
#gw_Si <- rnorm(N, mean=0.1, sd=0.1) 
#soil_Si <- rnorm(N, mean=-1.5, sd=0.1) 
#output_Si <- U$Si

## prior for raw data (2H)
gw_2H <- rnorm(N,mean=-72.88, sd=2.02) ##pre-event
event_2H <- rnorm(N, mean=-67.44, sd=19.4) ##event
output_2H <- db[,5]

## prior for raw data (Si )
gw_Si <- rnorm(N,mean=24.4, sd=3.5) ##pre-event
event_Si <- rnorm(N, mean=0.02, sd=0.03) ##event
output_Si <- db[,14]


d <- list(output_2H=output_2H,event_2H=event_2H, gw_2H=gw_2H,
          output_Si=output_Si,event_Si=event_Si, gw_Si=gw_Si,
          N=N)



m_code1 <- "
data{
int<lower=1> N;      //this defines the arrays for data, from 1:1000 for the output, x1, x2, and x3 data
int<lower=1> N_weeks;  //defines number of sampling days
real output_2H[N];
real event_2H[N];
real gw_2H[N];
real output_Si[N];
real event_Si[N];
real gw_Si[N];
}
parameters{
vector[22] soil_2H;
vector[22] soil_Si;
vector[3] a;
matrix[3,22] a_week;
real<lower=0> sigma;
vector<lower=0>[3] sigma_week;                          // variance of sampling across days
}

transformed parameters{
vector[3] p;
matrix[3,22] p_temp;
matrix[3,22] p_week;
p = softmax(a);         //this transforms the alpha vector to be a simplex, meaning they have to sum to 1! for main p
for ( k in 1:N_weeks ) {
p_temp[,k] = softmax(a + a_week[,k]);
p_week[,k] = p_temp[,k] - p;         // gets a p on the scale of real probability
}
}

model{
vector[N] mu_2H;
vector[N] mu_Si;
sigma ~ exponential( 1 );
sigma_week ~ exponential( 1 );
soil_2H ~ normal(-57.3, 20);
soil_Si ~ normal(6.1, 3);
a[3] ~ normal( 0 , 1);
a[2] ~ normal( 0 , 1 );
a[1] ~ normal( 0 , .001 );
a_week[3,] ~ normal(0, sigma_week[3]) ;
a_week[2,] ~ normal(0, sigma_week[2]) ;
a_week[1,] ~ normal(0, sigma_week[1]) ;

for ( i in 1:N ) {
mu_Si[i] = ((p[1] + p_week[1,i]) * event_Si[i]  + (p[2] + p_week[2,i]) * gw_Si[i] + (p[3] + p_week[3,i]) * soil_Si[i]);
mu_2H[i] = ((p[1] + p_week[1,i]) * event_2H[i]  + (p[2] + p_week[2,i]) * gw_2H[i] + (p[3] + p_week[3,i]) * soil_2H[i]);
}

output_2H ~ normal( mu_2H , sigma );
output_Si ~ student_t(20, mu_Si , sigma );
}

generated quantities{
vector[N] mu_2H;
vector[N] mu_Si;
for ( i in 1:N ) {
mu_Si[i] = ((p[1] + p_week[1,i]) * event_Si[i]  + (p[2] + p_week[2,i]) * gw_Si[i] + (p[3] + p_week[3,i]) * soil_Si[i]);
mu_2H[i] = ((p[1] + p_week[1,i]) * event_2H[i]  + (p[2] + p_week[2,i]) * gw_2H[i] + (p[3] + p_week[3,i]) * soil_2H[i]);
}
}
"



m2 <- stan( model_code=m_code1 , data=d , chains=4 , cores=4, control = list(adapt_delta=0.99))

# posterior ---------------------------------------------------------------
##script to convert wide to long data for annual proportions

p_matrix_long=matrix(NA,nrow=88000, ncol=3)
p_temp_1=read.xlsx('/Users/katie/Dropbox/J Hydrology/p_temp.xlsx', sheet = 1, colNames = TRUE)
p_temp_1=as.data.frame(p_temp_1)
p_temp_1_long <- gather(p_temp_1, date, factor_key=TRUE)
p_temp_2=read.xlsx('/Users/katie/Dropbox/J Hydrology/p_temp.xlsx', sheet = 2, colNames = TRUE)
p_temp_2=as.data.frame(p_temp_2)
p_temp_2_long <- gather(p_temp_2, date, factor_key=TRUE)
p_temp_3=read.xlsx('/Users/katie/Dropbox/J Hydrology/p_temp.xlsx', sheet = 3, colNames = TRUE)
p_temp_3=as.data.frame(p_temp_3)
p_temp_3_long <- gather(p_temp_3, date, factor_key=TRUE)
p_matrix_long[,1]=p_temp_1_long$value
p_matrix_long[,2]=p_temp_2_long$value
p_matrix_long[,3]=p_temp_3_long$value
p_matrix_long=as.data.frame(p_matrix_long)
mean(p_matrix_long[,1])
mean(p_matrix_long[,2])
mean(p_matrix_long[,3])
p_temp_4=read.xlsx('/Users/katie/Dropbox/J Hydrology/p_temp.xlsx', sheet = 4, colNames = TRUE)
p_temp_4=as.data.frame(p_temp_4)

# plot posterior for each mixture component
post <- extract.samples(m2)
p=as.data.frame(post$p)
mean(post$p[,1])
mean(post$p[,2])
mean(post$p[,3])

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")

mcmc_areas(p_temp_4,
           pars = c("V1","V3","V2"),
           prob = 0.95) + plot_title + xlim(0,1)

setEPS()
postscript('annual_p.eps', width=5, height=5)
boxplot( p_matrix_long, ylim=c(0,1), col=c('blue', 'red', 'purple'))
ggplot(data=p, aes(x=p$V1, col='blue'))+geom_density()+geom_density(data=p, aes(x=V2, col='purple'))+geom_density(data=p, aes(x=V3, col='red'))
d=density(p$V1)
d2=density(p$V2)
d3=density(p$V3)
plot(d2, xlim=0)
lines(d3)
lines(d1)
dev.off()

setEPS()
postscript('mu_2H.eps', height=4, width=8)
boxplot(post$mu_2H)
points(output_2H, col='red', pch=16)
dev.off()

setEPS()
postscript('mu_Si.eps', height=4, width=8)
boxplot(post$mu_Si)
points(output_Si, col='red', pch=16)
dev.off()

boxplot(post$soil_2H)
boxplot(post$soil_Si)

write.xlsx(post$mu_2H, '/Users/katherinemarkovich/Desktop/mu2H.xlsx')
write.xlsx(post$mu_Si, '/Users/katherinemarkovich/Desktop/muSi.xlsx')

# make p matrix

p_temp=post$p_temp
p_temp=as.data.frame(p_temp)

for (i in 1:66) {
  p_temp[4001, i]= mean(p_temp[1:4000,i])
  p_temp[4002, i]= sd(p_temp[1:4000,i])
}

write.xlsx(p_temp, '/Users/katie/Dropbox/J Hydrology/p_temp.xlsx')

# create proportion matrix ------------------------------------------------


rain_p=matrix(0,2,22)
rain_p[1,]=c(p_temp[4001,1],p_temp[4001,4],p_temp[4001,7],p_temp[4001,10],p_temp[4001,13],p_temp[4001,16],
             p_temp[4001,19],p_temp[4001,22],p_temp[4001,25],p_temp[4001,28],p_temp[4001,31],p_temp[4001,34],p_temp[4001,37],
             p_temp[4001,40],p_temp[4001,43],p_temp[4001,46],p_temp[4001,49],p_temp[4001,52],p_temp[4001,55],p_temp[4001,58],
             p_temp[4001,61],p_temp[4001,64])
rain_p=as.data.frame(rain_p)
rain_p[2,]=c(p_temp[4002,1],p_temp[4002,4],p_temp[4002,7],p_temp[4002,10],p_temp[4002,13],p_temp[4002,16],
             p_temp[4002,19],p_temp[4002,22],p_temp[4002,25],p_temp[4002,28],p_temp[4002,31],p_temp[4002,34],p_temp[4002,37],
             p_temp[4002,40],p_temp[4002,43],p_temp[4002,46],p_temp[4002,49],p_temp[4002,52],p_temp[4002,55],p_temp[4002,58],
             p_temp[4002,61],p_temp[4002,64])

gw_p=matrix(0,2,22)
gw_p[1,]=c(p_temp[4001,2],p_temp[4001,5],p_temp[4001,8],p_temp[4001,11],p_temp[4001,14],p_temp[4001,17],
           p_temp[4001,20],p_temp[4001,23],p_temp[4001,26],p_temp[4001,29],p_temp[4001,32],p_temp[4001,35],p_temp[4001,38],
           p_temp[4001,41],p_temp[4001,44],p_temp[4001,47],p_temp[4001,50],p_temp[4001,53],p_temp[4001,56],p_temp[4001,59],
           p_temp[4001,62],p_temp[4001,65])
gw_p=as.data.frame(gw_p)
gw_p[2,]=c(p_temp[4002,2],p_temp[4002,5],p_temp[4002,8],p_temp[4002,11],p_temp[4002,14],p_temp[4002,17],
           p_temp[4002,20],p_temp[4002,23],p_temp[4002,26],p_temp[4002,29],p_temp[4002,32],p_temp[4002,35],p_temp[4002,38],
           p_temp[4002,41],p_temp[4002,44],p_temp[4002,47],p_temp[4002,50],p_temp[4002,53],p_temp[4002,56],p_temp[4002,59],
           p_temp[4002,62],p_temp[4002,65])

soil_p=matrix(0,2,22)
soil_p[1,]=c(p_temp[4001,3],p_temp[4001,6],p_temp[4001,9],p_temp[4001,12],p_temp[4001,15],p_temp[4001,18],
             p_temp[4001,21],p_temp[4001,24],p_temp[4001,27],p_temp[4001,30],p_temp[4001,33],p_temp[4001,36],p_temp[4001,39],
             p_temp[4001,42],p_temp[4001,45],p_temp[4001,48],p_temp[4001,51],p_temp[4001,54],p_temp[4001,57],p_temp[4001,60],
             p_temp[4001,63],p_temp[4001,66])
soil_p=as.data.frame(soil_p)
soil_p[2,]=c(p_temp[4002,3],p_temp[4002,6],p_temp[4002,9],p_temp[4002,12],p_temp[4002,15],p_temp[4002,18],
             p_temp[4002,21],p_temp[4002,24],p_temp[4002,27],p_temp[4002,30],p_temp[4002,33],p_temp[4002,36],p_temp[4002,39],
             p_temp[4002,42],p_temp[4002,45],p_temp[4002,48],p_temp[4002,51],p_temp[4002,54],p_temp[4002,57],p_temp[4002,60],
             p_temp[4002,63],p_temp[4002,66])

all_p=matrix(0,22,6)
all_p=as.data.frame(all_p)
all_p[,1]=t(gw_p[1,])
all_p[,2]=t(rain_p[1,])
all_p[,3]=t(soil_p[1,])
all_p[,4]=t(gw_p[2,])
all_p[,5]=t(rain_p[2,])
all_p[,6]=t(soil_p[2,])


# write all_p -------------------------------------------------------------

all_p$Date=db$Date
all_p<-all_p[,c(7, 1:6)]
all_p=as.data.frame(all_p)
colnames(all_p)=c("Date", "gw", "event", "soil", "gw_sd", "event_sd", "soil_sd")
write.xlsx(all_p, '/Users/katherinemarkovich/Desktop/all_p.xlsx')
all_p=read.xlsx('/Users/katherinemarkovich/Desktop/all_p.xlsx', sheet=1, colNames = TRUE)
all_p$Date=as.Date(all_p$Date, origin='1899-12-30')
all_p=as.data.frame(all_p)

# plot average ps ---------------------------------------------------------


# plot posterior for each mixture component
boxplot( post$p , ylim=c(0,1), col=c('red',"green", 'blue'))
boxplot( post$mu)
points(output, col='red', pch=16)
# true values in red
#points( 1:3 , p_true , col="red" , pch=16 )


# plotting proportions over time ------------------------------------------


#plotting proportions over time
pd <- position_dodge(0.1)
time_p=ggplot(all_p, aes(x=Date, y=gw))+geom_line(colour='red', size=1)+ geom_errorbar(aes(ymin=all_p$gw-all_p$gw_sd, ymax=all_p$gw+all_p$gw_sd))+ylim(0,1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
#geom_errorbar(aes(ymin=all_p$gw-all_p$gw_sd, ymax=all_p$gw+all_p$gw_sd), width=.05, colour='black') +
#geom_point(colour='red')+ylim(-.1,1.1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
time_p

time_p2=ggplot(all_p, aes(x=Date, y=event))+geom_line(colour='blue', size=1)+geom_errorbar(aes(ymin=all_p$event-all_p$event_sd, ymax=all_p$event+all_p$event_sd))+ylim(0,1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
#geom_errorbar(aes(ymin=all_p$event-all_p$event_sd, ymax=all_p$event+all_p$event_sd), width=.05, colour='black') +
#geom_point(colour='blue')+ylim(-.1,1.1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
time_p2

time_p3=ggplot(all_p, aes(x=Date, y=soil))+geom_line(colour='purple', size=1)+geom_errorbar(aes(ymin=all_p$soil-all_p$soil_sd, ymax=all_p$soil+all_p$soil_sd))+ylim(0,1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
#geom_errorbar(aes(ymin=all_p$soil-all_p$soil_sd, ymax=all_p$soil+all_p$soil_sd), width=.05, colour='black') +
#geom_point(colour='purple')+ylim(-.1,1.1)+ scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y")
time_p3

setEPS()
postscript('temporal.eps', width=10, height=9)
multiplot(time_p, time_p2, time_p3)
dev.off()

all_p_no_error=all_p[,(1:4)]
melted = melt(all_p_no_error, id.vars = "Date")
area2= ggplot(melted, aes(x=Date, y=value))+geom_area(aes(colour=variable, fill=variable))+theme(legend.position='none')
area2


