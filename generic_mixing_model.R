#generic script for bayesian mixing model created by Katie Markovich with assistance from Brendan Barrett and Richard McElreath (Fall 2017)

###### Load Packages ######
library(rethinking)
library(rstan)
library(rstantools)

###### User Inputs #######
outlet_db=c(-9.40, -10.04, -9.31, -9.59, -9.75, -8.85, -9.59, -9.77, -9.80, -9.78, -9.85, -10.00, -9.83, -10.02, -10.21, -9.65, -9.67, -9.70, -9.75, -9.82, -9.30, -8.91)

N <- 22 # number of sampling points for endmember concentrations
N_weeks <- 22 # number of sampling events

gw_d18O <- rnorm(N,mean=-10.8, sd=0.61) ##groundwater -- change to endmember 1
precip_d18O <- rnorm(N, mean=-7.41, sd=1.65) ##rain -- change to endmember 2
outlet_d18O <- outlet_db

d <- list(outlet_d18O=outlet_d18O, gw_d18O=gw_d18O, precip_d18O=precip_d18O, N=N)



m_code1 <- "
data{
int<lower=1> N;      //this defines the arrays for data, from 1:1000 for the output, x1, x2, and x3 dara
int<lower=1> N_weeks;  //defines number of sampling days
real outlet_d18O[N];
real gw_d18O[N];
real precip_d18O[N];
}
parameters{
vector[2] a;
matrix[2,22] a_week;
real<lower=0> sigma;
vector<lower=0>[2] sigma_week;                          // variance of sampling across days
}

transformed parameters{
vector[2] p;
matrix[2,22] p_temp;
matrix[2,22] p_week;
p = softmax(a);         //this transforms the alpha vector to be a simplex, meaning they have to sum to 1! for main p
  for ( k in 1:N_weeks ) {
  p_temp[,k] = softmax(a + a_week[,k]);
  p_week[,k] = p_temp[,k]-p;         // gets a p on the scale of real probability
}
}

model{
vector[N] mu;
sigma ~ exponential( 1 );
sigma_week ~ exponential( 1 );
a[2] ~ normal( 0 , .01 ); //this
a[1] ~ normal( 0 , 1 );
a_week[2,] ~ normal(0, sigma_week[2]) ;
a_week[1,] ~ normal(0, sigma_week[1]) ;

for ( i in 1:N ) {
mu[i] = (p[1] + p_week[1,i]) * gw_d18O[i] + (p[2] + p_week[2,i]) * precip_d18O[i];
}
outlet_d18O ~ normal( mu , sigma );
}

generated quantities{
vector[N] mu;
for ( i in 1:N ) {
mu[i] = (p[1] + p_week[1,i]) * gw_d18O[i] + (p[2] + p_week[2,i]) * precip_d18O[i];
}
}
"



m1 <- stan( model_code=m_code1 , data=d , chains=4 , cores=4, control = list(adapt_delta=0.99))

# plot posterior for each mixture component
post <- extract.samples(m1)
boxplot( post$p , ylim=c(0,1) )
