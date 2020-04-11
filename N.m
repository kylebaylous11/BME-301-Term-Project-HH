function [N_prime] = N(n,v) 
theta=1.23138148;
alpha_n=0.01*(v+10)/(-1+exp((v+10)/10)); 
beta_n=0.125*exp(v/80);  
N_prime=(1/theta)*(alpha_n*(1-n)-beta_n*n); 
end