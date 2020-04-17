function [M_prime] = M(m,v)  
theta=1.23138148;
alpha_m=0.1*(v+25)/(-1+exp((v+25)/10)); 
beta_m=4*exp(v/18); 
M_prime=(1/theta)*(alpha_m*(1-m)-beta_m*m); 
end