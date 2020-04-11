function [H_prime] = H(h,v) 
theta=1.23138148; 

alpha_h=0.07*exp(v/20); 
beta_h=(1+exp((v+30)/10))^-1; 
H_prime=(1/theta)*(alpha_h*(1-h)-beta_h*h); 
end