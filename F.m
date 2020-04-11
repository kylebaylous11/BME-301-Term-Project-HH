function f=F(psi,v,n,m,h) 
R2_a=2974.8991; 
theta=1.23138148;  
C=0.001;
g_k=0.036; 
g_n=0.12; 
g_l=0.0003; 
v_k=12; 
v_n=-115; 
v_l=-10.5989;
N_f=N(n,v); 
M_f=M(m,v); 
H_f=H(h,v);

f=R2_a*C*theta*psi + (1/(theta*C))*(g_k*n^4 + g_n*m^3*h + g_l)*psi + ...
    g_k*(R2_a*n^4 +(4/(theta*C))*n^3*N_f)*(v-v_k) + ...
    g_n*(R2_a*m^3*h + (1/(theta*C))*(3*m^2*h*M_f + m^3*H_f))*(v-v_n) + ...
    g_l*R2_a*(v-v_l); 
end