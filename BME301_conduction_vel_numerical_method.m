clc;
close all;
clear all;
dt=0.0001; 
dx=1; 
s=dt/dx; 
t=0:dt:10; 
x=0:dx:10;

alpha_n0=0.01*(0+10)/(-1+exp((0+10)/10)); 
beta_n0=0.125*exp(0/80); 
alpha_m0=0.1*(0+25)/(-1+exp((0+25)/10)); 
beta_m0=4*exp(0/18); 
alpha_h0=0.07*exp(0/20); 
beta_h0=(1+exp((0+30)/10))^-1;
n=zeros(length(x),length(t)); 
m=zeros(length(x),length(t)); 
h=zeros(length(x),length(t)); 
v=zeros(length(x),length(t));
psi=zeros(length(x),length(t)); 
phi=zeros(length(x),length(t));

n(:,1)=alpha_n0/(alpha_n0+beta_n0); 
m(:,1)=alpha_m0/(alpha_m0+beta_m0); 
h(:,1)=alpha_h0/(alpha_h0+beta_h0);

%-15mV impulse
for i=1:1:30000 
    %v(1,i)= -15;
    v(1,i)= -15;
end

%below threshold 
%for i=1:1:30000
%v(1,i)=-5;
%end

for i=1:1:length(x)-1
    phi(i,1)=1.5; 
    v(i+1,1)=v(i,1)+dx*phi(i,1);
end

for  j=1:1:length(t)-1
%indicator for completion (in percent)
    j/(10/dt)*100
    for i=1:1:length(x)-1
        n(i,j+1)=n(i,j) + dt*N(n(i,j),v(i,j));
        m(i,j+1)=m(i,j) + dt*M(m(i,j),v(i,j));
        h(i,j+1)=h(i,j) + dt*H(h(i,j),v(i,j));
        v(i,j+1)=v(i,j) + dt*psi(i,j);
        
        phi(i+1,j)=phi(i,j) + s*(psi(i+1,j)-psi(i,j));
        psi(i,j+1)=psi(i,j) + s*(phi(i+1,j)-phi(i,j)) - ...
            dt*F(psi(i,j),v(i,j),n(i,j),m(i,j),h(i,j));
    end
    
end
v = -v;
 
figure(1)
plot(t,v(1,:))


figure(2) 
plot(t,v(1,:),t,v(2,:),t,v(3,:),t,v(4,:),t,v(5,:),t,v(6,:),t,v(7,:),t,v(8,:), t,v(9,:),t,v(10,:))
legend('1 cm','2 cm','3 cm','4 cm','5 cm','6 cm','7 cm','8 cm','9 cm','10cm')
xlabel('Time (ms)')
ylabel('Voltage (mV)')


