%BME 301 TERM PROJECT MATLAB CODE
%Reparameterized Hodgkin-Huxley Model

%THIS IS THE MASTER SIMULATION FILE: IT CALLS THE FUNCTION FOR SOLVING

%Ask for n-exponent first to know what model we are using...
%Accepts an array...n_exponent is passed to the function


clear all

n_exponent = input('Enter the exponent you would like to use for n: ');

for n = 1:length(n_exponent)

%Let's get the initial conditions when I_app = 0 (i.e. no applied current) 
%and the model runs for 5 msec:
 
I_app = 0;

%Run solver
[t,x] = ode15s(@BME301_TermProject_Reparameterized_HH_ode_function,[0 5],[0 0 0 0],[],I_app,n_exponent(n));

%Set/save the state variables x0 to the values of x at 30 msec:
x0(1) = x(end,1); %v
x0(2) = x(end,2); %m
x0(3) = x(end,3); %n 
x0(4) = x(end,4); %h

%Set the applied voltage to 6.2 mV.
I_app = 6.2;

%Run solver
[t,x] = ode15s(@BME301_TermProject_Reparameterized_HH_ode_function,[0 5],x0,[],I_app,n_exponent(n));


%Add 60mV from voltage variable and plot it vs time
%+60 mV is necessary to make the baseline of the AP around 0 mV.
%This is useful for comparing results with original paper


x_zero_baseline = x(:,1)+60;

plot(t,x_zero_baseline)

%Label plot
xlabel('Time (msec)')
ylabel('Voltage (mV)')
title('Action Potential Voltage vs. Time')

hold on 

legend('Original Hodgkin-Huxley Model', 'Reparameterized Hodgkin-Huxley Model')


%Use the function to solve for the gating capacitance and current at each time point given 
%the values of v, m, n and h that was just solved for.

if n_exponent(n)==6

for i = 1:length(x)
 [f varargout(i,:)] = BME301_TermProject_Reparameterized_HH_ode_function(0,x(i,:),I_app,n_exponent(n));
end

figure

plot(t,varargout(:,4).*10^6,'r')
title('Gating Capacitance of Reparameterized Model vs Time')
xlabel('Time (msec)')
ylabel('Gating capacitance (uF/cm^{2})')
legend('Gating Capacitance')

figure

plot(t,varargout(:,5).*10^(-6),'b')
title('Gating Current of Reparameterized Model vs Time')
xlabel('Time (msec)')
ylabel('Gating Current (mA/cm^{2})')
legend('Gating Current')
xlim([0 1*10^-11])


end


end

