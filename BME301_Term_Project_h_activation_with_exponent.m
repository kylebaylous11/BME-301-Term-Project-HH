%BME 301 TERM PROJECT MATLAB CODE
%Reparameterized Hodgkin-Huxley Model

%THIS IS THE MAIN SIMULATION FILE

%ASK FOR n-exponent

clear all

n_exponent = input('Enter the exponent you would like to use for n: ');
%b_h1 = input('Enter b_h1: ');
%b_h2 = input('Enter b_h2: ');

for n = 1:length(n_exponent)

%Let?s get the initial conditions when I_app = 0 (i.e. no applied current) 
%and the model runs for 30 msec:

I_app = 0;

%Run solver
[t,x] = ode15s(@BME301_TermProject_Reparameterized_HH_ode_function,[0 5],[0 0 0 0],[],I_app,n_exponent(n));

%Set the state variables x0 to the values of x at 30 msec:
%Note: These were hard-coded directly from the previous ode solution with I_app = 0


x0(1) = x(end,1);
x0(2) = x(end,2);
x0(3) = x(end,3);
x0(4) = x(end,4);

%Now let?s see what the neuron can do!
%Set the applied voltage to 6.2 mV.
%I_app = 6.2;
I_app = 6.2;

%Run solver
[t,x] = ode15s(@BME301_TermProject_Reparameterized_HH_ode_function,[0 5],x0,[],I_app,n_exponent(n));

%Subtract 60mV from voltage variable and plot it vs time

%DONT NEED THIS???
%x_new = x(:,1)-60;
x_zero_baseline = x(:,1)+60;

%Create new figure

%x(1,1) = x(20,1);
%figure
%t_updated = [0;t];
%x_updated = [x(1,1);x(:,1)]; 

%Plot voltage vs time
%plot(t,x_new)
%plot(t,x(:,1),'o')
plot(t,x_zero_baseline)

%plot(t_updated,x_updated)

%Label plot
xlabel('Time (msec)')
ylabel('Voltage (mV)')
title('Action Potential Plot')

%label = sprintf('Voltage vs Time Plot: n = %d and n = %f', n_exponent(1),n_exponent(2));
%title('Voltage vs Time Plot For Action Potential Simulation')
%title(label)

hold on 

%label1 = sprintf('n-exponent = %d',n_exponent(1));
%label2 = sprintf('n-exponent = %d',n_exponent(2));
%legend(label1,label2)
legend('Original Model', 'Reparameterized Model')


%Use your function to solve for the conductance at each time point given 
%the values of m, n and h that you just solved for.





if n_exponent==6

for i = 1:length(x)
 [f varargout(i,:)] = BME301_TermProject_Reparameterized_HH_ode_function(0,x(i,:),I_app,n_exponent(n));
end

figure


plot(t,varargout(:,4).*10^6,'r')
title('Gating Capacitance vs time')
xlabel('time (ms)')
ylabel('Gating capacitance (uF/cm^2)')
end
%{

%Plot the conductances with this code:
%Create new figure
figure
%Plot conductance vs time for potassium and sodium:
plot(t,varargout(:,1),'k','linewidth',1.5);
hold on
plot(t,varargout(:,2),'k--','linewidth',1.5);
%Label plot and format display
title('Action Potential: Conductance vs Time')
set(gca,'Fontsize',18);
hold off; box on; axis([0 30 0 45]);
xlabel('$t$ (msec)','interpreter','latex','fontsize',20);
ylabel('Conductance (mS$\cdot$cm$^{-2}$)','interpreter','latex','fontsize',20);
text(16.4,21,'$g_{Na}$','interpreter','latex','fontsize',20);
text(19,8,'$g_{K}$','interpreter','latex','fontsize',20);

%}

%}

end
