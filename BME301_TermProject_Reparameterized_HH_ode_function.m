%BME 301 TERM PROJECT MATLAB CODE
%Reparameterized Hodgkin-Huxley Model

%THIS IS THE FUNCTION FILE THAT GETS CALLED FROM THE MASTER FILE

function [f,varargout] = BME301_TermProject_Reparameterized_HH_ode_function(t,x,I_app,n_exponent)
%BME 301 Term Project - Model an entire action potential using
%Reparameterized HH model

%t is the dummy variable
%x is the state variable (holding v, m, n and h)
%I_app is the applied current (total current through membrane)

%Given:
%Nernst/Resting potentials, conductivities, and capacitance:

%V_Na = 115; Use this value if +65 is not already in model equations
V_Na = 50; %mV

%V_K = -12; Use this value if +65 is not already in model equations
V_K = -77; %mV

%V_L = 10.6; Use this value if +65 is not already in model equations
V_L = -69; %mV

g_Na = 120; %mS/cm^2

g_K = 36; %mS/cm^2

g_L = 0.3; %mS/cm^2

%Capacitance 
C_m = 8.8e-7; % 1 Micro_F = 1e-6 F
C_max = 1.3e-7;


%State Variables
v = x(1);
m = x(2);
n = x(3);
h = x(4);

%Constants for calculating alphas and betas: 

a_n1 = .01;
a_n2 = 10;
b_n1 = .125;
b_n2 = 80;
a_m1 = .1;
a_m2 = 25;
b_m1 = 4;
a_h1 = .07;

if n_exponent == 4
    b_h1 = 1;
    b_h2 = 30;
else
    b_h1 = 1.8;
    b_h2 = 49;
end 
    

%General Equations for forward and backward rate constants
a_n = (a_n1*(-(v+65)+a_n2))/(exp((-(v+65)+a_n2)/10)-1);
b_n = b_n1*exp(-(v+65)/b_n2);

a_m = (a_m1*(-(v+65)+a_m2))/(exp((-(v+65)+a_m2)/10)-1);
b_m = b_m1*exp(-(v+65)/18);

a_h = a_h1*exp(-(v+65)/20);
b_h = (b_h1)/(exp((-(v+65)+b_h2)/10)+1);

% Computing currents: For these, we can look at the HH equation including all the 
% currents, and we can separate them into three different terms (these
% terms are for sodium, potassium and leakage currents).
I_K = g_K.*(n.^(n_exponent)).*(v-V_K);         %n_exp usually 4
I_Na = g_Na.*(m.^3).*(h).*(v-V_Na);
I_L = g_L.*(v-V_L);

%Computing derivatives: 

if n_exponent==4 %Use these equations for original model
    f(1,1) = (1/C_m).*(I_app-(I_K)-(I_Na)-(I_L)); %dV/dt
    f(2,1) = a_m*(1-m)-(b_m*m); %dm/dt
    f(3,1) = (a_n*(1-n))-(b_n*n); %dn/dt
    f(4,1) = a_h*(1-h)-(b_h*h); %dh/dt
else %Use these equations for reparameterized model
    f(1,1) = (1/(C_max.*(1-m))).*(1/C_m).*(I_app-(I_K)-(I_Na)-(I_L)+(C_max*v*(a_m*(1-m)-(b_m*m)))); %dV/dt
    f(2,1) = a_m*(1-m)-(b_m*m); %dm/dt
    f(3,1) = (a_n*(1-n))-(b_n*n); %dn/dt
    f(4,1) = a_h*(1-h)-(b_h*h); %dh/dt 
end
    

%Output the conductivities and gating current/capacitance (if applicable) to varargout:
%We have the following...
%[Sodium conductance in terms of m and h, Potassium conductance in terms of n, g_L]

%varargout{1} = [g_Na.*m.^3.*h,g_K.*n.^(n_exponent), g_L];

%Gating current
I_gating = -v.*C_max.*(a_m*(1-m)-(b_m*m))+C_max.*(1-m).*((1/(C_max.*(1-m))).*(1/C_m).*(I_app-(I_K)-(I_Na)-(I_L)+(C_max*v*(a_m*(1-m)-(b_m*m)))));
%Gating capacitance
Gating_Capacitance_forplotting = C_max.*(1-m);

if n_exponent==4
varargout{1} = [g_Na.*m.^3.*h,g_K.*n.^(n_exponent), g_L];
else
varargout{1} = [g_Na.*m.^3.*h,g_K.*n.^(n_exponent), g_L,Gating_Capacitance_forplotting,I_gating];
end

end




