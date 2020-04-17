% BME301_conduction_vel_numerical_method.m
% BME 301 TERM PROJECT MATLAB CODE
%
% This is the code to perform a numerical analysis of the voltage along the
% membrane over both time and space using Euler's method
%
% Supplementary code needed: F.m, H.m, M.m, and N.m
% 
% Adapted from: "Mathematical Modeling of Action Potential with
% Transmission Equations and Hodgkin-Huxley Model", BENG 221 Problem
% Solving Report
% https://isn.ucsd.edu/courses/beng221/problems/2010/project6.pdf
%
% To run the code, a .mat file is needed that contains a 1x4 vector or a 
% structure with one field containing the initial values [n, m, h, I_app]
% n: Potassium n gate
% m: Sodium m gate
% h: Sodium h gate
% I_app: Applied voltage


%%% Set up and initialization of variables
clc;
close all;
clear all;
dt=0.0001; % Set time step, units of ms
dx=1; % Set distance step, units of cm
s=dt/dx; % Conversion between dt and dx
t=0:dt:10; % Total time of 10 ms
x=0:dx:10; % Total distance of 10 cm
n=zeros(length(x),length(t));  % Potassium n gates
m=zeros(length(x),length(t));  % Sodium m gates
h=zeros(length(x),length(t));  % Sodium h gates
v=zeros(length(x),length(t)); % Membrane voltage, units of mV
psi=zeros(length(x),length(t)); % dV/dt, units of mV/ms
phi=zeros(length(x),length(t)); %dV/dx, units of mV/cm


%%% Initial conditions

% Prompt the user to choose a .mat file
[filename, pathname] = uigetfile('*.mat', 'Please select the initial conditions');

% Check to make sure the user chose a file
if pathname == 0
    % If not, display an error message, and stop the callback function
    fprintf('Error: A file was not selected\n');
    return
end

% Create the full file name, and load it
pathfile = fullfile(pathname, filename);
initconds = load(pathfile);

% Check to make sure that the mat file is formated correctly
if length(initconds) == 4 && all(isnumeric(initconds))
    % This is a vector with four numerical values
elseif isstruct(initconds)
    % Check if the structure has one field, and see if that field has four
    % values. If not, alert the user and stop the script
    field = fieldnames(initconds);
    if length(field) ~= 1
        fprintf('Error: The MAT file selected is not compatible\n');
        help BME301_conduction_vel_numerical_method
        return
    elseif length(initconds.(field{1})) ~=4
        fprintf('Error: The MAT file selected is not compatible\n');
        help BME301_conduction_vel_numerical_method
        return
    end
    % Save the structure as the initial conditions
    initconds = initconds.(field{1});  
else
    % Otherwise, the file is not correct for the script to work
    % Alert the user, gove them help information, and stop the script
    fprintf('Error: The MAT file selected is not compatible\n');
    help BME301_conduction_vel_numerical_method
    return
end

% Initialize n, m, h, and v
n(:,1) = initconds(1);
m(:,1) = initconds(2);
h(:,1) = initconds(3);
v(1,:) = -initconds(4);

%{
n(:,1)=0.1556;
m(:,1)=0.0255;
h(:,1)=0.3357;
v(1,:)=-6.2;
%}

% Set the initial voltages based on the change in distance
for i=1:1:length(x)-1
    phi(i,1)=1.5; 
    v(i+1,1)=v(i,1)+dx*phi(i,1);
end


%%% Numerical Approximation

% Set the beginning of the percentage, and that the numerical approximation
% process is beginning
userper = 0;
fprintf('Performing numerical approximation...\n');

for  j=1:1:length(t)-1
    % Check if index is a new multiple of 10%
    if floor(j*10/(length(t)-1)) > userper
        % If so, alert the user of the progress, and update the percent
        userper = userper + 1;
        fprintf('%d%% done...\n', userper*10);
    end
    for i=1:1:length(x)-1
        % Update all variables based on Euler's method
        n(i,j+1)=n(i,j) + dt*N(n(i,j),v(i,j));
        m(i,j+1)=m(i,j) + dt*M(m(i,j),v(i,j));
        h(i,j+1)=h(i,j) + dt*H(h(i,j),v(i,j));
        v(i,j+1)=v(i,j) + dt*psi(i,j);
        phi(i+1,j)=phi(i,j) + s*(psi(i+1,j)-psi(i,j));
        psi(i,j+1)=psi(i,j) + s*(phi(i+1,j)-phi(i,j)) - ...
            dt*F(psi(i,j),v(i,j),n(i,j),m(i,j),h(i,j));
    end   
end

%%% Plotting and Conduction Velocity

% Plot the action potentials along the axon over time
v = -v;
figure(1) 
plot(t,v(1,:),t,v(2,:),t,v(3,:),t,v(4,:),t,v(5,:),t,v(6,:),t,v(7,:),...
    t,v(8,:), t,v(9,:),t,v(10,:))
legend('1 cm','2 cm','3 cm','4 cm','5 cm','6 cm','7 cm','8 cm',...
    '9 cm','10cm','Location','NW')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potentials along Neuron')

% Find the conduction velocity by considering the 6th and 8th positions

% Calculate the index of each peak
[~, index] = max(v');
index(end)=[];

% Find the time of each peak
peaktimes = t(index);
distances = x;
distances(1) = [];

% Calculate the conduction velocity 
Vc = 20 / (peaktimes(distances==8)-peaktimes(distances==6));
fprintf('The conduction velocity is %f m/s\n', Vc);
