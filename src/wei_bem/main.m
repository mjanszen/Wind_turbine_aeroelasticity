% this code is for the course: wind turbine aeroelasticity
% Delft university of technology

clc
close all;
clear all;

%% load key parameters
load 'STATE'
v0=WindSpeeds;    % read wind speed array
omega=RtSpeeds*2*pi/60;  %rotation angular velocity
pitch=PitchAngles;   % pitch angle: collective pitch wind turbine

%% loop: wind speed from 3m/s to 25m/s
for i=1:length(v0)
    % call BEM, outputs: Radius, Loads, and power outputs
    [Rx,FN,FT,P(i)]=BEM(v0(i),omega(i),pitch(i));
    
    % plot FN and FT
    figure(1)
    plot(Rx,FN,'r-o');
    hold on;
    plot(Rx,FT,'b-o');
    hold on
    grid on
    xlabel('Radius(m)');
    ylabel('Loads(N)');
    legend('Fn','Ft');
end

%% plot power curve regarding wind speed
figure(2)
plot(v0,P,'b-o','linewidth',1.5);
xlabel('Wind speed(m/s)');
ylabel('Power(W)');