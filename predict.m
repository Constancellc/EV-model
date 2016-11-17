function [E,MPGe] = predict(mass,targetA,targetC,x,cycle)
%{
inputs: mass - vehicle mass (kg)
        targetA/C - coastdown coefficients
        x - a vector containing model parameters
        cycle - a string containing the name of the drive cycle csv file

Output: Predicted energy consumption of vehicle over drive cycle
%}

% Load drive cycle
v = csvread(cycle,0,1)*0.277778;   % m/s
s = sum(v)/1000;                   % range - km
miles = s*0.621371192;                 % convert to miles
T = length(v);                     % length - s 
% Approximate acceleration
a = zeros(T,1);
for i = 1:T-1
    a(i) = v(i+1)-v(i);
end

k1 = x(1);
k2 = x(2);
%k3 = x(3);
efficiency = 0.99;%x(3);

efficiencyVector = zeros(T,1);
for j = 1:T
    if a(j) >= 0
        efficiencyVector(j) = 1/efficiency;
    else
        efficiencyVector(j) = efficiency;
    end
end

airDensity = 1.2;
numberOfWheels = 4;

CdA = 0.2;%k1*targetC;
rollingResistance = 0.014;%k2*targetA;
massWheels = 30;%k3*mass;

P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass+0.5*massWheels*numberOfWheels)+...
    rollingResistance*mass*9.81*v;

E = transpose(efficiencyVector)*P; %J
E = E*2.77778*10^-7; %kWh

gallonsEquivalent = E/33.7;
MPGe = miles/gallonsEquivalent;
% 

