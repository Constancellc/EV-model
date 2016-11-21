function [E,MPGe] = predict(mass,Ta,Tb,Tc,efficiency,cycle,slope)
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
T = length(v);                     % length - s 
theta = zero(T);
%theta = cvsread(slope)

% Approximate acceleration
a = zeros(T,1);
for i = 1:T-1
    a(i) = v(i+1)-v(i);
end


efficiencyVector = zeros(T,1);
for j = 1:T
    if a(j) >= 0
        efficiencyVector(j) = 1/efficiency;
    else
        efficiencyVector(j) = efficiency;
    end
end

% in N
F = mass*(a+9.81*sin(theta)) + 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);
P = F.*v;
E = transpose(efficiencyVector)*P;
MPGe = 75384669*s/E;

end 

