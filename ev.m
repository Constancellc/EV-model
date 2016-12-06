% Ok, let's imagine we've chosen the model, we want a function which
% predicts the energy consumption during a trip

function energy = ev(v,mVehicle,mPassenger,Ta,Tb,Tc)
% inputs: drive cycle, vehicle mass, passenger mass, coast-down
% coefficients, efficiency?

% v - drive cycle
% mVehicle - kerb weight
% mPassenger - load

% output: energy expenditure

% For simpicity let's use the constant efficiency model -- but assume that
% as our model runs it can self adjust?

mass = mVehicle + mPassenger;

n = 0.7;

dt = 1;
T = length(v); a = zeros(T,1); eff = ones(1,T);

for i = 1:T-1
    a(i) = (v(i+1)-v(i));
    if a(i) <= 0
        eff(i) = n;
    else
        eff(i) = 1/n;
    end
end

% in N
force = 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

P = force.*v + mass*a;
energy = eff*P*2.77778e-7; % kWh

end