function pure_ev

targetC = 0.01;
targetA = 40;
mass = 1000;
mpge = 100;

kmpge = mpge*1.60934;

v0 = csvread('UDDS.csv',0,1); %kmph
v0 = v0*0.277778; % m/s
s = sum(v0)/1000; %km

observedValue = 33.7*s/kmpge; %kW

airDensity = 1.2;
numberOfWheels = 4;

% vector of power req.
%P = mew_rr*mass*9.81 + rho_air*CdA/(2*eff)

T = length(v0)-1;

v = zeros(T,1);
a = zeros(T,1);

for i = 1:T
    v(i) = (v0(i)+v0(i+1))/2;
    a(i) = v0(i+1)-v0(i);
end

    function [f,g] = likelihood(x)
        k1 = x(1);
        k2 = x(2);
        k3 = x(3);
        var = x(5);
        
        CdA = k1*targetC;
        rollingResistance = k2*targetA;
        massWheels = k3*mass;
        efficiency = 1/(1+exp(-x(4));
        
        efficiencyVector = zeros(T,1);
        
        for j = 1:T
            if a(j) >= 0
                efficiencyVector(j) = 1/efficiency;
            else
                efficiencyVector(j) = efficiency;
            end
        end

        P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass+0.5*massWheels*numberOfWheels)+...
            rollingResistance*mass*9.81*v;
        
        predictedValue = transpose(efficiencyVector)*P;
        
        f = -0.5*log(2*var^2*pi)-0.5*(predictedValue-observedValue)^2/var^2;
    end





end
