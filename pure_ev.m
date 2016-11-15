function pure_ev

%[highway,udds] = getData;

v0 = csvread('UDDS.csv',0,1); %kmph
v0 = v0*0.277778; % m/s
s = sum(v0)/1000; %km

targetC = 0.01;
targetA = 40;
mass = 1000;
mpge = 100;

kmpge = mpge*1.60934;

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


options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');
fminunc(@likelihood,[log(20),log(0.2),log(0.2),5,1],options)

    function [f,g] = likelihood(x)
        k1 = exp(x(1));
        k2 = exp(x(2));
        k3 = exp(x(3));
        var = x(5);
        
        CdA = k1*targetC;
        rollingResistance = k2*targetA;
        massWheels = k3*mass;
        
        % constrain efficiency between 0 and 1
        efficiency = 1/(1+exp(-x(4)));
        
        efficiencyVector = zeros(T,1);
        efficiencyDerivatives = zeros(T,1);
        
        for j = 1:T
            if a(j) >= 0
                efficiencyVector(j) = 1/efficiency;
                efficiencyDerivatives(j) = 1/efficiency^2;
            else
                efficiencyVector(j) = efficiency;
                efficiencyDerivatives(j) = 1;
            end
        end

        P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass+0.5*massWheels*numberOfWheels)+...
            rollingResistance*mass*9.81*v;
        
        predictedValue = transpose(efficiencyVector)*P;
        d = predictedValue-observedValue;
        f = 0.5*log(2*var^2*pi)+0.5*transpose(d)*d/var^2;
        
        % gradients
        dk1 = 0.5*airDensity*targetC*transpose(efficiencyVector)*(v.^3);
        dk2 = 0.5*targetA*mass*9.81*transpose(efficiencyVector)*v;
        dk3 = 0.5*numberOfWheels*mass*transpose(efficiencyVector)*(a.*v);
        
        deff = transpose(efficiencyDerivatives)*P;
        dvar = 1/var*(1-(1/var^2)*(predictedValue-observedValue));
        
        dx1 = k1*dk1;
        dx2 = k2*dk2;
        dx3 = k3*dk3;
        
        dx4 = exp(-x(4))*((1+exp(-x(4)))^-2)*deff;
        
        g = [dx1;dx2;dx3;dx4;dvar];
        
    end


end
