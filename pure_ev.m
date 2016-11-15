function pure_ev
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

v0 = csvread('UDDS.csv',0,1); %kmph
v0 = v0*0.277778; % m/s
s = sum(v0)/1000; %km

[~,udds] = getData;

targetC = udds(:,4);
targetA = udds(:,3);
mass = udds(:,2);
mpge = udds(:,1);

N = length(targetC);

kmpge = mpge*1.60934;

observedValue = 1000*33.7*s/kmpge; %W

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

iSigmoid = 

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on');%,'DerivativeCheck','on');
fminunc(@likelihood,[log(20),log(0.2),log(0.2),5,1],options)



    function [f,g] = likelihood(x)
        
        sigmoid = @(a) (1+exp(-a))^-1;
        dSigmoid = @(a) exp(-a)*(1+exp(-a))^-2;
        
        k1 = exp(x(1));
        k2 = exp(x(2));
        k3 = sigmoid(x(3));
        efficiency = sigmoid(x(4));
        var = x(5);
        
        
        efficiencyVector = zeros(T,1);
        efficiencyDerivatives = zeros(T,1);
        
        for j = 1:T
            if a(j) >= 0
                efficiencyVector(j) = 1/efficiency;
                efficiencyDerivatives(j) = -1/efficiency^2;
            else
                efficiencyVector(j) = efficiency;
                efficiencyDerivatives(j) = 1;
            end
        end
        
        f = 0;
        dk1 = 0;
        dk2 = 0;
        dk3 = 0;
        deff = 0;
        dvar = 0;
        
        % Note: a for loop is a terrible way to do this
        for j = 1:N
            
            CdA = k1*targetC(j);
            rollingResistance = k2*targetA(j);
            massWheels = k3*mass(j);

            P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass(j)+0.5*massWheels*numberOfWheels)+...
                rollingResistance*mass(j)*9.81*v;

            predictedValue = transpose(efficiencyVector)*P;
            d = predictedValue-observedValue(j);
        
            f = 0.5*log(2*var^2*pi)+0.5*d^2/(var^2);

            % gradients
            dP = (1/var^2)*efficiencyVector*d;

            dk1 = dk1 + 0.5*airDensity*targetC(j)*transpose(dP)*(v.^3);
            dk2 = dk2 + targetA(j)*mass(j)*9.81*transpose(dP)*v;
            dk3 = dk3 + 0.5*numberOfWheels*mass(j)*transpose(dP)*(a.*v);

            deff = deff + (1/var^2)*d*transpose(P)*efficiencyDerivatives;
            dvar = dvar + 1/var*(1-(d/var)^2);
        end
        
        dx1 = k1*dk1;
        dx2 = k2*dk2;
        dx3 = dSigmoid(x(3))*dk3;
        dx4 = dSigmoid(x(4))*deff;
        
        g = [dx1;dx2;dx3;dx4;dvar];
        
    end


end
