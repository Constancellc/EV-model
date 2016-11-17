function pure_ev_mse
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

% GRADIENTS CURRENTLY NOT WORKING

v = csvread('highways.csv',0,1); %kmph
v = v*0.277778; % m/s
s = sum(v)/1000; %km
s = s/1.60934; % miles

% conversion from J to MPGe constant
K = 9.36112*10^6;

[mass,targetA,~,targetC,mpge] = getPureEVData(1);

N = length(targetC);

observedValue = mpge;%33.7*s./mpge; %kWh
d = zeros(N,1);

airDensity = 1.2;
numberOfWheels = 4;

T = length(v);

a = zeros(T,1);

for i = 1:T-1
    a(i) = v(i+1)-v(i);
end

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');

startPoints = 1;

latinHypercube = lhsdesign(startPoints,4);

% k1 between 1 and 100 => logk1 between 0 and 2
latinHypercube(:,1) = latinHypercube(:,1)*2;
% k2 between 10^-4 and 10^-2 => logk2 between -4 and -2
latinHypercube(:,2) = latinHypercube(:,2)*(-2)-2*ones(startPoints,1);
% k3 between 0.005 and 0.2 => x3 between -6 and -1
%latinHypercube(:,3) = latinHypercube(:,3)*(-5)-ones(startPoints,1);
% efficiency between 0.85 and 0.999 => x4 between 1 and 7
latinHypercube(:,3) = latinHypercube(:,3)*6+ones(startPoints,1);
% variance between 0 and 20 => x5 between 0 and 5
latinHypercube(:,4) = latinHypercube(:,4)*2;

best = 10^5;

for i = 1:startPoints
    [x,fval] = fminunc(@likelihood,[3.5,8,6.85,1.27],options);
    error(:,i) = d;
    sigma(i) = x(4);
    if fval <= best
        best = fval;
        parameters = x;
    end
end

exp(parameters(1))
exp(parameters(2))
sigmoid(parameters(3))
parameters(4)

obs = observedValue;
%error = 33.7*s./error;

plot([1:N],[obs,obs+error])
%{
figure
for i = 1:startPoints
    y = obs+error(:,i);
    l = y-2*sigma(i)*ones(N,1);
    u = y+2*sigma(i)*ones(N,1);
    
    
    subplot(4,3,i)
    plot([1:N],[obs])
    %hold on
    %errorbar([1:N],y,l,u)
end
%}
avError = sum(abs(error))/N
invSigmoid = @(a) log(a)-log(1-a);

    function [f,g] = likelihood(x)
        
        sigmoid = @(a) (1+exp(-a))^-1;
        dSigmoid = @(a) exp(-a)*(1+exp(-a))^-2;
        
        k1 = exp(x(1));
        k2 = exp(x(2));
        %k3 = sigmoid(x(3));
        efficiency = sigmoid(x(3));
        var = x(4);
             
        massWheels = 30;
        
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
        
        f = ((k1-33.5)/12)^2+((k2-3200)/4000)^2+((efficiency-0.707)/0.1)^2;
        f = 0.5*f;
        
        dk1 = (k1-33.5)/(12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ;
        dk2 = (k2-3200)/4000;
        %dk3 = 0;
        deff = (efficiency-0.707)/0.1;
        dvar = 0;
        
        % Note: a for loop is a terrible way to do this
        for j = 1:N
            
            CdA = k1*targetC(j);
            rollingResistance = k2*targetA(j);
            %massWheels = 30;%k3*mass(j);

            P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass(j)+0.5*massWheels*numberOfWheels)+...
                rollingResistance*mass(j)*9.81*v;

            E = transpose(efficiencyVector)*P;
            predictedValue = K*s/E;
            d(j) = predictedValue-observedValue(j)
        
            f = f + 0.5*log(2*var^2*pi)+0.5*d(j)^2/(var^2);

            % gradients
            dP = -(1/var^2)*d(j)*K*s/(E^2);
            %dP = (1/var^2)*efficiencyVector*d(j);

            dk1 = dk1 + 0.5*airDensity*targetC(j)*dP*transpose(efficiencyVector)*(v.^3);
            dk2 = dk2 + targetA(j)*mass(j)*9.81*dP*transpose(efficiencyVector)*v;
            %dk3 = dk3 + 0.5*numberOfWheels*mass(j)*transpose(dP)*(a.*v);

            deff = deff + dP*transpose(P)*efficiencyDerivatives;
            dvar = dvar + (1/var)*(1-(d(j)/var)^2);
        end
        
        dx1 = k1*dk1;
        dx2 = k2*dk2;
        %dx3 = dSigmoid(x(3))*dk3*2.77778*10^-7;
        dx3 = dSigmoid(x(3))*deff;
        
        g = [dx1;dx2;dx3;dvar];
        
    end
%}


end
