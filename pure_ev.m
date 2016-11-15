function pure_ev
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

v0 = csvread('UDDS.csv',0,1); %kmph
v0 = v0*0.277778; % m/s
s = sum(v0)/1000; %km

[highway,udds] = getData;

targetC = udds(:,4);
targetA = udds(:,3);
mass = udds(:,2);
mpge = udds(:,1);

N = length(targetC);

kmpge = mpge*1.60934;

observedValue = 1000*33.7*s./kmpge; %W
d = zeros(N,1);

airDensity = 1.2;
numberOfWheels = 4;

T = length(v0);

v = zeros(T,1);
a = zeros(T,1);

for i = 1:T-1
    a(i) = v0(i+1)-v0(i);
end


options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on');%,'DerivativeCheck','on');
%fminunc(@squareError,[0,0,-3.7906,6.85],options)
%fminunc(@likelihood,[0.09,-10.7522,-3.7906,6.85,4.27],options)

x =  [4,-2,3.1798,8.4927];

figure(1) 
plot([1:N],d,'x')
ylabel('Training Error')

% Now let's get the testing data

targetC = highway(:,4);
targetA = highway(:,3);
mass = highway(:,2);
mpge = highway(:,1);

N = length(targetC);
kmpge = mpge*1.60934;

vtest = csvread('UDDS.csv',0,1); %kmph
vtest = vtest*0.277778; % m/s
s = sum(vtest)/1000; %km

T = length(vtest);

v = abs(vtest);
a = zeros(T,1);

for i = 1:T-1
    a(i) = vtest(i+1)-vtest(i);
end


observedValue = 1000*33.7*s./kmpge 

E = predict(x);

figure(2)
plot([1:N],[E,observedValue],'x')

    function E = predict(x)
        
        sigmoid = @(a) (1+exp(-a))^-1;
        
        k1 = 10;%exp(x(1));
        k2 = 0.1;%exp(x(2));
        k3 = 0.2;%sigmoid(x(3));
        efficiency = 0.9;%sigmoid(x(4));
        
        efficiencyVector = zeros(T,1);
        for j = 1:T
            if a(j) >= 0
                efficiencyVector(j) = 1/efficiency;
            else
                efficiencyVector(j) = efficiency;
            end
        end
        
        E = zeros(N,1);
        
        for j = 1:N
            
            CdA = k1*targetC(j);
            rollingResistance = k2*targetA(j);
            massWheels = k3*mass(j);

            P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass(j)+0.5*massWheels*numberOfWheels)+...
                rollingResistance*mass(j)*9.81*v;

            E(j) = transpose(efficiencyVector)*P;
                        
        end
    end


    function [f,g] = squareError(x)
        
        sigmoid = @(a) (1+exp(-a))^-1;
        dSigmoid = @(a) exp(-a)*(1+exp(-a))^-2;
        
        k1 = exp(x(1));
        k2 = exp(x(2));
        k3 = sigmoid(x(3));
        efficiency = sigmoid(x(4));
        
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
        
        % Note: a for loop is a terrible way to do this
        for j = 1:N
            
            CdA = k1*targetC(j);
            rollingResistance = k2*targetA(j);
            massWheels = k3*mass(j);

            P = 0.5*airDensity*CdA*v.^3 + a.*v*(mass(j)+0.5*massWheels*numberOfWheels)+...
                rollingResistance*mass(j)*9.81*v;

            predictedValue = transpose(efficiencyVector)*P;
            d(j) = predictedValue-observedValue(j);
            
            % gradients
            dP = efficiencyVector*d(j);

            dk1 = dk1 + 0.5*airDensity*targetC(j)*transpose(dP)*(v.^3);
            dk2 = dk2 + targetA(j)*mass(j)*9.81*transpose(dP)*v;
            dk3 = dk3 + 0.5*numberOfWheels*mass(j)*transpose(dP)*(a.*v);
            deff = deff + d(j)*transpose(P)*efficiencyDerivatives;
            
        end
        
        f = 0.5*transpose(d)*d;
        
        dx1 = k1*dk1;
        dx2 = k2*dk2;
        dx3 = dSigmoid(x(3))*dk3;
        dx4 = dSigmoid(x(4))*deff;
        
        g = [dx1;dx2;dx3;dx4];
        
    end
        
%{


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
            d(j) = predictedValue-observedValue(j);
        
            f = 0.5*log(2*var^2*pi)+0.5*d(j)^2/(var^2);

            % gradients
            dP = (1/var^2)*efficiencyVector*d(j);

            dk1 = dk1 + 0.5*airDensity*targetC(j)*transpose(dP)*(v.^3);
            dk2 = dk2 + targetA(j)*mass(j)*9.81*transpose(dP)*v;
            dk3 = dk3 + 0.5*numberOfWheels*mass(j)*transpose(dP)*(a.*v);

            deff = deff + (1/var^2)*d(j)*transpose(P)*efficiencyDerivatives;
            dvar = dvar + 1/var*(1-(d(j)/var)^2);
        end
        
        dx1 = k1*dk1;
        dx2 = k2*dk2;
        dx3 = dSigmoid(x(3))*dk3;
        dx4 = dSigmoid(x(4))*deff;
        
        g = [dx1;dx2;dx3;dx4;dvar];
        
    end
%}


end
