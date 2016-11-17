function pure_ev_mse
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

v = csvread('highways.csv',0,1); %kmph
v = v*0.277778; % m/s
s = sum(v)/1000; %km
s = s/1.60934; % miles

[mass,targetA,~,targetC,mpge] = getPureEVData(1)
mass = mass/2.2;

N = length(targetC);

observedValue = 33.7*s./mpge; %kWh

d = zeros(N,1);

airDensity = 1.2;
numberOfWheels = 4;

T = length(v);

a = zeros(T,1);

for i = 1:T-1
    a(i) = v(i+1)-v(i);
end

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on');%,'DerivativeCheck','on');

[x,fval] = fminunc(@squareError,[3,-5,0.5,10],options)

error = 33.7*s./d;
plot([1:N],[mpge,mpge+error])

avError = sum(abs(error))/N



%fminunc(@likelihood,[0.09,-10.7522,-3.7906,6.85,4.27],options)


    function [f,g] = squareError(x)
        
        sigmoid = @(a) (1+exp(-a))^-1;
        dSigmoid = @(a) exp(-a)*(1+exp(-a))^-2;
        % {
        k1 = exp(x(1));
        k2 = exp(x(2));
        k3 = sigmoid(x(3));
        efficiency = sigmoid(x(4));
        %}
        %{
        k1 = x(1);
        k2 = x(2);
        k3 = x(3);
        efficiency = x(4);
        %}
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

            predictedValue = transpose(efficiencyVector)*P*2.77778*10^-7;
            d(j) = predictedValue-observedValue(j);
            
            % gradients
            dP = efficiencyVector*d(j);

            dk1 = dk1 + 0.5*airDensity*targetC(j)*transpose(dP)*(v.^3);
            dk2 = dk2 + targetA(j)*mass(j)*9.81*transpose(dP)*v;
            dk3 = dk3 + 0.5*numberOfWheels*mass(j)*transpose(dP)*(a.*v);
            deff = deff + d(j)*transpose(P)*efficiencyDerivatives;

        end
        
        f = 0.5*transpose(d)*d;
        % {
        dx1 = k1*dk1*2.77778*10^-7;
        dx2 = k2*dk2*2.77778*10^-7;
        dx3 = dSigmoid(x(3))*dk3*2.77778*10^-7;
        dx4 = dSigmoid(x(4))*deff*2.77778*10^-7;
        %}
        %{
        dx1 = dk1*2.77778*10^-7;
        dx2 = dk2*2.77778*10^-7;
        dx3 = dk3*2.77778*10^-7;
        dx4 = deff*2.77778*10^-7;
        %}
        g = [dx1;dx2;dx3;dx4];
        
    end

end
