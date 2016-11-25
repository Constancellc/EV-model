function test
% evidence calculation for one test vehicle 
v_h = csvread('highways.csv',0,1)*0.277778; % m/s
s_h = sum(v_h)/1000; %km
T_h = length(v_h); a_h = zeros(T_h,1); 

for i = 1:T_h-1
    a_h(i) = v_h(i+1)-v_h(i);
end

[mass,targetA,targetB,targetC,hwys,udds] = getPureEVData;

% converting from MPGe to energy consumption
hwys = 75384669*s_h./hwys;

force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

F = force(v_h,targetA(1),targetB(1),targetC(1)) + mass(1)*a_h;
P = F.*v_h;

pd = makedist('Gamma','a',0.5,'b',1);

tot = 0;
No = 100000;

for j = 1:No
    e = rand;
    sd = random(pd);
    tot = tot + N(e,sd,hwys(1));
end

p = tot/No

function f = N(e,sd,y)
%e = 0.7263;
    eff = zeros(1,T_h);
    for i = 1:T_h-1
        if a_h(i) < 0
            eff(i) = e;
        else
            eff(i) = 1/e;
        end
    end
    pred = eff*P;
    
    f = (1/sqrt(2*pi*sd^2))*exp(-0.5*(y-pred)^2/(sd^2));
    
end

end