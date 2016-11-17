function pure_ev
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

v = csvread('highways.csv',0,1); %kmph
v = v*0.277778; % m/s
s = sum(v)/1000; %km

T = length(v);
a = zeros(T,1);

for i = 1:T-1
    a(i) = v(i+1)-v(i);
end

[mass,targetA,targetB,targetC,mpge] = getPureEVData(1);

N = length(targetC);

efficiency = 0.7;

% in N
force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707));

efficiencyVector = zeros(T,1);
for i = 1:T
    if a(i) <= 0
        efficiencyVector(i) = 1/efficiency;
    else
        efficiencyVector(i) = efficiency;
    end
end

prediction = zeros(N,1);

for i = 1:N
    f = mass(i)*a;
    
    Ta = targetA(i);
    Tb = targetB(i);
    Tc = targetC(i);
    
    f = f + force(v,Ta,Tb,Tc);
    
    P = f.*v;
    E = transpose(efficiencyVector)*P; %J
    prediction(i) = 9.36112*10^6*s/E;
end

plot([1:N],[mpge,prediction])


end
