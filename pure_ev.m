function efficiency = pure_ev
% First let's try building a model based on only one drive cycle
%[highway,udds] = getData;

v_h = csvread('highways.csv',0,1)*0.277778; % m/s
s_h = sum(v_h)/1000; %km
T_h = length(v_h); a_h = zeros(T_h,1);

for i = 1:T_h-1
    a_h(i) = v_h(i+1)-v_h(i);
end

v_u = csvread('udds.csv',0,1)*0.277778; % m/s
s_u = sum(v_u)/1000; %km
T_u = length(v_u); a_u = zeros(T_u,1);

for i = 1:T_u-1
    a_u(i) = v_u(i+1)-v_u(i);
end

[mass,targetA,targetB,targetC,hwys,udds] = getPureEVData;

N = length(targetC);

% in N
force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

efficiency = zeros(N,1);
prediction = zeros(N,2);

% First learn efficiency from the highways drive cycle
for i = 1:N
    efficiency(i) = fmincon(@squareError,0.65,[1;-1],[1;0]);
end

% Then test the results on the UDDS cycle

%{
for i = 1:N
    efficiencyVector = zeros(T_h,1);
   
    for k = 1:T_h
        % During braking the power required is less
        if a_h(k) < 0
            efficiencyVector(k) = efficiency(i);
        % While during acceleration more power is required to account for
        % losses
        else
            efficiencyVector(k) = 1/efficiency(i);
        end
    end

    f = mass(i)*a_h;

    Ta = targetA(i);
    Tb = targetB(i);
    Tc = targetC(i);

    f = f + force(v_h,Ta,Tb,Tc);

    P = f.*v_h;
    E = transpose(efficiencyVector)*P; %J
    prediction2(i) = 75384669*s_h/E;
end
    
%}

figure(1)
subplot(3,1,1)
plot([1:N],[hwys,prediction(:,1)])
xlabel('Vehicle No.')
ylabel('predicted vs. actual MPGe on highways')
subplot(3,1,2)
plot([1:N],[udds,prediction(:,2)])
xlabel('Vehicle No.')
ylabel('predicted vs. actual MPGe on UDDS')
subplot(3,1,3)
plot([1:N],efficiency)
xlabel('Vehicle No.')
ylabel('efficiency')


    function [f,g] = squareError(efficiency)
        
        efficiencyVectorh = zeros(T_h,1);
        for j = 1:T_h
            if a_h(j) < 0
                efficiencyVectorh(j) = efficiency;
            else
                efficiencyVectorh(j) = 1/efficiency;
            end
        end

        efficiencyVectoru = zeros(T_u,1);
        for j = 1:T_u
            if a_u(j) < 0
                efficiencyVectoru(j) = efficiency;
            else
                efficiencyVectoru(j) = 1/efficiency;
            end
        end

        Ta = targetA(i);
        Tb = targetB(i);
        Tc = targetC(i);

        F_h = mass(i)*a_h + force(v_h,Ta,Tb,Tc);
        P_h = F_h.*v_h;
        E_h = transpose(efficiencyVectorh)*P_h; %J
        prediction(i,1) = 75384669*s_h/E_h;
        error_h = prediction(i,1) - hwys(i);

        F_u = mass(i)*a_u + force(v_u,Ta,Tb,Tc);
        P_u = F_u.*v_u;
        E_u = transpose(efficiencyVectoru)*P_u; %J
        prediction(i,2) = 75384669*s_u/E_u;
        error_u = prediction(i,2) - udds(i);

        f = 0.5*error_u^2 + 0.5*error_h^2;
        
    end

        
%{        
efficiencyVector = zeros(T,1);
for i = 1:T
    % During braking the power required is less
    if a(i) < 0
        efficiencyVector(i) = efficiency;
    % While during acceleration more power is required to account for
    % losses
    else
        efficiencyVector(i) = 1/efficiency;
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
    prediction(i) = 75384669*s/E;
end


%}



end
