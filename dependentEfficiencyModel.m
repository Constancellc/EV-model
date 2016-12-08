function dependentEfficiencyModel()

v_h = csvread('highways.csv',0,1)*0.277778; % m/s
s_h = sum(v_h)/1000; %km
T_h = length(v_h); a_h = zeros(T_h,1); 
vh_av = s_h*1000/T_h;

for j = 1:T_h-1
    a_h(j) = v_h(j+1)-v_h(j);
end

v_u = csvread('udds.csv',0,1)*0.277778; % m/s
s_u = sum(v_u)/1000; %km
T_u = length(v_u); a_u = zeros(T_u,1); 
vu_av = s_u*1000/T_u;

for j = 1:T_u-1
    a_u(j) = v_u(j+1)-v_u(j);
end

[mass,targetA,targetB,targetC,hwys,udds] = getPureEVData;
data = csvread('../EPA-Code/ev_data.csv',0,2);
horsepower = data(:,2);

% converting from MPGe to energy consumption
hwys = 75384669*s_h./hwys;
udds = 75384669*s_u./udds;

N = length(targetC);

% in N
force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

P_h = zeros(N,T_h); P_u = zeros(N,T_u);
Th_rms = zeros(N,1); Tu_rms = zeros(N,1); 

for j = 1:N
    Ta = targetA(j); Tb = targetB(j); Tc = targetC(j);

    F_h = force(v_h,Ta,Tb,Tc)+mass(j)*a_h;
    F_u = force(v_u,Ta,Tb,Tc)+mass(j)*a_u;

    P_h(j,:) = F_h.*v_h; P_u(j,:) = F_u.*v_u;
    
    Th_rms(j) = sum(F_h.^2)/N; Tu_rms(j) = sum(F_u.^2)/N;
end
% {
options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on');%,'DerivativeCheck','on');
%options.StepTolerance = 10e-10;
%options.MaxFunctionEvaluations = 300;
%options.MaxIterations = 400;
%options.FunctionTolerance = 1e-6;
x = fminunc(@likelihood,[1;10^-4],options)
%x=fmincon(@likelihood,[1;1;1])
exp(x(1))
exp(x(2))
%exp(x(3))
%x = fmincon(@likelihood,[1;(1/vu_av);(1/sum(P_u(1,:)))])


% {
% I really want to covert back to MPGe
hwys = 75384669*s_h./hwys;
udds = 75384669*s_u./udds;
% {
h = 75384669*s_h./predictions(:,1);
%h_l = 75384669*s_h./(predictions(:,1)+error);
%h_u = 75384669*s_h./(predictions(:,1)-error);

u = 75384669*s_u./predictions(:,2);
%u_l = 75384669*s_u./(predictions(:,2)+error);
%u_u = 75384669*s_u./(predictions(:,2)-error);
%}

% {
errors_h = h-hwys; errors_u = u-udds;

h_u = zeros(N,1); h_l = zeros(N,1); u_u = zeros(N,1); u_l = zeros(N,1);

for j = 1:N
    if errors_h(j) <= 0
        h_u(j) = errors_h(j);
    else
        h_l(j) = errors_h(j);
    end
    
    if errors_u(j) <= 0
        u_u(j) = errors_u(j);
    else
        u_l(j) = errors_u(j);
    end
end
av1 = sum(abs(errors_h))/N
av2 = sum(abs(errors_u))/N
%}
% {
figure(1)
subplot(2,1,1)
errorbar([1:N],h,h_l,h_u)
%plot([1:N],h,'o')
hold on
plot([1:N],hwys,'x')
title('Predicted vs. Observed MPGe on Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')

subplot(2,1,2)
errorbar([1:N],u,u_l,u_u)
%plot([1:N],u,'o')
hold on
plot([1:N],udds,'x')
title('Predicted vs. Observed MPGe on Urban Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')
%}

%{
% plot varied efficiency
x1 = 0:0.01:1;
for k = 1:length(x1)
    y1(k) = 1/likelihood([x1(k);x(2)]);
end

% and then vary variance
x2 = 0:0.1:40;
for k = 1:length(x2)
    y2(k) = 1/likelihood([efficiency;x2(k)]);
end

figure(2)
subplot(2,1,1)
plot(x1,y1)
subplot(2,1,2)
plot(x2,y2)

BELOW ARE THE PLOTS FOR ENERGY USED RATHER THAN MPGe
----------------------------------------------------

figure(1)
subplot(2,1,1)
errorbar([1:N],predictions(:,1),error);
hold on
plot([1:N],hwys,'x')
subplot(2,1,2)
errorbar([1:N],predictions(:,2),error);
hold on
plot([1:N],udds,'x')
%}


function [f,g] = likelihood(x)
    k1 = x(1); k2 = x(2);
    
    d_h = zeros(N,1); d_u = zeros(N,1);
    
    g = zeros(2,1);
    
    for j = 1:N
        eff = k1*(1 + k2*horsepower(j));
        
        eff_h = zeros(T_h,1); eff_u = zeros(T_u,1); 
        deff_h = zeros(T_h,1); deff_u = zeros(T_u,1); 
        
        for i = 1:T_u-1
            if a_u(i) < 0
                eff_u(i) = eff; deff_u(i) = 1;
            else
                eff_u(i) = 1/eff; deff_u(i) = -1/eff^2;
            end
        end

        for i = 1:T_h-1
            if a_h(i) < 0
                eff_h(i) = eff; deff_h(i) = 1;
            else
                eff_h(i) = 1/eff; deff_h(i) = -1/eff^2;
            end
        end
    
        predictions(j,:) = [P_h(j,:)*eff_h,P_u(j,:)*eff_u];
        d_h(j) = (P_h(j,:)*eff_h)-hwys(j);
        d_u(j) = (P_u(j,:)*eff_u)-udds(j);
        
        g(1) = g(1) + 2*(d_h(j)*P_h(j,:)*deff_h + ...
            d_u(j)*P_u(j,:)*deff_u)*(1+k2*horsepower(j));
        g(2) = g(2) + k1*horsepower(j)*(2*d_h(j)*P_h(j,:)*deff_h + ...
            2*d_u(j)*P_u(j,:)*deff_u);

        
    end
    
    %g(1) = g(1)*k1; g(2) = g(2)*k2; g(3) = g(3)*k3; 
    d = [d_h;d_u];

    f = transpose(d)*d;

end
end