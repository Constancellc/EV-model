function vehicleSpecificEfficiencyModel()

v_h = csvread('highways.csv',0,1)*0.277778; % m/s
s_h = sum(v_h)/1000; %km
T_h = length(v_h); a_h = zeros(T_h,1); 

for j = 1:T_h-1
    a_h(j) = v_h(j+1)-v_h(j);
end

v_u = csvread('udds.csv',0,1)*0.277778; % m/s
s_u = sum(v_u)/1000; %km
T_u = length(v_u); a_u = zeros(T_u,1); 

for j = 1:T_u-1
    a_u(j) = v_u(j+1)-v_u(j);
end

[mass,targetA,targetB,targetC,hwys,udds] = getPureEVData;

% converting from MPGe to energy consumption
hwys = 75384669*s_h./hwys;
udds = 75384669*s_u./udds;

N = length(targetC);

% in N
force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

% constructing hessian and grad x = [eff,var]

P_h = zeros(N,T_h); P_u = zeros(N,T_u);
% eff, deff, d2eff, P, d

for j = 1:N
    Ta = targetA(j); Tb = targetB(j); Tc = targetC(j);

    F_h = force(v_h,Ta,Tb,Tc)+mass(j)*a_h;
    F_u = force(v_u,Ta,Tb,Tc)+mass(j)*a_u;

    P_h(j,:) = F_h.*v_h; P_u(j,:) = F_u.*v_u;
end

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');
[x,fval] = fminunc(@likelihood2,[ones(N,1)],options)

%csvwrite('efficiencies.csv',x(2:end))
%predictions = [P_h*eff_h,P_u*eff_u];
%error = sd*ones(N,1);

% {
% I really want to covert back to MPGe
%hwys = 75384669*s_h./hwys;
%udds = 75384669*s_u./udds;
% {
h = predictions(:,1);
%h = 75384669*s_h./predictions(:,1);

%u = 75384669*s_u./predictions(:,2);
u = predictions(:,2);
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
av1 = sum(abs(errors_h))/N;
av2 = sum(abs(errors_u))/N;
%}
% {
figure(1)
subplot(2,1,1)
b = bar([hwys,h]);
b(1).FaceColor = [0.6 0.9 1];
b(2).FaceColor = [0 0.7 0.7];
%errorbar([1:N],h,h_l,h_u)
%plot([1:N],h,'o')
%hold on
%plot([1:N],hwys,'x')
title('Predicted vs. Observed MPGe on Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')

subplot(2,1,2)
bar([udds,u])
%errorbar([1:N],u,u_l,u_u)
%plot([1:N],u,'o')
%hold on
%plot([1:N],udds,'x')
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
    
    eff = x;
    %k = 0.93;
    
    predictions = zeros(N,2);
    g = zeros(N,1);

    
    for v = 1:N
        eff_h = zeros(T_h,1); eff_u = zeros(T_u,1);
        deff_h = zeros(T_h,1); deff_u = zeros(T_u,1);
        %{
        for i = 1:T_u-1
            if a_u(i) < 0
                eff_u(i) = k*eff(v); deff_u(i) = k;
            else
                eff_u(i) = 1/(k*eff(v)); deff_u(i) = -1/(k*eff(v)^2);
            end
        end
        %}
        for i = 1:T_u-1
            if a_u(i) < 0
                eff_u(i) = eff(v); deff_u(i) = 1;
            else
                eff_u(i) = 1/(eff(v)); deff_u(i) = -1/(eff(v)^2);
            end
        end

        for i = 1:T_h-1
            if a_h(i) < 0
                eff_h(i) = eff(v); deff_h(i) = 1;
            else
                eff_h(i) = 1/eff(v); deff_h(i) = -1/eff(v)^2;
            end
        end
        
        predictions(v,:) = [P_h(v,:)*eff_h,P_u(v,:)*eff_u];
        dH = predictions(v,1)-hwys(v); dU = predictions(v,2)-udds(v);
        
        g(v) = 2*(dH*P_h(v,:)*deff_h+dU*P_u(v,:)*deff_u);
    end
    
    d = predictions-[hwys,udds];
    d = [d(:,1);d(:,2)];

    f = transpose(d)*d;

end
function [f,g] = likelihood2(x)
    
    eff = x;
    %k = 0.93;
    
    predictions = zeros(N,2);
    g = zeros(N,1);

    
    for v = 1:N
        eff_h = zeros(T_h,1); eff_u = zeros(T_u,1);
        deff_h = zeros(T_h,1); deff_u = zeros(T_u,1);


        for i = 1:T_u-1
            if a_u(i) < 0
                eff_u(i) = eff(v); deff_u(i) = 1;
            else
                eff_u(i) = 1/(eff(v)); deff_u(i) = -1/(eff(v)^2);
            end
        end
        
        
        for i = 1:T_h-1
            if a_h(i) < 0
                eff_h(i) = eff(v); deff_h(i) = 1;
            else
                eff_h(i) = 1/eff(v); deff_h(i) = -1/eff(v)^2;
            end
        end
        
        predictions(v,:) = [P_h(v,:)*eff_h,P_u(v,:)*eff_u];
        dH = predictions(v,1)-hwys(v);% dU = predictions(v,2)-udds(v);
        
        g(v) = 2*(dH*P_h(v,:)*deff_h);%+dU*P_u(v,:)*deff_u);
    end
    
    d = predictions-[hwys,udds];
    d = [d(:,1)];%;d(:,2)];

    f = transpose(d)*d;

end


%plot(y(:,1),y(:,2))

%{

function f = likelihood(efficiency)%,var)
    % this will store the MPGe predictions
    E = zeros(N,2);

    for k = 1:N
        Ta = targetA(k); Tb = targetB(k); Tc = targetC(k);

        eff_h = zeros(1,T_h); eff_u = zeros(1,T_u);

        for l = 1:T_h-1
            if a_h(l) <= 0
                eff_h(l) = efficiency;
            else
                eff_h(l) = 1/efficiency;
            end
        end

        for l = 1:T_u-1
            if a_u(l) <= 0
                eff_u(l) = efficiency;
            else
                eff_u(l) = 1/efficiency;
            end
        end        

        F_h = force(v_h,Ta,Tb,Tc)+mass(k)*a_h;
        F_u = force(v_u,Ta,Tb,Tc)+mass(k)*a_u;

        P_h = F_h.*v_h; P_u = F_u.*v_u;

        E(k,1) = 75384669*s_h/(eff_h*P_h); E(k,2) = 75384669*s_u/(eff_u*P_u);

    end
    %errors = E-[hwys,udds];
    errors = [E(:,1);E(:,2)]-[hwys;udds];

    %f = 1/(transpose(errors)*errors);
    %f = exp(-transpose(errors)*errors/(2*var^2))/(2*pi*var^2);
    f = -0.5*log(2*pi)-0.5*transpose(errors)*errors;
end
% {
x = 0.1:0.001:1;
y = zeros(length(x),1);

for i = 1:length(y)
    y(i) = likelihood(x(i));
end

figure(1)
plot(x,y)
 
%} 
%{
x = 0.3:0.005:1;
y = 1:0.1:40;

z = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        z(i,j) = likelihood(x(i),y(j));
    end
end

figure(1)
surf(z)

%}
end