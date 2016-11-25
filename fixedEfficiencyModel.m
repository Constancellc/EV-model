function [f,g,H] = fixedEfficiencyModel()

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
x = fminunc(@likelihood,[0.3;1],options);

efficiency = x(1);
sd = exp(x(2));

predictions = [P_h*eff_h,P_u*eff_u];
error = sd*ones(N,1);


% I really want to covert back to MPGe
hwys = 75384669*s_h./hwys;
udds = 75384669*s_u./udds;

h = 75384669*s_h./predictions(:,1);
h_l = 75384669*s_h./(predictions(:,1)+error);
h_u = 75384669*s_h./(predictions(:,1)-error);

u = 75384669*s_h./predictions(:,2);
u_l = 75384669*s_h./(predictions(:,2)+error);
u_u = 75384669*s_h./(predictions(:,2)-error);

figure(1)
subplot(2,1,1)
errorbar([1:N],h,h_l,h_u)
hold on
plot([1:N],hwys,'x')

subplot(2,1,2)
errorbar([1:N],u,u_l,u_u)
hold on
plot([1:N],udds,'x')


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

function [f,g,H] = likelihood(x)
    eff = x(1); 
    var = exp(x(2));
    
    eff_h = zeros(T_h,1); eff_u = zeros(T_u,1); deff_h = zeros(T_h,1);
    deff_u = zeros(T_u,1); d2eff_h = zeros(T_h,1); d2eff_u = zeros(T_u,1);
        
    for i = 1:T_u-1
        if a_u(i) < 0
            eff_u(i) = eff; deff_u(i) = 1; d2eff_u(i) = 0;
        else
            eff_u(i) = 1/eff; deff_u(i) = -1/eff^2; d2eff_u(i) = 2/eff^3;
        end
    end

    for i = 1:T_h-1
        if a_h(i) < 0
            eff_h(i) = eff; deff_h(i) = 1; d2eff_h(i) = 0;
        else
            eff_h(i) = 1/eff; deff_h(i) = -1/eff^2; d2eff_h(i) = 2/eff^3;
        end
    end
    
    d_h = (P_h*eff_h)-hwys;
    d_u = (P_u*eff_u)-udds; d = [d_h;d_u];

    f = N*log(2*pi)+2*N*log(var)+0.5*(1/var^2)*transpose(d)*d;

    % find gradient and hessian of objective fn
    g = zeros(2,1);
    H = zeros(2,2);

    g(2) = (2*N/var)-(1/var^3)*(transpose(d)*d); % CHECKED
    g(2) = g(2)*var;

    H(2,2) = -2*N/(var^2)+(3/var^4)*transpose(d)*d;

    for i = 1:N
        g(1) = g(1) + d_h(i)*P_h(i,:)*deff_h + d_u(i)*P_u(i,:)*deff_u;
        H(1,2) = H(1,2) + d_h(i)*P_h(i,:)*deff_h + d_u(i)*P_u(i,:)*deff_u;
        H(1,1) = H(1,1) + (P_h(i,:)*deff_h)^2 + (P_u(i,:)*deff_u)^2 +...
            d_h(i)*P_h(i,:)*d2eff_h + d_u(i)*P_u(i,:)*d2eff_u;

    end

    g(1) = g(1)/(var^2); H(1,2) = -2*H(1,2)/(var^3); H(1,1) = H(1,1)/(var^2);

    H(2,1) = H(1,2);

    %x = x-0.001*H\g;
    %y(i+1,:) = x;
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