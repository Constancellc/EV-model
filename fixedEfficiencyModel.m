function E = fixedEfficiencyModel()

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

% constructing hessian and grad x = [eff,var]

P_h = zeros(N,T_h); P_u = zeros(N,T_u);
eff_h = zeros(T_h,1); eff_u = zeros(T_u,1); deff_h = zeros(T_h,1);
deff_u = zeros(T_u,1); d2eff_h = zeros(T_h,1); d2eff_u = zeros(T_u,1);
% eff, deff, d2eff, P, d

x = [0.7;1];
y = zeros(11,2); y(1,:) = x;
for m = 1:10
    eff = x(1); var = x(2);
    
    for i = 1:N
        Ta = targetA(i); Tb = targetB(i); Tc = targetC(i);
        
        F_h = force(v_h,Ta,Tb,Tc)+mass(i)*a_h;
        F_u = force(v_u,Ta,Tb,Tc)+mass(i)*a_u;

        P_h(i,:) = F_h.*v_h; P_u(i,:) = F_u.*v_u;
    end
    
    P = [P_h;P_u];
        
    for i = 1:T_u-1
        if a_u(i) < 0
            eff_u(i) = eff; deff_u(i) = 1; d2eff_u(i) = 0;
        else
            eff_u(i) = 1/eff; deff_u(i) = -1/eff^2; d2eff_u = 2/eff^3;
        end
    end

    for i = 1:T_h-1
        if a_h(i) < 0
            eff_h(i) = eff; deff_h(i) = 1; d2eff_h(i) = 0;
        else
            eff_h(i) = 1/eff; deff_h(i) = -1/eff^2; d2eff_h = 2/eff^3;
        end
    end
    
    deff = [deff_h;deff_u]; d2eff = [d2eff_h;d2eff_u];
    
    d_h = P_h*eff_h-hwys; d_u = P_u*eff_u-udds; d = [d_h;d_u];   

    % find gradient and hessian of objective fn
    g = zeros(2,1);
    H = zeros(2,2);

    g(2) = -(2*N/var)+(1/var^3)*(transpose(d)*d);

    H(1,1) = (1/var^2)*(transpose(P_h*d_h)*d2eff_h+...
        transpose(P_u*d_u)*d2eff_u);
    H(1,2) = -(2/var^3)*(transpose(d_h)*(P_h*deff_h)+...
        transpose(d_u)*(P_u*deff_u));
    H(2,1) = H(1,2);
    H(2,2) = 2*N/(var^2)-(3/var^4)*transpose(d)*d;

    for i = 1:n
        g(1) = g(1) + (1/var^2)*(d_h(i)*transpose(deff_h)*P_h(i,:)+...
            d_u(i)*transpose(deff_u)*P_u(i,:));
        % HAVEN'T FIGURED OUT ONE BELOW
        H(1,1) = H(1,1) - (1/var^2)*((transpose(deff)*P(i,:))^2;
    end
    
    x = x+H\g;
    y(i+1,:) = x;
end

plot(y(:,1),y(:,2))

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