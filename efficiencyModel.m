function efficiencyModel()

% Getting the efficiencies
y = csvread('efficiencies.csv');

% Now getting the drive cycles
v_h = csvread('highways.csv',0,1)*0.277778; % m/s
s_h = sum(v_h)/1000; %km
T_h = length(v_h); a_h = zeros(T_h,1); vh_av = s_h*1000/T_h;

for j = 1:T_h-1
    a_h(j) = v_h(j+1)-v_h(j);
end

v_u = csvread('udds.csv',0,1)*0.277778; % m/s
s_u = sum(v_u)/1000; %km
T_u = length(v_u); a_u = zeros(T_u,1); vu_av = s_u*1000/T_u;

for j = 1:T_u-1
    a_u(j) = v_u(j+1)-v_u(j);
end

% Now getting the vehicle data
[mass,targetA,targetB,targetC,~,~] = getPureEVData;

N = length(targetC);

a_rms = [sqrt(transpose(a_h)*a_h/T_h);sqrt(transpose(a_u)*a_u/T_u)];

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');
options.MaxFunctionEvaluations = 3000;
options.MaxIterations = 4000;
[x,fval] = fminunc(@squareError,[ones(2,1)],options)

x(1)
x(2)
x(3)

plot([1:N],[d(1:N),d(N+1:end)])

    function [f,g] = squareError(x)
        
        k1 = x(1); k2 = x(2); %k3 = x(3);
        g = zeros(2,1);
        d = zeros(2*N,1);
        
        for i = 1:N
            %{
            loss_h = k1 + k2*vh_av + k3*mass(i)*a_rms(1);
            loss_u = k1 + k2*vu_av + k3*mass(i)*a_rms(2);
            
            d(i) = y(i)-1/loss_h; d(i+N) = y(i+N)-1/loss_u;
            
            g(1) = g(1) + 2*(d(i)/loss_h^2 + d(i+N)/loss_u^2);
            g(2) = g(2) + 2*(d(i)*vh_av/loss_h^2 + d(i+N)*vu_av/loss_u^2);
            g(3) = g(3) + 2*mass(i)*(d(i)*a_rms(1)/loss_h^2 +...
                d(i+N)*a_rms(2)/loss_u^2);
            %}
            
            eff_h = k1 + k2/(mass(i)*a_rms(1));%/vh_av + k3/(mass(i)*a_rms(1));
            eff_u = k1 + k2/(mass(i)*a_rms(2));%/vu_av + k3/(mass(i)*a_rms(2));
            
            d(i) = y(i)-eff_h; d(i+N) = y(i+N)-eff_u;
            
            g(1) = g(1) - 2*(d(i) + d(i+N));
            %g(2) = g(2) - 2*(d(i)/vh_av + d(i+N)/vu_av);
            g(2) = g(2) - (2/mass(i))*(d(i)/a_rms(1) + d(i+N)/a_rms(2));
            
        end
        
        f = transpose(d)*d;
    end
end
