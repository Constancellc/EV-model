function rmse = modelTesting
%testingDataSet
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

data = csvread('../EPA-Code/pureEVdata.csv',0,3);

%{
 1 |   horsepower    |
 2 | # driven wheels |    -
 3 |   kerb weight   |   lbs
 4 |    axle ratio   |    -
 5 |    n/v ratio    |    -
 6 |  Target coeff A |   lbf
 7 |  Target coeff B |  lbf/mph
 8 |  Target coeff C | lbf/mph^2
 9 |   Highways fe   |   MPGe
10 |     UDDS fe     |   MPGe

%}
M = 10;
%testing = [3,6,12,16,20,23,24,29,35,38];
testing = [];
training = [];

for i = 1:length(data)
    if rand <= M/length(data)
        testing = [testing;i];
    else
        training = [training;i];
    end
end

N = length(training)
M = length(testing);

%horsepower = data(:,1); mass = data(:,3)/2.2; axleRatio = data(:,4);
%targetA = data(:,6); targetB = data(:,7); targetC = data(:,8);
hwys0 = data(:,9); udds0 = data(:,10);

T_training = data(training,[6,7,8]); T_testing = data(testing,[6,7,8]);
X_training = data(training,[1,3,4]); X_testing = data(testing,[1,3,4]);

mTrain = X_training(:,2); mTest = X_testing(:,2);

% converting from MPGe to energy consumption
hwys = 75384669*s_h./hwys0;
udds = 75384669*s_u./udds0;

y_h = hwys(training); y_u = udds(training);

% in N
force =@(v,Ta,Tb,Tc) 4.44822*(Ta + Tb*(v./0.44704)+Tc*(v./0.44707).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   THIS SECTION DEALS WITH TRAINING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_h = zeros(N,T_h); P_u = zeros(N,T_u);

for j = 1:N
    Ta = T_training(j,1); Tb = T_training(j,2); Tc = T_training(j,3);

    F_h = force(v_h,Ta,Tb,Tc)+mTrain(j)*a_h;
    F_u = force(v_u,Ta,Tb,Tc)+mTrain(j)*a_u;

    P_h(j,:) = F_h.*v_h; P_u(j,:) = F_u.*v_u;
end

% Start with single efficiency case

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');
x = fminunc(@singleEfficiency,[0.3],options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   THIS SECTION DEALS WITH TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now let's predict for the test case

P_hT = zeros(M,T_h); P_uT = zeros(M,T_u);

for j = 1:M
    Ta = T_testing(j,1); Tb = T_testing(j,2); Tc = T_testing(j,3);

    F_h = force(v_h,Ta,Tb,Tc)+mTest(j)*a_h;
    F_u = force(v_u,Ta,Tb,Tc)+mTest(j)*a_u;

    P_hT(j,:) = F_h.*v_h; P_uT(j,:) = F_u.*v_u;
end
    
effh = zeros(T_h,1); effu = zeros(T_u,1);

for i = 1:T_u-1
    if a_u(i) < 0
        effu(i) = x(1);
    else
        effu(i) = 1/x(1);
    end
end

for i = 1:T_h-1
    if a_h(i) < 0
        effh(i) = x(1);
    else
        effh(i) = 1/x(1);
    end
end

predH = P_hT*effh; predU = P_uT*effu;

% converting back to MPGe
h = 75384669*s_h./predH; u = 75384669*s_u./predU;

error_h = h-hwys0(testing); error_u = u-udds0(testing);

avError1 = sqrt(transpose(error_h)*error_h/M);
avError2 = sqrt(transpose(error_u)*error_u/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   THIS SECTION PLOTS THE GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1)
subplot(2,1,1)
plot([1:M],h,'x')
hold on
plot([1:M],hwys0(testing),'o')
textbox1 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.25 0.65 .3 .05]);
set(textbox1, 'String', ['rms error: ' num2str(avError1) ' MPGe']);
title('Predicted vs. Observed MPGe on Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')

subplot(2,1,2)
plot([1:M],u,'x')
hold on
plot([1:M],udds0(testing),'o')
textbox2 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.25 0.15 .3 .05]);
set(textbox2, 'String', ['rms error: ' num2str(avError2) ' MPGe']);
title('Predicted vs. Observed MPGe on Urban Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NOW LOOKING AT THE VEHICLE SPECIFIC EFFICIENCY CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('fminunc','Algorithm','trust-region','GradObj','on','DerivativeCheck','on');
x = fminunc(@vehicleSpecificEfficiency,[ones(N,1)],options);

% and testing, now we need a way of running KNN to estimate the
% efficiencies of the test cases. Not sure I can bother to do this now..

% I want a function which I can give the training inputs and efficiencies
% as well as the test inputs and it will retun the testing efficiencies

effPred = KNN(X_training,X_testing,x);%(2:end));
predH = zeros(M,1); predU = zeros(M,1);

for j = 1:M
    effh = zeros(T_h,1); effu = zeros(T_u,1);
    for i = 1:T_u-1
        if a_u(i) < 0
            effu(i) = effPred(j);
        else
            effu(i) = 1/effPred(j);
        end
    end

    for i = 1:T_h-1
        if a_h(i) < 0
            effh(i) = effPred(j);
        else
            effh(i) = 1/effPred(j);
        end
    end

    predH(j) = P_hT(j,:)*effh; predU(j) = P_uT(j,:)*effu;
end

% converting back to MPGe
h2 = 75384669*s_h./predH; u2 = 75384669*s_u./predU;

error_h = h2-hwys0(testing); error_u = u2-udds0(testing);

avError3 = sqrt(transpose(error_h)*error_h/M);
avError4 = sqrt(transpose(error_u)*error_u/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   THIS SECTION PLOTS THE GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,1,1)
bar([1:M],[h,h2])
hold on
scatter([1:M],hwys0(testing),40,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
%plot([1:M],hwys0(testing),'x')
title('Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')
legend('show')
legend('Model 1','Model 2','Observed')

subplot(2,1,2)
bar([1:M],[u,u2])
hold on
scatter([1:M],udds0(testing),40,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
%plot([1:M],udds0(testing),'x')
title('Urban Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')
legend('show')
legend('Model 1','Model 2','Observed')

textbox1 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.91 0.74 .09 .1]);
set(textbox1, 'String', ['Model 1 RMSE: ' num2str(round(avError1)) ' MPGe']);
textbox2 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.91 0.62 .09 .1]);
set(textbox2, 'String', ['Model 2 RMSE: ' num2str(round(avError3)) ' MPGe']);
textbox3 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.91 0.27 .09 .1]);
set(textbox3, 'String', ['Model 1 RMSE: ' num2str(round(avError2)) ' MPGe']);
textbox4 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.91 0.15 .09 .1]);
set(textbox4, 'String', ['Model 2 RMSE: ' num2str(round(avError4)) ' MPGe']);

rmse = [avError1;avError2;avError3;avError4];

%{
subplot(2,2,1)
bar([1:M],[h,hwys0(testing)])
title('Homogenous Efficiency on Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')


subplot(2,2,2)
bar([1:M],[h2,hwys0(testing)])
title('Vehicle Specific Efficiency on Highways Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')

subplot(2,2,3)
bar([1:M],[u,udds0(testing)])
title('Homogenous Efficiency on Urban Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')

subplot(2,2,4)
bar([1:M],[u2,udds0(testing)])
title('Vehicle Specific Efficiency on Urban Drive Cycle')
ylabel('MPGe')
xlabel('Vehicle No.')


textbox1 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.15 0.5 .3 .025]);
set(textbox1, 'String', ['rms error: ' num2str(avError1) ' MPGe']);
textbox2 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.58 0.025 .3 .025]);
set(textbox2, 'String', ['rms error: ' num2str(avError2) ' MPGe']);
textbox3 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.15 0.025 .3 .025]);
set(textbox3, 'String', ['rms error: ' num2str(avError3) ' MPGe']);
textbox4 = uicontrol('Style', 'text', 'Units', 'norm','Position',[0.58 0.5 .3 .025]);
set(textbox4, 'String', ['rms error: ' num2str(avError4) ' MPGe']);
%}

    function [f,g] = singleEfficiency(x)
        var = 1;
        eff = x(1); %var = exp(x(2));
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
    
        d_h = (P_h*eff_h)-y_h;
        d_u = (P_u*eff_u)-y_u; d = [d_h;d_u];

        f = N*log(2*pi)+2*N*log(var)+0.5*(1/var^2)*transpose(d)*d;

        % find gradient and hessian of objective fn
        %g = zeros(2,1);
        g = zeros(1,1);
        

        %g(2) = (2*N/var)-(1/var^3)*(transpose(d)*d); % CHECKED
        %g(2) = g(2)*var;

        for i = 1:N
            g(1) = g(1) + d_h(i)*P_h(i,:)*deff_h + d_u(i)*P_u(i,:)*deff_u;

        end

        g(1) = g(1)/(var^2); 

    end

    function [f,g] = vehicleSpecificEfficiency(x)
        
    var = 1;%exp(x(1));
    eff = x;%(2:end);
    
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
        dH = predictions(v,1)-y_h(v); dU = predictions(v,2)-y_u(v);
        
        %g(v+1) = (1/var^2)*(dH*P_h(v,:)*deff_h+dU*P_u(v,:)*deff_u);
        g(v) = (1/var^2)*(dH*P_h(v,:)*deff_h+dU*P_u(v,:)*deff_u);

    end
    
    d = predictions-[y_h,y_u];
    d = [d(:,1);d(:,2)];

    f = N*log(2*pi)+2*N*log(var)+0.5*(1/var^2)*transpose(d)*d;

    %g(1) = (2*N/var)-(1/var^3)*(transpose(d)*d); % CHECKED
    %g(1) = g(1)*var;

    end
end
