% THIS VERSION PLOTS THE CHARGING REQUIREMENT

pPH = 6; % points per hour
l = 24*pPH;

% create a sample journey
jan = zeros(l,7); feb = zeros(l,7); mar = zeros(l,7); apr = zeros(l,7);
may = zeros(l,7); jun = zeros(l,7); jul = zeros(l,7); aug = zeros(l,7);
sep = zeros(l,7); oct = zeros(l,7); nov = zeros(l,7); dec = zeros(l,7);

% pick a number of journeys to simulate
N = 377831;

% load your vehicle
mVehicle = 3625/2.2;
mPassenger = 100;
Ta = 29.97;
Tb = 0.0713;
Tc = 0.02206;

% load distance lengths
lengths = csvread('NTS/purpose.csv',1,1); 
lengths = round(lengths(:,1)*1609.34); % convert from miles to m

% Then calculate energy expenditure of each trip
energy = zeros(7,1);
times = zeros(7,1);
%power = zeros(7,1);
chargeTime = zeros(7,1);

for i = 1:7
    v = scaleArtemis(lengths(i));
    times(i) = length(v)/60; % minutes
    times(i) = round(times(i)*pPH/60);
    energy(i) = ev(v,mVehicle,mPassenger,Ta,Tb,Tc);
    chargeTime(i) = round(energy(i)*pPH/4); % 4kW charging
    %power(i) = energy(i)*3600/length(v);
    %energy(i) = 10^5*energy(i);
end

% first sample from the month, purpose
mp = csvread('NTS/purpose_month.csv',1,1);
dp = csvread('NTS/purpose_day.csv',1,1);
tp = csvread('NTS/purpose_time.csv',1,1);

for q = 1:N
    cdf = zeros(7*12,1);
    cdf(1) = mp(1,1);

    for i = 1:7
        for j = 1:12
            if i ==1 && j ==1
                continue
            else
                cdf(12*(i-1)+j) = cdf(12*(i-1)+j-1)+mp(i,j);
            end
        end
    end

    sample = rand;

    c = 1;
    while cdf(c) <= sample
        c = c+1;
    end
    
    if c >= length(cdf)
        i = 7; j = 12;
    else
        r = rem(c,12); i = ((c-r)/12)+1;
        if r == 0
            j = 12;
        else
            j = r;
        end
    end

    month = j; purpose = i;
    % month = j, purpose = i

    % now we want to find the day of week and time of day

    % lets go for day of the week first


    % we already know the purpose is i so we can select the relevant row 

    pdf = dp(purpose,:); 
    pdf = pdf/sum(pdf); % normalise

    cdf = zeros(length(pdf),1); cdf(1) = pdf(1);

    for k = 2:length(cdf)
        cdf(k) = cdf(k-1)+pdf(k);
    end

    sample = rand;

    c = 1;
    while cdf(c)<= sample
        c = c+1;
    end

    if c >= 7
        day = 7;
    else
        day = c;
    end

    % lastly we want to get the time of day


    pdf = tp(purpose,:); 
    pdf = pdf/sum(pdf); % normalise

    cdf = zeros(length(pdf),1); cdf(1) = pdf(1);

    for k = 2:length(cdf)
        cdf(k) = cdf(k-1)+pdf(k);
    end

    sample = rand;

    c = 1;
    while cdf(c)<= sample
        c = c+1;
    end

    if c >= 24
        hour = 24;
    else
        hour = c;
    end
    
    offset = round(rand*pPH);
    
    if month == 1
        for m = 1:chargeTime(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            jan(t,day) = jan(t,day) + 4*10^5;
        end
    end


    %{
    if month == 1
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            jan(t,day) = jan(t,day) + power(purpose);
        end
    end
    
    if month == 2
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            feb(t,day) = feb(t,day) + power(purpose);
        end
    end
    
    if month == 3
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            mar(t,day) = mar(t,day) + power(purpose);
        end
    end
    
    if month == 4
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            apr(t,day) = apr(t,day) + power(purpose);
        end
    end
    
    if month == 5
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            may(t,day) = may(t,day) + power(purpose);
        end
    end
    
    if month == 6
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            jun(t,day) = jun(t,day) + power(purpose);
        end
    end
    
    if month == 7
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            jul(t,day) = jul(t,day) + power(purpose);
        end
    end
    
    if month == 8
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            aug(t,day) = aug(t,day) + power(purpose);
        end
    end
    
    if month == 9
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            sep(t,day) = sep(t,day) + power(purpose);
        end
    end
    
    if month == 10
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            oct(t,day) = oct(t,day) + power(purpose);
        end
    end
    
    if month == 11
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            nov(t,day) = nov(t,day) + power(purpose);
        end
    end
    
    if month == 12
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            dec(t,day) = dec(t,day) + power(purpose);
        end
    end
    %}
            
    %    jan(hour,day) = jan(hour,day) + energy(purpose);
        
        %{
    elseif month == 2
        feb(hour,day) = feb(hour,day) + energy(purpose);
    elseif month == 3
        mar(hour,day) = mar(hour,day) + energy(purpose);
    elseif month == 4
        apr(hour,day) = apr(hour,day) + energy(purpose);
    elseif month == 5
        may(hour,day) = may(hour,day) + energy(purpose);
    elseif month == 6
        jun(hour,day) = jun(hour,day) + energy(purpose);
    elseif month == 7
        jul(hour,day) = jul(hour,day) + energy(purpose);
    elseif month == 8
        aug(hour,day) = aug(hour,day) + energy(purpose);
    elseif month == 9
        sep(hour,day) = sep(hour,day) + energy(purpose);
    elseif month == 10
        oct(hour,day) = oct(hour,day) + energy(purpose);
    elseif month == 11
        nov(hour,day) = nov(hour,day) + energy(purpose);
    elseif month == 12
        dec(hour,day) = dec(hour,day) + energy(purpose);
    end
    %}
end

figure(1)
plot(jan)
%legend('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')

%{
figure(1)
subplot(6,2,1)
plot(jan)
title('January')
%xlabel('Hour')
%ylabel('kWh')
legend('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')

subplot(6,2,2)
plot(feb)
title('Febuary')

subplot(6,2,3)
plot(mar)
title('March')

subplot(6,2,4)
plot(apr)
title('April')

subplot(6,2,5)
plot(may)
title('May')

subplot(6,2,6)
plot(jun)
title('June')

subplot(6,2,7)
plot(jul)
title('July')

subplot(6,2,8)
plot(aug)
title('August')

subplot(6,2,9)
plot(sep)
title('September')

subplot(6,2,10)
plot(oct)
title('October')

subplot(6,2,11)
plot(nov)
title('November')

subplot(6,2,12)
plot(dec)
title('December')
%}
%{
for i = 1:7
    figure(i)
    subplot(3,4,1)
    bar(jan(:,i))
    subplot(3,4,2)
    bar(feb(:,i))
    subplot(3,4,3)
    bar(mar(:,i))
    subplot(3,4,4)
    bar(apr(:,i))
    subplot(3,4,5)
    bar(may(:,i))
    subplot(3,4,6)
    bar(jun(:,i))
    subplot(3,4,7)
    bar(jul(:,i))
    subplot(3,4,8)
    bar(aug(:,i))
    subplot(3,4,9)
    bar(sep(:,i))
    subplot(3,4,10)
    bar(oct(:,i))
    subplot(3,4,11)
    bar(nov(:,i))
    subplot(3,4,12)
    bar(dec(:,i))
end

%}
