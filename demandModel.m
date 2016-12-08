function demandModel()

% create a sample journey

jan = zeros(24,7); feb = zeros(24,7); mar = zeros(24,7); apr = zeros(24,7);
may = zeros(24,7); jun = zeros(24,7); jul = zeros(24,7); aug = zeros(24,7);
sep = zeros(24,7); oct = zeros(24,7); nov = zeros(24,7); dec = zeros(24,7);

% pick a number of journeys to simulate
N = 3778310;

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
for i = 1:7
    v = scaleArtemis(lengths(i));
    energy(i) = 10^4*ev(v,mVehicle,mPassenger,Ta,Tb,Tc);
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

    if month == 1
        jan(hour,day) = jan(hour,day) + energy(purpose);
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
end


figure(1)
subplot(3,4,1)
plot(jan)
title('January')
xlabel('Hour')
ylabel('kWh')
legend('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')

subplot(3,4,2)
plot(feb)
title('Febuary')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,3)
plot(mar)
title('March')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,4)
plot(apr)
title('April')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,5)
plot(may)
title('May')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,6)
plot(jun)
title('June')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,7)
plot(jul)
title('July')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,8)
plot(aug)
title('August')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,9)
plot(sep)
title('September')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,10)
plot(oct)
title('October')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,11)
plot(nov)
title('November')
xlabel('Hour')
ylabel('kWh')

subplot(3,4,12)
plot(dec)
title('December')
xlabel('Hour')
ylabel('kWh')
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
end