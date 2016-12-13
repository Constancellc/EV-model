function untitled10
% THIS VERSION PLOTS THE CONTINUOUS-ISH POWER REQUIREMENT

pPH = 6; % points per hour
l = 24*pPH;

% First initialise a structure for each of the regions
UC = struct('jan',zeros(l,7),'feb',zeros(l,7),'mar',zeros(l,7),'apr',...
    zeros(l,7),'may',zeros(l,7),'jun',zeros(l,7),'jul',zeros(l,7),'aug',...
    zeros(l,7),'sep',zeros(l,7),'oct',zeros(l,7),'nov',zeros(l,7),'dec',...
    zeros(l,7));
UCT = struct('jan',zeros(l,7),'feb',zeros(l,7),'mar',zeros(l,7),'apr',...
    zeros(l,7),'may',zeros(l,7),'jun',zeros(l,7),'jul',zeros(l,7),'aug',...
    zeros(l,7),'sep',zeros(l,7),'oct',zeros(l,7),'nov',zeros(l,7),'dec',...
    zeros(l,7));
RT = struct('jan',zeros(l,7),'feb',zeros(l,7),'mar',zeros(l,7),'apr',...
    zeros(l,7),'may',zeros(l,7),'jun',zeros(l,7),'jul',zeros(l,7),'aug',...
    zeros(l,7),'sep',zeros(l,7),'oct',zeros(l,7),'nov',zeros(l,7),'dec',...
    zeros(l,7));
RV = struct('jan',zeros(l,7),'feb',zeros(l,7),'mar',zeros(l,7),'apr',...
    zeros(l,7),'may',zeros(l,7),'jun',zeros(l,7),'jul',zeros(l,7),'aug',...
    zeros(l,7),'sep',zeros(l,7),'oct',zeros(l,7),'nov',zeros(l,7),'dec',...
    zeros(l,7));

data = [UC UCT RT RV];

% matrix with rows: purpose and cols: region
lengths = csvread('NTS/FINALregionTypePurposeLength.csv',1,1)*1609.34;
% pick a number of journeys to simulate
N = round(sum(sum(csvread('NTS/FINALregionTrips.csv',4,1)))*10^-5);

% load your vehicle
mVehicle = 3625/2.2;
mPassenger = 100;
Ta = 29.97;
Tb = 0.0713;
Tc = 0.02206;

% load distance lengths
%lengths = csvread('NTS/purpose.csv',1,1); 
%lengths = round(lengths(:,1)*1609.34); % convert from miles to m

% Then calculate energy expenditure of each trip
energy = zeros(7,9);
times = zeros(7,9);
power = zeros(7,9);

for i = 1:7
    for j = 1:4
        v = scaleArtemis(lengths(i,j));
        times(i,j) = length(v)/60; % minutes
        times(i,j) = round(times(i,j)*pPH/60);
        energy(i,j) = 10^5*ev(v,mVehicle,mPassenger,Ta,Tb,Tc);
        power(i,j) = energy(i,j)*3600/length(v);
    end
end

%{
samples = [12938;13829;2982;2934];
samples = samples/sum(samples);

Nuc = samples(1)*N; Nuct = samples(2)*N; Nrt = samples(3)*N; Nrv = 
%}
% first sample from the month, purpose
mp = csvread('NTS/FINALpurposeMonth.csv',1,1);
dp = csvread('NTS/FINALpurposeDay.csv',1,1);
tp = csvread('NTS/FINALpurposeStartHour.csv',1,1);
rp = csvread('NTS/FINALregionTypePurpose.csv',1,1);

function out = sampleGivenPurpose(joint,purpose,max)
    pdf = joint(purpose,:); 
    pdf = pdf/sum(pdf); % normalise

    cdf_ = zeros(length(pdf),1); cdf_(1) = pdf(1);

    for k = 2:length(cdf_)
        cdf_(k) = cdf_(k-1)+pdf(k);
    end

    sample_ = rand;

    c_ = 1;
    while cdf_(c_)<= sample_
        c_ = c_+1;
    end

    if c_ >= max
        out = max;
    else
        out = c_;
    end
end

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
    while cdf(c) <= sample && c<length(cdf)
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
    
    % lets find the region
    
    region = sampleGivenPurpose(rp,purpose,4);

    day = sampleGivenPurpose(dp,purpose,7);

    % lastly we want to get the time of day

    hour = sampleGivenPurpose(tp,purpose,24);
    
    offset = round(rand*pPH);%/pointsPerHour;
    
    if month == 1
        for m = 1:times(purpose,region)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).jan(t,day) = data(region).jan(t,day) + power(purpose,region);
            %jan(t,day) = jan(t,day) + power(purpose);
        end
    end
    
    if month == 2
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).feb(t,day) = data(region).feb(t,day) + power(purpose,region);
            %feb(t,day) = feb(t,day) + power(purpose);
        end
    end
    
    if month == 3
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).mar(t,day) = data(region).mar(t,day) + power(purpose,region);
            %mar(t,day) = mar(t,day) + power(purpose);
        end
    end
    
    if month == 4
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).apr(t,day) = data(region).apr(t,day) + power(purpose,region);
            %apr(t,day) = apr(t,day) + power(purpose);
        end
    end
    
    if month == 5
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).may(t,day) = data(region).may(t,day) + power(purpose,region);
            %may(t,day) = may(t,day) + power(purpose);
        end
    end
    
    if month == 6
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).jun(t,day) = data(region).jun(t,day) + power(purpose,region);
            %jun(t,day) = jun(t,day) + power(purpose);
        end
    end
    
    if month == 7
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).jul(t,day) = data(region).jul(t,day) + power(purpose,region);
            %jul(t,day) = jul(t,day) + power(purpose);
        end
    end
    
    if month == 8
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).aug(t,day) = data(region).aug(t,day) + power(purpose,region);
            %aug(t,day) = aug(t,day) + power(purpose);
        end
    end
    
    if month == 9
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).sep(t,day) = data(region).sep(t,day) + power(purpose,region);
            %sep(t,day) = sep(t,day) + power(purpose);
        end
    end
    
    if month == 10
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).oct(t,day) = data(region).oct(t,day) + power(purpose,region);
            %oct(t,day) = oct(t,day) + power(purpose);
        end
    end
    
    if month == 11
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).nov(t,day) = data(region).nov(t,day) + power(purpose,region);
            %nov(t,day) = nov(t,day) + power(purpose);
        end
    end
    
    if month == 12
        for m = 1:times(purpose)
            t = (hour-1)*(pPH)+offset+m;
            if t > 24*pPH
                t = t-24*pPH;
            end
            data(region).dec(t,day) = data(region).dec(t,day) + power(purpose,region);
            %dec(t,day) = dec(t,day) + power(purpose);
        end
    end
end

titles = {'Urban Conurbation','Urban City and Town','Rural Town and Fringe'...
    'Rural Village, Hamlet and Isolated Dwelling'};
figure(1)
for i = 1:4
    subplot(4,1,i)
    plot(data(i).feb)
    title(titles(i))
end

%figure(1)
%plot(jan)
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
end