% lets scale artemis
function v = scaleArtemis(s)%also t?
% inputs: time in s, distance in m

% output: new drive cycle

v0 = csvread('artemis_urban.csv')*0.277778; %m/s

T = length(v0);

dist = sum(v0);

v = [];
runningDisp = 0;

if s <= dist
    i = 1;
    while runningDisp <= s
        v = [v;v0(i)];
        runningDisp = runningDisp + v0(i);
        i = i+1;
    end
else
    N = floor(s/dist);
    for i = 1:N
        v = [v;v0];
        runningDisp = runningDisp + dist;
    end
    
    i = 1;
    while runningDisp <= s
        v = [v;v0(i)];
        runningDisp = runningDisp + v0(i);
        i = i+1;
    end
end
        
% the code below would also scale for length but I don't think this is a
% good idea
        
%{
if t/T <= 1
    v = v0(1:t);
else
    v = []; c = 0;
    for i = 1:floor(t/T)
        v = [v;v0];
        c = c+1;
    end
    v = [v;v0(1:t-T*c)];
end

v = (s/sum(v))*v;
%}
end

