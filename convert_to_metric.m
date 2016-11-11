cycle = 'highways.csv';

a = csvread(cycle);

for i = 1:length(a)
    a(i,2) = a(i,2)*1.60934;
end

csvwrite(cycle,a)
