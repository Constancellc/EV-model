
data = csvread('../EPA-Code/pureEVdata.csv',0,3);
efficiency = pure_ev;

file = 'efficiencies.csv';
data = [data,efficiency];
csvwrite(file,data)