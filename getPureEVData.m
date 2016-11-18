function [mass,targetA,targetB,targetC,hwys,udds] = getPureEVData

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

mass = data(:,3)/2.2; %kg
targetA = data(:,6);
targetB = data(:,7);
targetC = data(:,8);

hwys = data(:,9);
udds = data(:,10);

size(data)
end

