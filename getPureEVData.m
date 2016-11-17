function [mass,targetA,targetB,targetC,mpge] = getPureEVData(cycle)

data = csvread('../EPA-Code/ev_data2.csv',0,2);

%{
 1 |   horsepower    |
 2 | # driven wheels |    -
 3 |   kerb weight   |   lbs
 4 |    axle ratio   |    -
 5 |    n/v ratio    |    -
 6 |   drive cycle   |    *
 7 | fuel efficiency |   MPGe
 8 |  Target coeff A |   lbf
 9 |  Target coeff B |  lbf/mph
10 |  Target coeff C | lbf/mph^2

* 1 => CD Highway     2 => CD UDDS     0 => Other
%}

indicies = [];

for i = 1:length(data)
    if data(i,6) == cycle
        indicies = [indicies,i];
    end
end

mass = data(indicies,3)/2.2;
mpge = data(indicies,7);
targetA = data(indicies,8);
targetB = data(indicies,9);
targetC = data(indicies,10);

min(targetA)
max(targetA)
sum(targetA)/length(indicies)
end

