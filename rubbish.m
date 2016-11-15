data = csvread('../EPA-Code/ev_data.csv',0,2);

%{
 1 |   engine size   |    l
 2 |   horsepower    |
 3 | # driven wheels |    -
 4 |   kerb weight   |   lbs
 5 |    axle ratio   |    -
 6 |    n/v ratio    |    -
 7 |   drive cycle   |    *
 8 | fuel efficiency |   MPGe
 9 |  Target coeff A |   lbf
10 |  Target coeff B |  lbf/mph
11 |  Target coeff C | lbf/mph^2

* 1 => CD Highway     2 => CD UDDS     0 => Other
%}

mass = data(:,4)/2.2;
mpge = data(:,8);
targetA = data(:,9);
targetC = data(:,11);