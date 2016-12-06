v = csvread('udds.csv',0,1)*0.277778; % m/s

mVehicle = 4250/2.2;
mPassenger = 100;
Ta = 33.83;
Tb = 0.1242;
Tc = 0.03132;

energy = ev(v,mVehicle,mPassenger,Ta,Tb,Tc);