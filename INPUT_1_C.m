%INPUT FILE
%Change input data here
%-----------------------------
%Input parameters
%k=thermal conductivity (W/mK)
k=30;
%qn=heat flux (W/m^2)
qn=2e5;
%T1=prescribed temperature at the boundary (Celsius)
T1=200;
%Q=heat generation (W/m^3)
Q2=1e6; Q1=0;
%h=convective heat transfer coefficient (W/m^2K)
h=60;
%Tinf=Surrounding temperature (Celsius)
Tinf=25;
%-----------------------------
%-----------------------------
%Coordinate data (1st column-x coordinate, 2nd column-y coordinate)
coord=[0.0 0.0;
    0.5 0.0;
    0.5 0.6;
    0.0 0.6;
    0.5 0.3;
    0.8 0.3;
    0.8 0.0];
%------------------------------
%------------------------------
%Connectivity data (Each row represents the nodes of that element)
connect=[1 2 3 4;
    2 7 6 5];
%------------------------------
%------------------------------
%Gauss points for 3x3 domain integration
xi=[-sqrt(0.6) sqrt(0.6) 0]; 
eta=[-sqrt(0.6) sqrt(0.6) 0]; 
xiw=[5/9 5/9 8/9];
etaw=[5/9 5/9 8/9];
%------------------------------
%-----------------------------

