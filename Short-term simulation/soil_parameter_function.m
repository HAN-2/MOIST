function [X2]=soil_parameter_function(X)
poi=X.poi;
if poi==0
    poi=[];
end
n=X.n;
[multi_row,~]=size(poi);
layered_he=                [-0.03067*1; -0.03067*1.5; -0.03067*0.8];%X.hexx1;%; % 3.07,2.1,                      %[-1/3.07;-1/2.07]   [-1/3.07;-1/2.83;-1/3.07];
layered_lamta=             [0.1494*0.9;  0.1494*1.1;  0.1494*0.8];%X.lamtaxx2;%;    %0.15,0.12,                        % [0.1523;0.12;0.15];[0.123;0.17]; [0.1523] [0.1523;0.14;0.1523];
layered_yeeta=             [2/layered_lamta(1,1)+3;2/layered_lamta(2,1)+3;2/layered_lamta(3,1)+3];%2/X.lamtaxx2+3;%%[9.14;9.14]; []
layered_ksat=              [3.5e-5;    3.5e-5*1.5;   3.5e-5*1];%X.ksatxx3;% %[1.16e-5;7.2e-6;9.16e-6];  %8.16e-5                       %[6e-6;1.35e-5] [1.1e-5] [1.1e-5;5.3e-6;1.1e-5]; [1.1e-5;5.3e-6;0.6e-5];
layered_thetasat=          [0.3*1;     0.3*0.9;      0.3-0.02    ];%X.thetasatxx5;%    %[0.48;0.48;0.46];                           %[0.471] [0.451;0.5;0.461];
layered_thetares=          [0.005;     0.005*10;     0.005*6];%X.thetaresxx4;%                              %[0.01;0.01]; %[0.04] [0.04;0.03;0.04];


he=      zeros(n,1); % pressure head at air entry, m
lamda=   zeros(n,1); % shape parameter of soil retention curve
yeeta=   zeros(n,1); % shape parameter of soil retention curve
ksat=    zeros(n,1); % saturated hydrolic conductivity, m/s
theta_res=zeros(n,1);% volumetric residence water content, m3/m3
theta_sat=zeros(n,1);% saturated volumetric water content, m3/m3

he(1:end,1)=      layered_he(1,1);
lamda(1:end,1)=   layered_lamta(1,1);
yeeta(1:end,1)=   layered_yeeta(1,1);
ksat(1:end,1)=    layered_ksat(1,1);
theta_res(1:end,1)=layered_thetares(1,1);
theta_sat(1:end,1)=layered_thetasat(1,1);

for i=1:multi_row
    he(poi(i)+1:end,1)=      layered_he(i+1,1);
    lamda(poi(i)+1:end,1)=   layered_lamta(i+1,1);
    yeeta(poi(i)+1:end,1)=   layered_yeeta(i+1,1);
    ksat(poi(i)+1:end,1)=    layered_ksat(i+1,1);
    theta_res(poi(i)+1:end,1)=layered_thetares(i+1,1);
    theta_sat(poi(i)+1:end,1)=layered_thetasat(i+1,1);
end

thetac=0.2;
thetam=0.393;% content of other material
thetaq=0.243;% content of quatze
thetan=0.529;% content of soild
thetao=0;% content of organic matter

rho = 1000;
rhob = 1300;
%mim

mass_transfer_coefficient = 0;%0.15/86400;
f = 0.4;
kad = 0;%0.5;
coe = theta_res(1) + (1 - f) * rhob/1000. * kad;

ra=(log((1+0.001)/0.001)+0)^2/0.1681/0.185;                                % air resistance to water vapor, s/m  from hydrus B5  %%%%1110
rbh=ra;
X2=struct('lamda',lamda,'yeeta',yeeta,'he',he,'theta_sat',theta_sat,'theta_res',theta_res,'ksat',ksat,...
    'thetan',thetan,'thetao',thetao,'thetac',thetac,'thetam',thetam,'thetaq',thetaq,...
    'rho',rho,'rhob',rhob,...
    'mass_transfer_coefficient',mass_transfer_coefficient,'coe',coe,'f',f,'kad',kad);
end