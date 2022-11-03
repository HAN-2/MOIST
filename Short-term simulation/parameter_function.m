function [X1,X3,X4] = parameter_function(X,i)
global surface_flag
persistent qls_pre qvs_pre qTs_pre Ts_pre hrs_pre cvs_pre thetas_pre G_pre rbh_pre rbw_pre qevap_pre Ep_flag_pre
r0=10;
i=i/i+1;
if X.theta(1,i-1)<0.15
    rs(1,1)=r0*exp(35.65*(0.15-X.theta(1,i-1)));
else
    rs(1,1)=r0;
end
[X2]=soil_parameter_function(X);

[X4]=isotope_parameter_function(X);


Cp=1200;
CV=1.8*10^6;
CW=4.18*10^6;
g=9.8;                                                                     % gravitational acceleration, m/s2
M=0.018;                                                                   % molar weight of water, kg/mol
R=8.314;
Dv0=2.12e-5;
Tcoe = 0;
i=X.i;
i=i/i+1;
S=(X.h(:,i-1)./X2.he).^(-X2.lamda);
k=X2.ksat.*S.^X2.yeeta;
C= (X2.theta_sat - X2.theta_res) .* (-X2.lamda) .* (X.h(:,i-1) ./ X2.he) .^ (-X2.lamda - 1) ./ X2.he;
X.theta(:,i-1)=S.*(X2.theta_sat-X2.theta_res)+X2.theta_res;


sat_index = find(X.h(:,i-1) >= X2.he(:,1));
X.theta(sat_index, i-1) = X2.theta_sat(sat_index, 1);
S(sat_index, 1) = (X.theta(sat_index, i-1) - X2.theta_res(sat_index, 1)) ./ (X2.theta_sat(sat_index, 1) - X2.theta_res(sat_index, 1));
k(sat_index, 1) = X2.ksat(sat_index, 1) .* S(sat_index, 1) .^ X2.yeeta(sat_index, 1);

C(sat_index, 1) = 0;%S(sat_index, 1) * 100;
cvsat = 610.78 .* exp(17.27 .* X.T(:,i-1) ./ (X.T(:,i-1) + 237.3)) .* M ./ R ./ X2.rho ./ (X.T(:,i-1) + 273.15);
dcvsatdT = -(610.78 * exp((17.27 .* X.T(:,i-1)) ./ (237.3 + X.T(:,i-1))) .* M .*...
    (-1.0631 .* 10 .^ 6 - 3623.57 .* X.T(:,i-1) + X.T(:,i-1) .^ 2)) ./...
    (R .* X2.rho .* (237.3 + X.T(:,i-1)) .^ 2 .* (273.15 + X.T(:,i-1)).^2);
s = dcvsatdT;
hr = exp(M .* X.h(:,i-1) .* g ./ R ./ (X.T(:,i-1)+273.15));

dhrdh = exp(M .* X.h(:,i-1) .* g ./ R ./ (X.T(:,i-1)+273.15)) .* M .* g ./ R ./ (X.T(:,i-1)+273.15);
dhrdT = -((exp((g .* X2.he .* M .* S .^ (-1 ./ X2.lamda)) ./...
    (R .* (273.15 + X.T(:,i-1)))) .* g .* X2.he .* M .* S .^ (-1 ./ X2.lamda)) ./(R .* (273.15 + X.T(:,i-1)) .^ 2));
dhrdh(sat_index, 1) = 0;
dhrdT(sat_index, 1) = 0;
hr(sat_index, 1) = 1;

lamdaE = (2.501 - 2.361 .* 10 ^ (-3) .* X.T(:,i-1)) .* 10 .^ 6;
Csoil = (1.92 .* X2.thetan + 2.51 .* X2.thetao + 4.18 .* X.theta(:,i-1)) .* 10 ^ 6;
kH = (0.65 - 0.78 * X2.rhob / 1000 + 0.6 * (X2.rhob / 1000) .^ 2) + (1.06 * (X2.rhob / 1000)) * X.theta(:,i-1) -...
    ((0.65 - 0.78 * X2.rhob / 1000 + 0.6 * (X2.rhob / 1000) .^ 2) - (0.03 + 0.1 * (X2.rhob / 1000) .^ 2)) *...
    exp(-((1 + 2.6 * X2.thetac .^ (-0.5)) * X.theta(:,i-1)) .^ 4);
thao = ((X2.theta_sat-X.theta(:,i-1)).^(7/3))./((X2.theta_sat).^2);
thaol=((X.theta).^(7/3))./((X2.theta_sat).^2);%0.66*(theta./thetasat).^(8/3);




Dv=Dv0.*10.^5./101325.*((X.T(:,i-1)+273.15)./273.15).^1.88.*thao.*(X2.theta_sat-X.theta(:,i-1));

cva=610.78*exp(17.27*X.Ta/(X.Ta+237.3))*M/R/X2.rho/(X.Ta+273.15)*X.hra;

X3=struct('h',X.h,'T',X.T,'cvsat',cvsat,'hr',hr,'s',s,'Ep',X.Ep,...
    'frozen_flag',X.frozen_flag,'hcria',X.hcria,'Dv',Dv,'R',R,'k',k,...
    'M',M,'g',g,'dz',X.dz,'rho',X2.rho,'kH', kH,'cva', cva,'Cp', Cp,...
    'rs',rs,'Ta',X.Ta,'hra',X.hra,'u',X.u,'Rnet',X.Rnet,'root_flag',X.root_flag,'he',X2.he,'i',i,...
    'theta_sat',X2.theta_sat,'theta_res',X2.theta_res,'lamda',X2.lamda,'yeeta',X2.yeeta,'ksat',X2.ksat,'itm',X.itm,'qprec',X.qprec);

if surface_flag
    [qls,qvs,qTs,...
        Ts,hrs,cvs,thetas,...
        G,rbh,rbw,...
        qevap,Ep_flag]=surface_function(X3,i,X2,X);

    qls_pre=qls;qvs_pre=qvs;qTs_pre = qTs;Ts_pre=Ts;hrs_pre=hrs;cvs_pre=cvs;
    thetas_pre=thetas;G_pre=G;rbh_pre=rbh;rbw_pre=rbw;qevap_pre=qevap;Ep_flag_pre=Ep_flag;
    surface_flag = 0;
else
    qls=qls_pre;qvs=qvs_pre;qTs=qTs_pre;Ts=Ts_pre;hrs=hrs_pre;cvs=cvs_pre;
    thetas=thetas_pre;G=G_pre;rbh=rbh_pre;rbw=rbw_pre;qevap=qevap_pre;Ep_flag=Ep_flag_pre;
end
if thetas > X2.theta_sat(1)
    thetas = X2.theta_sat(1);
end
alpha = exp(-(X4.aca./(X.T(:,i-1)+273.15).^2+X4.acb./(X.T(:,i-1)+273.15)+X4.acc));
alphas= exp(-(X4.aca./(Ts+273.15).^2+X4.acb./(Ts+273.15)+X4.acc));%1;%
cvsats = 610.78 .* exp(17.27 .* Ts ./ (Ts + 237.3)) .* M ./ R ./ X2.rho ./ (Ts + 273.15);
if X.pTest==1||X.pTest==2
    alpha =1;
    alphas=1;
end
dalphadT = exp(-(X4.aca ./ (X.T(:,i-1)+273.15) ./ (X.T(:,i-1)+273.15) + -X4.acb ./ (X.T(:,i-1) + 273.15) + X4.acc)) .*...
    (2 .* X4.aca .* (X.T(:,i-1)+273.15) ./ ((X.T(:,i-1)+273.15) .^ 4) - X4.acb ./ ((X.T(:,i-1) + 273.15) .^ 2));

thaos = ((X2.theta_sat(1)-thetas(1)).^(7/3))./((X2.theta_sat(1)).^2);
thaols= ((thetas(1)).^(7/3))./((X2.theta_sat(1)).^2);
Dilo= 1e-7.*exp(-577./(X.T(:,i-1)+273.15-145))/X4.Dil_coe; %1.026 for 18O;1.013
Dil = Dilo.*thaol.*X.theta(:,i-1);
Dil0=10^(-7)*exp(-577/(Ts(1)+273.15-145))/X4.Dil_coe*thaols(1)*thetas(1);


nD = 0.67+0.33.*exp(1-X.theta(1,i-1)./X2.theta_res(1,1));
nDs = 0.67+0.33.*exp(1-thetas(1,1)./X2.theta_res(1,1));
Div = Dv0.*10.^5./101325.*((X.T(:,i-1)+273.15)./273.15).^1.88.*(X4.alphadiff).^nD.*(X2.theta_sat-X.theta(:,i-1)).*thao;
Div0 = Dv0.*10.^5./101325.*((Ts(1)+273.15)./273.15).^1.88.*(X4.alphadiff).^nDs.*(X2.theta_sat(1)-thetas(1)).*thaos;
cv= cvsat.* hr;

nk = ((thetas-X2.theta_res(1))*0.5+(X2.theta_sat(1)-thetas)*1)/(X2.theta_sat(1)-X2.theta_res(1));
alphak=(1/X4.alphadiff)^nk;

if X.pTest~=5&&X.pTest~=6
    Dil=0;
    Dil0=0;
end

X1=struct('k',k,'S', S,'C', C,'theta', X.theta,'cvsat', cvsat,'dcvsatdT', dcvsatdT,'hr', hr,...
    'dhrdh',dhrdh,'dhrdT', dhrdT,'theta_sat', X2.theta_sat,'theta_res',X2.theta_res,...
    'lamdaE',lamdaE,'Csoil', Csoil,'CW', CW,'CV',CV,'Tcoe', Tcoe,'kH',kH,'rho', X2.rho, ...
    'mass_transfer_coefficient',X2.mass_transfer_coefficient,'coe', X2.coe, 'f',X2.f,'kad', X2.kad,...
    'rhob',X2.rhob, ...
    'alpha',alpha,'alphas',alphas,'dalphadT',dalphadT,'Dil',Dil,'Dil0',Dil0,'Div',Div,...
    'Div0',Div0,'cv',cv,'nD',nD,'alphak',alphak,'thao',thao,'qls',qls,...
    'qvs',qvs,'G',G,'alphadiff',X4.alphadiff,'qevap',qevap,'cvs',cvs,'Ts',Ts,'hrs',hrs,'cvsats',cvsats,'rbw',rbw);

end

