function D_BA2=BA2_fun(qevap)
global out_X1 out_X2 out_X4 out_X3 out_X
Div = out_X1.Div;
cvsat = out_X1.cvsat;
Dil = out_X1.Dil;
hr = out_X1.hr;
alphak = out_X1.alphak;
alpha = out_X1.alpha;
alphas = out_X1.alphas;
dz = out_X.dz;
cils = out_X1.cils;
ref = out_X4.ref;
n = out_X.n;
Ta = out_X.Ta;
M = out_X3.M;
R = out_X3.R;
rho = out_X3.rho;
hrs = out_X1.hrs;
Ts = out_X1.Ts;
T = out_X.T;


cvsata=610.78.*exp(17.27.*Ta./(Ta+237.3)).*M./R./rho./(Ta+273.15);
cvsats=610.78.*exp(17.27.*Ts(1)./(Ts(1)+237.3)).*M./R./rho./(Ts(1)+273.15);
zrho=1e-3.*exp(31.3716-6014.79./(T+273.15)-7.92495.*1e-3.*(T+273.15))./(T+273.15);
B2=zeros(n,1);

qevap_s=qevap;

alphak_(1:n,1)=alphak; % 10166 for 2H and 1.0324 for 18O when nk=1. THese val;ues fropm TT zhou.
alphaks=alphak;%1.0166;
zv=Div.*cvsat./qevap_s;
zl=Dil/qevap_s;

B1=hr.*zv./(zl+hr.*zv).*(alphak_-alpha);
B2(1,1)=(log(hr(1).*cvsat(1).*1000.*(alphak_(1)-alpha(1)))-log(hrs(1)*cvsats(1).*1000.*(alphaks(1)-alphas(1))))/(dz/2);
B2(2:end,1)=(log(cvsat(2:end).*hr(2:end).*1000.*(alphak_(2:end)-alpha(2:end)))-log(cvsat(1:end-1).*hr(1:end-1).*1000.*(alphak_(1:end-1)-alpha(1:end-1))))/dz;

B=B1.*B2;
Ds=((cils/1000/0.02*0.018)-ref)/ref*1000/1000;
D(1,1)=(B(1)+Ds/(dz/2))/(1/(zl(1)+zv(1)*hr(1))+1/(dz/2));

for i=2:n
    D(i,1)=(B(i)+D(i-1)/(dz))/(1/(zl(i)+zv(i)*hr(i))+1/(dz));
end
D_BA2=D.*1000;
end