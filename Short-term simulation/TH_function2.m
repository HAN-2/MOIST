function [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_function2(Tinitial,hrinitial,X3,hrcria)
i=X3.i;

Jmatrix=zeros(2,2);

x0=[Tinitial;hrinitial];

of=zeros(2,1);
iterationtimes=1;
maxiteration=100;

errtarget=10^(-13);
fllag=true;

Ss=(hs(1)/X3.he(1))^(-X3.lamta(1));
ks=X3.ksat(1)*(Ss)^(X3.yeeta(1));
ke=(X3.k(1)+ks)/2;

while false
    [rbh,rbw]=AeroRes(x0(1),X3.rs);

    Jmatrix(1,1)=(X3.hr(1)*X3.s(1))/rbw(1)+(2*X3.Dv(1)*X3.hr(1)*X3.s(1))/X3.dz;
    Jmatrix(1,2)=X3.cvsat(1)/rbw(1)+(2*X3.Dv(1)*X3.cvsat(1))/X3.dz;

    %-----------------------------------------------------------------------------------------------------------


    Jmatrix(2,1)=(1000000*X3.hr(1)*X3.rho*X3.s(1)*((1361023836188383*x0(1))/576460752303423488-2501/1000))/rbw(1)-(2*X3.kH(1))/X3.dz-...
        (21265997440443484375*X3.rho*(X3.cva-X3.cvsat(1)*x0(2)+X3.hr(1)*X3.s(1)*(X3.T(1,i-1)-x0(1))))/(9007199254740992*rbw(1))-X3.Cp/rbh;
    Jmatrix(2,2)=(1000000*X3.cvsat(1)*X3.rho*((1361023836188383*x0(1))/576460752303423488 - 2501/1000))/rbw(1);

    of(1)=1/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1)*X3.hr(1)*(X3.T(1,i-1)-x0(1)))-X3.cva)...
        -X3.Dv(1)/(X3.dz/2)*X3.cvsat(1)*(X3.hr(1)-x0(2))-X3.Dv(1)/(X3.dz/2)*X3.s(1)*X3.hr(1)*(X3.T(1,i-1)+0-(x0(1)+0))...
        -(-(-ke*(X3.h(1)-hs(1))/(X3.dz/2)+ke));

    of(2)=X3.Rnet...
        -X3.Cp/rbh*((x0(1)+273.15)-(X3.Ta+273.15))-...
        X3.rho*(2.501-2.361*10^(-3)*x0(1))*10^6/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1,i-1)-x0(1)))-X3.cva)...
        +X3.kH(1)/(X3.dz/2)*((X3.T(1,i-1)+273.15)-(x0(1)+273.15));


    TTHH=Jmatrix\(-of)+x0;
    err1=abs(TTHH(1)-x0(1));
    err2=abs(TTHH(2)-x0(2));

    x0=TTHH;

    iterationtimes=iterationtimes+1;
    if err1<errtarget&&err2<errtarget
        fllag=false;
    end
    if iterationtimes>maxiteration/2&&X3.root_flag&&iterationtimes<maxiteration/2+2
        if X3.Rnet==0
            X3.Rnet=-15;
        else
            X3.Rnet=0;
        end
    end
    if iterationtimes>maxiteration
        fllag=false;
    end
    if x0(2,1)>1

        x0(2,1)=hr(1);
    end
end
[rbh,rbw]=AeroRes(X3.Ta,X3.rs);
if isempty(hrcria)
    hrcria=exp(hs(1)*X3.M*X3.g/X3.R/(X3.T(1,i-1)+273.15));
end
qvs_hr=-X3.Dv(1)/(X3.dz(1)/2)*X3.cvsat(1)*(X3.hr(1)-hrcria);
qls=-(-ke*(X3.h(1)-hs(1))/(X3.dz/2)+ke);
qvs_T=X3.Dv(1)/(X3.dz/2)*X3.s(1)*X3.hr(1);



% syms x01
y1=@(x01) X3.Rnet-X3.Cp/X3.rbh*(x01+273.15)+X3.Cp/rbh*(X3.Ta+273.15)-X3.rho*10^6*(-(qvs_hr-qvs_T*(X3.T(1,i-1)-(x01))-(X3.qls)))*2.501+X3.rho*10^6*(-(qvs_hr-qvs_T*(X3.T(1,i-1)-(x01))-(X3.qls)))*(2.361*10^(-3)*x01)+X3.kH(1)/(X3.dz(1)/2)*(X3.T(1,i-1))-X3.kH(1)/(X3.dz(1)/2)*(x01);

sol1=-((0.00021177467174925878*(-1*(9444*X3.dz*qvs_T*rbh(1)*X3.rho*(X3.Cp*X3.dz*X3.Ta-2.501*10^6*X3.dz*X3.qls*rbh(1)*X3.rho+2.501*10^6*X3.dz*qvs_hr*rbh(1)*X3.rho-2.501*10^6*X3.dz*qvs_T*rbh(1)*X3.rho*X3.T(1,i-1)+X3.dz*rbh(1)*X3.Rnet+...
    2*X3.kH(1)*rbh(1)*X3.T(1,i-1))+(-1*X3.Cp*X3.dz+2361*X3.dz*X3.qls*rbh(1)*X3.rho-2361*X3.dz*qvs_hr*rbh(1)*X3.rho+2361*X3.dz*qvs_T*rbh(1)*X3.rho*X3.T(1,i-1)+2.501*10^6*X3.dz*qvs_T*rbh(1)*X3.rho-...
    2*X3.kH(1)*rbh(1))^2)^0.5+X3.Cp*X3.dz-2361*X3.dz*X3.qls*rbh(1)*X3.rho+2361*X3.dz*qvs_hr*rbh(1)*X3.rho-2361*X3.dz*qvs_T*rbh(1)*X3.rho*X3.T(1,i-1)-2.501*10^6*X3.dz*qvs_T*rbh(1)*X3.rho+...
    2*X3.kH(1)*rbh(1)))/(X3.dz*qvs_T*rbh(1)*X3.rho));

qevap=-(qvs_hr-qvs_T*(X3.T(1,i-1)-(sol1))-(X3.qls));
T11=sol1;

lamtaE_s=(2.501-2.361*10^(-3)*T11)*10^6;
Sensibleheat=X3.Cp/rbh*((T11+273.15)-(X3.Ta+273.15));
latentheat=X3.rho*lamtaE_s*qevap;
G=-X3.kH(1,1)/(X3.dz/2)*((X3.T(1,i-1)+273.15)-(T11+273.15));
TTHH(1,1)=T11;
TTHH(2,1)=hrcria;

end


