function [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_Ep_function2(Tinitial,hrinitial,X3)

qevap=X3.Ep;%levenberg-marquardt,'StepTolerance',1e-10
i = X3.i;
i = i/i+1;
T1=Tinitial;
x0=T1;

while true

    [rbh,rbw]=AeroRes(x0,X3.Ta,X3.u,X3.Cp,X3.root_flag,X3.rs(1));

    T11=(X3.Rnet+X3.Cp/rbh*(X3.Ta+273.15)+...
        X3.kH(1)/(X3.dz/2)*(X3.T(1,i-1)+273.15)-...
        X3.Cp/rbh*273.15-X3.rho*10^6*qevap*2.501-X3.kH(1)/(X3.dz/2)*273.15)/...
        (X3.Cp/rbh-X3.rho*10^6*qevap*2.361*10^(-3)+X3.kH(1)/(X3.dz/2));

    hs(1)=(-qevap-X3.k(1))/X3.k(1)*X3.dz/2+X3.h(1);%hr(1);%(hra+hr(1))/2;%hr(1);
    hrs(1,1)=exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R);%log(hrs(1,j))*R*(Ts(1,j)+273.15)/g/M;

    roa=0.001*exp(31.3716-6014.79/(X3.Ta+273.15)-0.00792495*(X3.Ta+273.15))/(X3.Ta+273.15)*X3.hra;

    x01=(X3.Rnet+X3.Cp/rbh*X3.Ta-...
        1e6*X3.rho*2.501*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R))-roa)/(rbh+rbw)/1000)+...
        X3.kH(1)/(X3.dz/2)*X3.T(1,i-1))/...
        (X3.Cp/rbh-X3.rho*1e6*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R))-roa)/(rbh+rbw)/1000)*2.631*1e-3+X3.kH(1)/(X3.dz/2));

    rov=0.001*exp(31.3716-6014.79/(x01+273.15)-0.00792495*(x01+273.15))/(x01+273.15)*(exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R));
    qevap1=(rov-roa)/(rbh(1)+rbw(1))/1000;
    if X3.Ep==0
        qevap1=0;
        x01=(X3.Rnet+X3.Cp/rbh*X3.Ta-...
            1e6*X3.rho*2.501*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R))-roa)*0/(rbh+rbw)/1000)+...
            kH(1)/(dz/2)*T(1))/...
            (X3.Cp/rbh-X3.rho*1e6*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*X3.M*X3.g/(x0+273.15)/X3.R))-roa)*0/(rbh+rbw)/1000)*2.631*1e-3+X3.kH(1)/(X3.dz/2));
    end

    if abs(T11-T1)<1e-1&&abs(x01-x0)<1e-1
        break
    else
        T1=T11;
        x0=x01;
    end
end
T11=x01;
qevap=qevap1;



TTHH(1)=T11;
TTHH(2)=hrs(1,1);%hr1;
of=1;%f22;

lamtaE_s=(2.501-2.361*10^(-3)*T11)*10^6;

Sensibleheat=X3.Cp/rbh*((T11+273.15)-(X3.Ta+273.15));%
latentheat=X3.rho*lamtaE_s*qevap;
G=-X3.kH(1,1)/(X3.dz/2)*((X3.T(1,i-1)+273.15)-(T11+273.15)); %heat flux into soil
err=X3.Rnet-(Sensibleheat(1)+latentheat(1)+G(1));

if abs(err)>10
    if ~err<0
        Sensibleheat=Sensibleheat-Sensibleheat/(Sensibleheat(1)+latentheat(1)+G(1))*err;
        latentheat=latentheat-latentheat/(Sensibleheat(1)+latentheat(1)+G(1))*err;
        G=G-G/(Sensibleheat(1)+latentheat(1)+G(1))*err;
    else
        Sensibleheat=Sensibleheat+Sensibleheat/(Sensibleheat(1)+latentheat(1)+G(1))*err;
        latentheat=latentheat+latentheat/(Sensibleheat(1)+latentheat(1)+G(1))*err;
        G=G+G/(Sensibleheat(1)+latentheat(1)+G(1))*err;
    end

end
err=X3.Rnet-(Sensibleheat(1)+latentheat(1)+G(1));
qevap=X3.Ep;
hs_err=1;
hs=X3.h(1);
it=0;
flag=true;
while flag
    thetas=(hs(1)/X3.he(1))^(-X3.lamda(1))*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1);
    theta1=(X3.h(1)/X3.he(1))^(-X3.lamda(1))*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1);

    ks(1)=X3.ksat(1)*((thetas(1)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
    ke=(ks(1)+X3.k(1))/2; %BA2
    h_s(1)=(-qevap-(ke))/(ke)*X3.dz/2+X3.h(1);

    if abs(h_s-hs)<1e-5
        thetas=(hs(1)/X3.he(1))^(-X3.lamda(1))*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1);
        theta1=(X3.h(1)/X3.he(1))^(-X3.lamda(1))*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1);
        ks(1)=X3.ksat(1)*((thetas(1)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=(ks(1)+X3.k(1))/2;%BA2
        flag=false;
    else
        hs=h_s;
    end

    it=it+1;

end
hrs(1,1)=exp(hs(1)*X3.M*X3.g/(T11+273.15)/X3.R);%log(hrs(1,j))*R*(Ts(1,j)+273.15)/g/M;
TTHH(2)=hrs(1,1);
end




