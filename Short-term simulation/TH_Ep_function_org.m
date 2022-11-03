function [TH,fval_Ep,qevap,Sensibleheat,latentheat,G,rbh,rbw,it] = TH_Ep_function_org(Tinitial,hrinitial,X3)


qevap=Ep;

T1=Tinitial;
x0=T1;

while true

    [rbh,rbw]=AeroRes(x0,rs(1));
    T11=(Rnet+Cp/rbh*(Ta+273.15)+...
        kH(1)/(dz(1)/2)*(T(1)+273.15)-...
        Cp/rbh*273.15-rho*10^6*qevap*2.501-kH(1)/(dz(1)/2)*273.15)/(Cp/rbh-rho*10^6*qevap*2.361*10^(-3)+kH(1)/(dz(1)/2));

    hs(1)=(-qevap-k(1))/k(1)*dz/2+h(1);
    hrs(1,1)=exp(hs(1)*M*g/(x0+273.15)/R);




    roa=0.001*exp(31.3716-6014.79/(Ta+273.15)-0.00792495*(Ta+273.15))/(Ta+273.15)*hra;

    x01=(Rnet+Cp/rbh*Ta-...
        1e6*rho*2.501*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*M*g/(x0+273.15)/R))-roa)/(rbh+rbw)/1000)+...
        kH(1)/(dz/2)*T(1))/...
        (Cp/rbh-rho*1e6*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*M*g/(x0+273.15)/R))-roa)/(rbh+rbw)/1000)*2.631*1e-3+kH(1)/(dz/2));



    rov=0.001*exp(31.3716-6014.79/(x01+273.15)-0.00792495*(x01+273.15))/(x01+273.15)*(exp(hs(1)*M*g/(x0+273.15)/R));
    qevap1=(rov-roa)/(rbh(1)+rbw(1))/1000;
    if Ep==0
        qevap1=0;
        x01=(Rnet+Cp/rbh*Ta-...
            1e6*rho*2.501*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*M*g/(x0+273.15)/R))-roa)*0/(rbh+rbw)/1000)+...
            kH(1)/(dz/2)*T(1))/...
            (Cp/rbh-rho*1e6*((0.001*exp(31.3716-6014.79/(x0+273.15)-0.00792495*(x0+273.15))/(x0+273.15)*(exp(hs(1)*M*g/(x0+273.15)/R))-roa)*0/(rbh+rbw)/1000)*2.631*1e-3+kH(1)/(dz/2));

    end



    if abs(T11-T1)<1e-1&&abs(x01-x0)<1e-1
        break
    else
        T1=T11;
        x0=x01;
    end
end


T11=x01; %1220
qevap=qevap1; %1220


TTHH(1)=T11;
TTHH(2)=hrs(1,1);hr(1);%hr1;
of=1;%f22;

lamtaE_s=(2.501-2.361*10^(-3)*T11)*10^6;

Sensibleheat=Cp/rbh*((T11+273.15)-(Ta+273.15));%
latentheat=rho*lamtaE_s*qevap;
G=-kH(1,1)/(dz/2)*((T(1,1)+273.15)-(T11+273.15)); %heat flux into soil
err=Rnet-(Sensibleheat(1)+latentheat(1)+G(1));


err=Rnet-(Sensibleheat(1)+latentheat(1)+G(1));

qevap=Ep;
hs_err=1;
hs=h(1);
it=0;
flag=true;
while flag
    thetas=(hs(1)/he(1))^(-lamta(1))*(thetasat(1)-thetares(1))+thetares(1);
    ks(1)=ksat(1)*((thetas(1)-thetares(1))/(thetasat(1)-thetares(1)))^(yeeta(1));
    ke=(ks(1)+k(1))/2;
    h_s(1)=(-qevap-(ke))/(ke)*dz/2+h(1);%hr(1);%(hra+hr(1))/2;%hr(1);

    if abs(h_s-hs)<1e-5
        thetas=(hs(1)/he(1))^(-lamta(1))*(thetasat(1)-thetares(1))+thetares(1);
        ks(1)=ksat(1)*((thetas(1)-thetares(1))/(thetasat(1)-thetares(1)))^(yeeta(1));
        qevap2=-(h_s(1)-h(1))/(dz/2)*(ks(1)+k(1))/2-(ks(1)+k(1))/2;
        flag=false;
    else
        hs=h_s;
    end

    it=it+1;
end

hrs(1,1)=exp(hs(1)*M*g/(T11+273.15)/R);
TTHH(2)=hrs(1,1);
end

