function [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] =TH_function_ke(Tinitial,hrinitial,X3)

Jmatrix=zeros(2,2);

x0=[Tinitial;hrinitial];

of=zeros(2,1);
iterationtimes=1;
maxiteration=100;

errtarget=10^(-13);
fllag=true;
i = X3.i;
i = i/i+1;
while fllag

    [rbh,rbw]=AeroRes(x0(1),X3.Ta,X3.u,X3.Cp,X3.root_flag,X3.rs(1));

    Jmatrix(1,1) = (1800793563611553*X3.hr(1)*X3.s(1))/18014398509481984 + (41570*log(x0(2))*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2))/(441*X3.dz) + (2*X3.Dv(1)*X3.hr(1)*X3.s(1))/X3.dz -...
        (X3.R*X3.k(1)*X3.ksat(1)*X3.lamda(1)*X3.yeeta(1)*log(x0(2)))/(2*X3.M*X3.g*X3.he(1)*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2)*((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1) + 1)) +...
        (X3.R*X3.k(1)*X3.ksat(1)*X3.lamda(1)*X3.yeeta(1)*log(x0(2))*(X3.h(1) - (20785*log(x0(2))*(x0(1) + 5463/20))/441))/(X3.M*X3.dz*X3.g*X3.he(1)*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2)*((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1) + 1));

    Jmatrix(1,2) = (1800793563611553*X3.cvsat(1))/18014398509481984 + (2*X3.Dv(1)*X3.cvsat(1))/X3.dz(1) + (41570*(x0(1) + 5463/20)*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2))/(441*X3.dz*x0(2)) -...
        (X3.R*X3.k(1)*X3.ksat(1)*X3.lamda(1)*X3.yeeta(1)*(x0(1) + 5463/20))/(2*X3.M*X3.g*X3.he(1)*x0(2)*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2)*((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1) + 1)) +...
        (X3.R*X3.k(1)*X3.ksat(1)*X3.lamda(1)*X3.yeeta(1)*(x0(1) + 5463/20)*(X3.h(1) - (20785*log(x0(2))*(x0(1) + 5463/20))/441))/(X3.M*X3.dz*X3.g*X3.he(1)*x0(2)*((X3.k(1)*X3.ksat(1))/((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1)))^(1/2)*((X3.R*log(x0(2))*(x0(1) + 5463/20))/(X3.M*X3.g*X3.he(1)))^(X3.lamda(1)*X3.yeeta(1) + 1));

    Jmatrix(2,1)=(1000000*X3.hr(1)*X3.rho*X3.s(1)*((1361023836188383*x0(1))/576460752303423488-2501/1000))/rbw(1)-(2*X3.kH(1))/X3.dz-(21265997440443484375*X3.rho*(X3.cva-X3.cvsat(1)*x0(2)+X3.hr(1)*X3.s(1)*(X3.T(1)-x0(1))))/(9007199254740992*rbw(1))-X3.Cp/rbh;
    Jmatrix(2,2)=(1000000*X3.cvsat(1)*X3.rho*((1361023836188383*x0(1))/576460752303423488 - 2501/1000))/rbw(1);

    if isnan(Jmatrix(1,1))
        Jmatrix(1,1)=(X3.hr(1)*X3.s(1))/rbw + (41570*log(x0(2))*((3423715700080493*X3.ksat(1))/295147905179352825856)^(1/2))/(441*X3.dz) + (2*X3.Dv(1)*X3.hr(1)*X3.s(1))/X3.dz;
        Jmatrix(1,2)=X3.cvsat(1)/rbw + (2*X3.Dv(1)*X3.cvsat(1))/X3.dz + (41570*((3423715700080493*X3.ksat(1))/295147905179352825856)^(1/2)*(x0(1) + 5463/20))/(441*X3.dz*x0(2));
    end



    ks=(log(x0(2))*X3.R*(x0(1)+273.15) /(X3.M*X3.g)/(X3.he(1)))^(-X3.lamda(1)*X3.yeeta(1))*X3.ksat(1);
    if isnan(ks)||ks>X3.ksat(1)
        ks = X3.ksat(1);
    end
    k=(ks*X3.k(1))^0.5;
    of(1)=1/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1)-x0(1)))-X3.cva)-...
        X3.Dv(1)/(X3.dz/2)*X3.cvsat(1)*(X3.hr(1)-x0(2))-X3.Dv(1)/(X3.dz/2)*X3.s(1)*X3.hr(1)*(X3.T(1)+0-(x0(1)+0))-...
        (X3.k(1)*(X3.h(1)-log(x0(2))*X3.R*(x0(1)+273.15) /(X3.M*X3.g))/(X3.dz/2)-X3.k(1));

    of(1)=1/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1)-x0(1)))-X3.cva)-...
        X3.Dv(1)/(X3.dz/2)*X3.cvsat(1)*(X3.hr(1)-x0(2))-X3.Dv(1)/(X3.dz/2)*X3.s(1)*X3.hr(1)*(X3.T(1)+0-(x0(1)+0))-...
        (k(1)*(X3.h(1)-log(x0(2))*X3.R*(x0(1)+273.15) /(X3.M*X3.g))/(X3.dz/2)-k(1));


    of(2)=X3.Rnet-X3.Cp/rbh*((x0(1)+273.15)-(X3.Ta+273.15))-...
        X3.rho*(2.501-2.361*10^(-3)*x0(1))*10^6/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1)-x0(1)))-X3.cva)+...
        X3.kH(1)/(X3.dz/2)*((X3.T(1)+273.15)-(x0(1)+273.15));




    TTHH=Jmatrix\(-of)+x0;
    if isnan(TTHH)
        error('TTHH is nan ')
    end
    err1=abs(TTHH(1)-x0(1));
    err2=abs(TTHH(2)-x0(2));

    x0=TTHH;

    iterationtimes=iterationtimes+1;
    if err1<errtarget&&err2<errtarget
        fllag=false;
    end
    if iterationtimes>maxiteration
        fllag=false;
    end
    if x0(2,1)>1

        x0(2,1)=0.999999999;
    end
    if ~isreal(x0)
        error("complex upper boundary")
    end
end
qevap=1/rbw(1)*((x0(2)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1)-x0(1)))-X3.cva);

T11=TTHH(1,1);
H22=TTHH(2,1);

[rbh,rbw]=AeroRes(T11,X3.Ta,X3.u,X3.Cp,X3.root_flag,X3.rs(1));
lamtaE_s=(2.501-2.361*10^(-3)*T11)*10^6;

Sensibleheat=X3.Cp/rbh*((T11+273.15)-(X3.Ta+273.15));
latentheat=X3.rho*lamtaE_s*qevap;
G=-X3.kH(1,1)/(X3.dz/2)*((X3.T(1,1)+273.15)-(T11+273.15));

energy_balance=X3.Rnet-(Sensibleheat+latentheat+G);

if abs(energy_balance)>20
    error("surface energy is not balanced in isofun_hrs");
end

end
