function [qls,qvs,qTs,...
    Ts,hrs,cvs,thetas,...
    G,rbh,rbw,...
    qevap,Ep_flag]=surface_function(X3,i,X2,X)


i=i/i+1;
cvsat = X3.cvsat;
s = X3.s;
hr= X3.hr;
Tinitial=X3.T(1,i-1);
hrinitial=X3.hr(1,1);
hs_flag=true;
if X3.itm
    if hs_flag
        Ep_flag=1;
        [TH,fval_Ep,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_Ep_function2(Tinitial,hrinitial,X3);
        Ts(1,1)=TH(1,1);
        hrs(1,1)=TH(1,2);
        hs(1,1)=log(hrs(1,1))*X3.R*(Ts(1,1)+273.15)/X3.g/X3.M;
        thetas(1)=((hs(1)/X3.he(1))^-X3.lamta(1)*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1));
        ks(1)=X3.ksat(1)*((thetas(1)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=X3.ksat(1)*((((thetas(1)*X3.theta(1,i-1))^0.5)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=(ks(1)+X3.k(1))/2;

        if hrs(1)>hra
            hrcria=hrs(1);
        end
    end

    if  hs(1)<=X3.hcria

        Ep_flag=0;
        hs(1)=hcria;
        [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_function2(Tinitial,hrinitial,X3,hrcria);%isofun_hrs2(Tinitial,hrinitial,Ep,rs(1));
        Ts(1,gj)=TTHH(1,1);
        hrs(1,gj)=TTHH(2,1);
        thetas(1)=((hs(1)/X3.he(1))^-X3.lamta(1)*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1));
        ks(1)=X3.ksat(1)*((thetas(1)-X3.thetares(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=X3.ksat(1)*((((thetas(1)*X3.theta(1,i-1))^0.5)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=(X3.k(1)+ks(1))/2;

        hs_flag=false;
    end

    if BA1
        [TH,fval_Ep,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_Ep_function2(Tinitial,hrinitial,X3);
        Ts(1,gj)=TH(1,1);
        hrs(1,gj)=TH(1,2);
        hs(1,1)=log(hrs(1,1))*X3.R*(Ts(1,1)+273.15)/X3.g/X3.M;
        thetas(1)=((hs(1)/X3.he(1))^-X3.lamta(1)*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1));
        ks(1)=X3.ksat(1)*((thetas(1)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        ke=X3.ksat(1)*((((thetas(1)*X3.theta(1,i-1))^0.5)-X3.theta_res(1))/(X3.theta_sat(1)-X3.theta_res(1)))^(X3.yeeta(1));
        qevap=Ep;
        Ep_flag=1;
    end
else
    if ~X3.frozen_flag
        hs_flag=true;
        if hs_flag
            Ep_flag=1;
            [TH,fval_Ep,qevap,Sensibleheat,latentheat,G,rbh,rbw,it,fail_flag] = TH_Ep_function(Tinitial,hrinitial,X3);

            Ts=TH(1,1);
            hrs=TH(1,2);
            hs=log(hrs)*X3.R*(Ts+273.15)/X3.g/X3.M;
            thetas=((hs/X3.he(1))^-X3.lamda(1)*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1));
            ks=X3.ksat(1)*((thetas-X3.theta_res(1,1))/(X3.theta_sat(1,1)-X3.theta_res(1,1))).^(X3.yeeta(1,1));
            ke=X3.k(1);


            if hrs(1)>X3.hra
                hrcria=hrs;
            end

            if it>100
                hs=X3.hcria;
            else
                energy_balance=X3.Rnet-(Sensibleheat(1)+latentheat(1)+G(1));
                if abs(energy_balance)>20
                    error("surface energy is not balanced--line 86");
                end
            end
        end

        if  hs<=X3.hcria || fail_flag == 0
            Ep_flag=0;
            isofun_hrs3_flag=0;
            if hs<=X3.hcria
                hs=X3.hcria;
            end
            if  0
                [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_Ep_function(Tinitial,hrinitial,X3);%isofun_hrs2(Tinitial,hrinitial,Ep,rs(1));
                Ts=TTHH(1,1);
                hrs=TTHH(1,2);
            else

                [TTHH,of,qevap,Sensibleheat,latentheat,G,rbh,rbw] = TH_function(Tinitial,hrinitial,X3);
                Ts=TTHH(1,1);
                hrs=TTHH(2,1);
            end

            thetas=((hs/X3.he(1,1))^-X3.lamda(1,1)*(X3.theta_sat(1,1)-X3.theta_res(1,1))+X3.theta_res(1,1));
            ks=X3.ksat(1,1)*((thetas-X3.theta_res(1,1))/(X3.theta_sat(1,1)-X3.theta_res(1,1)))^(X3.yeeta(1,1));
            ke=X3.k(1);
            hs_flag=false;
            energy_balance=X3.Rnet-(Sensibleheat(1)+latentheat(1)+G(1));

            phaie = X2.ksat.*X2.he./(1-X2.lamda.*X2.yeeta);
            S = (X.theta-X2.theta_res)./(X2.theta_sat-X2.theta_res);

        end

    else
        Ep_flag=0;
        [TH,fval_Ep,qevap,Sensibleheat,latentheat,G,rbh,rbw,it] = TH_frozen_function(Tinitial,hrinitial,X3);
        Ts=TH(1,1);
        hrs=TH(1,2);
        hs=log(hrs)*X3.R*(Ts+273.15)/X3.g/X3.M;
        thetas=((hs/X3.he(1))^-X3.lamda(1)*(X3.theta_sat(1)-X3.theta_res(1))+X3.theta_res(1));
    end

    if ~X3.frozen_flag
        j=1;
        gj=1;

        cvs=(hrs(1)*X3.cvsat(1,gj)-X3.s(1,gj)*X3.hr(1)*(X3.T(1,i-1)-Ts(1)));
        qvs=-X3.Dv(1,j)/(X3.dz/2)*(X3.cvsat(1,j)*X3.hr(1,j)-cvs(1,j));
        qls=-(ke*(X3.h(1,i-1)-hs)/(X3.dz/2)-ke);
        qTs=G;
    else
        cvs=(hrs(1)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1,1)-Ts(1)));
        qls=0;
        qvs=0;
        qevap=0;
        qTs=G;
    end
    ql0=(X3.qprec-qevap(1,1));%0401
    QEerr=qls+qvs-(-qevap); %1112
    qls=qls-QEerr;
    qls=qls+X3.qprec;%1212
    QEerr=qls+qvs-(ql0); %1112

    if abs(QEerr)>1e-12
        qls
        qvs
        qevap
        error('qls, qvs and qevap unbalanced')
    end

end

