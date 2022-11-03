function cils = cils_function_test(i,X,X1,X3,X4)
global out_X1
i=i/i+1;
qvs_rc =1;
if X.iso_spe

    catmos=(-112/1000+1)*X4.ref*1000*0.02/0.018;  %atmosphere isotopic composition

else

    catmos=(-15/1000+1)*X4.ref*1000*0.02/0.018;  %atmosphere isotopic composition

end


civa=X3.cva*catmos;%*alpha(1);%alpha(Ta); cilb means the concentration of Ratmos
cvsats=610.78.*exp(17.27*X1.Ts./(X1.Ts+237.3))*X3.M./X3.R/X3.rho./(X1.Ts+273.15);
cvs=(X1.hrs*X1.cvsat(1,1)-X3.s(1,1)*X1.hr(1)*(X.T(1,i-1)-X1.Ts));

alphas= exp(-(X4.aca./(X1.Ts+273.15).^2+X4.acb./(X1.Ts+273.15)+X4.acc));%1;%
if X.pTest==1||X.pTest==2
    alphas=1;
end

rbw_re(1)=((cvs(1)-X3.cva)/X1.qevap); %1212
if X1.qevap==0
    rbw_re(1)=abs((cvs(1)-X3.cva)/1e-20); %1212
end
if X.qprec==0%abs(qls+qvs-(-qevap))<10^(-15)%Rain_lost_flag

    cils=(X1.Div(1)/(X.dz/2)*X1.hr(1)*X1.cvsat(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1))/...
        (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(1)+alphas*cvs(1)*X1.Div(1)/(X.dz/2)+X1.Dil(1)/(X.dz/2));
else

    if X.qprec>X1.qevap%theta(1)==thetasat(1)
        cils=cil_rain;
    else
        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(0)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(1-0)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
    end
end
if X.frozen_flag
    cils=X.cil(1,i-1);
end
if cils<0%||cils>1.8
    error('cils is smaller than 0')
end
out_X1.cils = cils;
end



