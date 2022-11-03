function [cils,ws] = cils_function(i,X,X1,X3,X4)
global out_X1 out_X4 month
i=i/i+1;
qvs_rc =1;
if X.iso_spe
    if ~X4.itm
        alpha_atmos=exp(-(X4.aca./(X.Ta+273.15).^2+X4.acb./(X.Ta+273.15)+X4.acc));
        e_=(1/alpha_atmos-1)*1000;
        qprec_amount=X.qprec_Day*24*3600*1000;
        [cil_rain_qianfenzhi,cil_rain]=IVA_H2(qprec_amount,'H',X);
    end

    if X.qprec~=0
        catmos_qianfenzhi=(cil_rain_qianfenzhi-e_)/(1/alpha_atmos);
        catmos=(catmos_qianfenzhi/1000+1)*X4.ref*1000*0.02/0.018;
    else
        catmos=(-130/1000+1)*X4.ref*1000*0.02/0.018;

    end
else
    if ~X4.itm
        alpha_atmos=exp(-(X4.aca./(X.Ta+273.15).^2+X4.acb./(X.Ta+273.15)+X4.acc));
        e_=(1/alpha_atmos-1)*1000;
        qprec_amount=X.qprec_Day*24*3600*1000;
        [cil_rain_qianfenzhi,cil_rain]=IVA_H2(qprec_amount,'O',X);
    end

    if X.qprec~=0
        catmos_qianfenzhi=(cil_rain_qianfenzhi-e_)/(1/alpha_atmos);
        catmos=(catmos_qianfenzhi/1000+1)*X4.ref*1000*0.02/0.018;
    else
        catmos=(-20/1000+1)*X4.ref*1000*0.02/0.018;
    end
end


civa=X3.cva*catmos;
cvsats=610.78.*exp(17.27*X1.Ts./(X1.Ts+237.3))*X3.M./X3.R/X3.rho./(X1.Ts+273.15);
cvs=(X1.hrs*X1.cvsat(1,1)-X3.s(1,1)*X1.hr(1)*(X.T(1,i-1)-X1.Ts));

alphas= exp(-(X4.aca./(X1.Ts+273.15).^2+X4.acb./(X1.Ts+273.15)+X4.acc));%1;%
if X.pTest==1||X.pTest==2
    alphas=1;
end

rbw_re(1)=abs((cvs(1)-X3.cva)/X1.qevap); %1212
rbw_re(1)=((cvs(1)-X3.cva)/X1.qevap); %1212

if X1.qevap==0
    rbw_re(1)=abs((cvs(1)-X3.cva)/1e-20); %1212
    rbw_re(1)=((cvs(1)-X3.cva)/1e-20); %1212
end
if out_X1.ql(1)<0
    ws = 0;
else
    ws = 1;
end

if out_X1.ql(2)<0
    w1 = 0;
else
    w1 = 1;
end


if X.qprec==0

    cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
        (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
    D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    if abs(D_cils)>200
        ws = 1-ws;
        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    end
    if abs(D_cils)>200 && rbw_re <0
        rbw_re = abs(rbw_re);
        cils1=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        ws = 1-ws;
        cils2=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        cils = min(cils1,cils2);
    end

else

    cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
        (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
    D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    if abs(D_cils)>200
        ws = 1-ws;
        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    end
    if abs(D_cils)>200
        cils=cil_rain;
    end


end
if X.frozen_flag
    cils=X.cil(1,i-1);
end

if X.frozen_flag
    if X.qprec>0
        cils = cil_rain;
    else

        cils = X.cil(1);
    end

end


if cils<0
    error('cils is smaller than 0')
end
out_X1.cils = cils;
D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
if any(abs(D_cils)>200)&&~X.frozen_flag
    error('wrong surface delta values')
end
end



