function [cils,ws,cil_E,cil_rain,qis] = cils_function(i,X,X1,X3,X4)
global out_X1 out_X4
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
        catmos=(-200/1000+1)*X4.ref*1000*0.02/0.018;  %atmosphere isotopic composition

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
cvsata=610.78.*exp(17.27*X.Ta./(X.Ta+237.3))*X3.M./X3.R/X3.rho./(X.Ta+273.15);
hra_norm = X.hra * cvsata/cvsats;
cvs=(X1.hrs*X1.cvsat(1,1)-X3.s(1,1)*X1.hr(1)*(X.T(1,i-1)-X1.Ts));
k = cvs/cvsats/X1.hrs;
alphas= exp(-(X4.aca./(X1.Ts+273.15).^2+X4.acb./(X1.Ts+273.15)+X4.acc));%1;%
if X.pTest==1||X.pTest==2
    alphas=1;
end
rbw_re(1)=abs((cvs(1)-X3.cva)/X1.qevap);
rbw_re(1)=((cvs(1)-X3.cva)/X1.qevap);

if X1.qevap==0
    rbw_re(1)=abs((cvs(1)-X3.cva)/1e-20);
    rbw_re(1)=((cvs(1)-X3.cva)/1e-20);
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
    if cils < 0

        rbw_re = -rbw_re;
        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
    end



else

    if X.qprec>X1.qevap %&& X1.qevap > 0%theta(1)==thetasat(1)
        cils=cil_rain;
    else

        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(0)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(1-0)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
            (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
        if cils < 0

            ws = 1- ws;
            %rbw_re = -rbw_re;
            cils=(X1.Div(1)*qvs_rc/(X.dz/2)*X1.cvsat(1)*X1.hr(1)*X1.alpha(1)*X.cil(1,i-1)+civa/X1.alphak(1)/rbw_re(1)+X1.Dil(1)/(X.dz/2)*X.cil(1,i-1)-(X1.qls)*(1-ws)*X.cil(1,i-1)+X.qprec*cil_rain)/...
                (alphas/X1.alphak*cvs(1)/rbw_re(1)+(X1.qls)*(ws)+alphas*cvs(1)*X1.Div(1)*qvs_rc/(X.dz/2)+X1.Dil(1)/(X.dz/2));
        end



    end
end




if X.frozen_flag
    cils=X.cil(1,i-1);
end
if cils<0%||cils>1.8
    error('cils is smaller than 0')
end


out_X1.cils = cils;
D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
if any(abs(D_cils)>400)&&~X.frozen_flag
    %error('wrong surface delta values')
end
D_cils = (cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
D_atmos = (catmos./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
epsilon_star = (1 - alphas)*1000*X1.hrs*k;
epsilon_k = (out_X1.alphak-1)*1000*(X1.hrs*k - hra_norm);
epsilon = epsilon_star + epsilon_k;


Delta_E_cg = (X1.hrs*k*alphas*D_cils - hra_norm*D_atmos - epsilon)/( epsilon_k/1000 + (X1.hrs*k - hra_norm));
Delta_E_cg_cil1 = (X1.hrs*k*alphas*((X.cil(1)./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000) - hra_norm*D_atmos - epsilon)/( epsilon_k/1000 + (X1.hrs*k - hra_norm));
out_X1.Delta_E_cg = Delta_E_cg;
out_X1.cilp = cil_rain_qianfenzhi;
out_X1.D_atmos = D_atmos;
cil_E = 1/out_X1.alphak*(cvs(1)*cils*alphas-civa)/(cvs(1)-X3.cva);


out_X1.cil_E = cil_E;
out_X1.Delta_E_cg_cil1 = Delta_E_cg_cil1;
out_X1.acok  = k;
if X.qprec == 0
    qis = - X1.qevap * cil_E;
else
    qis = X.qprec * cil_rain - X1.qevap * cil_E;
end
qis = ((1-ws)*X.cil(1) + ws.*cils)*X1.qls -X1.Dil0(1) * (X.cil(1,i-1) - cils) / (X.dz/2)+...
    ((1-ws)*X.cil(1) + ws.*cils)*X1.qvs.*X4.alphadiff.*X1.alphas-...
    X1.Div0(1).*cvs.*X1.hrs.*X1.alphas.*(X.cil(1,i-1) - cils) / (X.dz/2)-...
    X1.Div0(1).*cvs.*X1.hrs.*((1-ws)*X.cil(1) + ws.*cils).*(X1.alpha(1,i-1) - X1.alphas) / (X.dz/2);

end



