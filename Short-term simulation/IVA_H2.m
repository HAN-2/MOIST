function [cil_rain_qianfenzhi,cil_rain]=IVA_H2(qprec_amount,H,X)

if strcmp(H,'H')
    ref=155.76e-6;
else
    ref=2005.2e-6;
end
qprec_amount = qprec_amount/86400/1000;
if qprec_amount > 0 && qprec_amount < 888

    if strcmp(H,'H')
        cil_rain_qianfenzhi=X.IVA_array(X.climate_i_day,3);

        cil_rain=(cil_rain_qianfenzhi/1000+1)*1000*ref*0.02/0.018;

        if isnan(cil_rain_qianfenzhi)
            cil_rain_qianfenzhi = 0;
            cil_rain = 0;
        end

    else
        cil_rain_qianfenzhi=X.IVA_array(X.climate_i_day,3);
        if isnan(cil_rain_qianfenzhi)
            ind = find(X.IVA_array(:,1)==qprec_amount);
            ind2 = ind(find(min(abs(ind-X.climate_i_day))));
            cil_rain_qianfenzhi = X.IVA_array(ind2,3);
        end
        cil_rain=(cil_rain_qianfenzhi/1000+1)*1000*ref*0.02/0.018;
    end

end
if qprec_amount==0
    cil_rain_qianfenzhi=0;
    cil_rain=0;
end
end