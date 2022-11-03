function [X5op,X5ip2_update]=climate_function(X5ip2,X5ip1,X,i)
persistent X5op_pre Month Day Year Tmax Tmin Sunshinehour Daily_hr...
    hrmax hrmin u Tp_Day Ep_Day ETp qprec_Day
climate_i_start_clock_orignal=0;
if ~X.itm
    if X5ip2.CalTT>=15*60*(X5ip2.climate_i_day-1)+(15)*60
        X5ip2.climate_i_day=X5ip2.climate_i_day+1;
        X5ip2.climate_i_day=ceil(X5ip2.CalTT/900);
        X5ip2.day_flag=true;
        dtflag=true;
    end
    if X5ip2.day_flag
        Ta=X5ip1.Daily_climate_record(X5ip2.climate_i_day,6);
        hra=X5ip1.Daily_climate_record(X5ip2.climate_i_day,8)/100;
        Sd=X5ip1.Daily_climate_record(X5ip2.climate_i_day,10);
        ET_ignored=X5ip1.Daily_climate_record(X5ip2.climate_i_day,11);
        u=X5ip1.Daily_climate_record(X5ip2.climate_i_day,12);
        zu=2;                                                                     % height of weed speed measurement, m
        canopy_height = 2.5;
        dis=canopy_height*0.67;

        zom=0.123*canopy_height;
        zoh=0.1*zom;
        ra=(log((zu-dis+zom)/zom)+0)*(log((zu-dis+zoh)/zoh+0))/0.1681/u;
        rbh=ra;
        qprec_Day=X5ip1.Rain_record(X5ip2.climate_i_day,3);

        SCF=1-exp(-0.463*X5ip2.LAI);

        Tp_Day=ET_ignored*SCF; % potential transpiration rate, m/s
        Ep_Day=ET_ignored*(1-SCF); % potential evaporation rate, m/s
        Ep=Ep_Day;
        Tp=Tp_Day;
        qprec=qprec_Day;
        if hra>0.99999
            hra=0.99;
        end
        Rnet=0.685*Sd-47.5;
        X5ip2.day_flag=false;
        %%%%%%
        if 0
            u=X5ip1.Daily_climate_record(X5ip2.climate_i_day,9);
            es_max=0.6108*exp((17.27*Tmax)/(Tmax+237.3));
            es_min=0.6108*exp((17.27*Tmin)/(Tmin+237.3));

            hrmax=min(2*Daily_hr*es_max/(es_max+es_min),0.999999999);
            hrmin=2*Daily_hr*es_min/(es_max+es_min);
            ra=(log((1+0.001)/0.001)+0)^2/0.1681/u;                          % air resistance to water vapor, s/m  from hydrus B5  %%%%1110
            ra=(log((X5ip2.zu+0.001)/0.001)+0)*(log((X5ip2.zT+0.001)/0.001)+0)/0.1681/u;

            ra=(log((X5ip2.zu-X5ip2.dis+X5ip2.zom)/X5ip2.zom)+0)*(log((X5ip2.zu-X5ip2.dis+X5ip2.zoh)/X5ip2.zoh+0))/0.1681/u;

            rbh=ra;                                                   % boundary resistance to heat transfer, s/m
            X5ip2.day_flag=false;

            [DR,ET_ignored]=Rnetfun("Daily", Month, Day, Year, X5ip2.degree, X5ip2.minute, X5ip2.South, ...
                Sunshinehour, X5ip2.elevation, Tmax, Tmin, Daily_hr, u, X5ip2.zu, X5ip2.LZ, X5ip2.LM, 0, X5ip2.root_flag, X.theta, X.T, i);

            DailyRnet(X5ip2.climate_i_day,1)=DR;
            DailyET(X5ip2.climate_i_day,1)=ET_ignored;
            qprec_Day=X5ip1.Rain_record(X5ip2.climate_i_day,3);
            SCF=1-exp(-0.463*X5ip2.LAI);

            Tp_Day=ET_ignored*SCF/1000/86400; % potential transpiration rate, m/s
            Ep_Day=ET_ignored*(1-SCF)/1000/86400; % potential evaporation rate, m/s
        else

        end
    end

else
    Ta=30;
    hra=0.2;
    Rnet=0;
end

if i==2
    X5op=struct('Ta',Ta,'hra',hra,'u',u,'Rnet',Rnet,'qprec',qprec,'qprec_Day',qprec_Day,...
        'climate_i_day', X5ip2.climate_i_day,'Ep',Ep,'Tp',Tp,'ETp',ETp,'frozen_flag', X.frozen_flag);
    X5op_pre=X5op;
else
    try
        X5op=struct('Ta',Ta,'hra',hra,'u',u,'Rnet',Rnet,'qprec',qprec,'qprec_Day',qprec_Day,...
            'climate_i_day', X5ip2.climate_i_day,'Ep',Ep,'Tp',Tp,'ETp',ETp,'frozen_flag', X.frozen_flag);
        X5op_pre=X5op;
    catch
        X5op=X5op_pre;
    end
end

X5ip2_update=X5ip2;
end