function [X5op,X5ip2_update]=climate_function(X5ip2,X5ip1,X,i)
persistent X5op_pre Month Day Year Tmax Tmin Sunshinehour Daily_hr...
    hrmax hrmin u Tp_Day Ep_Day ETp qprec_Day LAI
global month
climate_i_start_clock_orignal=0;
if ~X.itm
    if X5ip2.CalTT>=3600*24*(X5ip2.climate_i_day-1)+(24)*3600
        X5ip2.climate_i_day=X5ip2.climate_i_day+1;
        X5ip2.climate_i_day=ceil(X5ip2.CalTT/86400);
        X5ip2.day_flag=true;
        dtflag=true;
    end
    if X5ip2.day_flag
        Month=X5ip1.Daily_climate_record(X5ip2.climate_i_day,1);
        month = Month;
        Day=X5ip1.Daily_climate_record(X5ip2.climate_i_day,2);
        Year=X5ip1.Daily_climate_record(X5ip2.climate_i_day,3);
        Tmax=X5ip1.Daily_climate_record(X5ip2.climate_i_day,4);
        Tmin=X5ip1.Daily_climate_record(X5ip2.climate_i_day,5);
        Tavg=X5ip1.Daily_climate_record(X5ip2.climate_i_day,6);
        Sunshinehour=X5ip1.Daily_climate_record(X5ip2.climate_i_day,7);
        Daily_hr=X5ip1.Daily_climate_record(X5ip2.climate_i_day,8);


        LAI = X5ip1.Daily_climate_record(X5ip2.climate_i_day,10);
        X5ip2.treeheight = X5ip1.Daily_climate_record(X5ip2.climate_i_day,11);
        pET= X5ip1.Daily_climate_record(X5ip2.climate_i_day,12);
        X5ip2.rootlength = X5ip1.Daily_climate_record(X5ip2.climate_i_day,13);

        X.frozen_flag=false;
        if Daily_hr<888
            u=X5ip1.Daily_climate_record(X5ip2.climate_i_day,9);
            es_max=0.6108*exp((17.27*Tmax)/(Tmax+237.3));
            es_min=0.6108*exp((17.27*Tmin)/(Tmin+237.3));
            %         hrmax=Daily_climate_record(climate_i_day,10);
            %         hrmin=Daily_climate_record(climate_i_day,11);
            hrmax=min(2*Daily_hr*es_max/(es_max+es_min),0.999999999);
            hrmin=2*Daily_hr*es_min/(es_max+es_min);
            ra=(log((1+0.001)/0.001)+0)^2/0.1681/u;
            ra=(log((X5ip2.zu+0.001)/0.001)+0)*(log((X5ip2.zT+0.001)/0.001)+0)/0.1681/u;

            ra=(log((X5ip2.zu-X5ip2.dis+X5ip2.zom)/X5ip2.zom)+0)*(log((X5ip2.zu-X5ip2.dis+X5ip2.zoh)/X5ip2.zoh+0))/0.1681/u;

            rbh=ra;                                                   % boundary resistance to heat transfer, s/m

            X5ip2.day_flag=false;
            qprec_Day=X5ip1.Rain_record(X5ip2.climate_i_day,3);
            SCF=1-exp(-0.4*LAI);
            ET_ignored = pET;


            Tp_Day=ET_ignored*SCF/1000/86400; % potential transpiration rate, m/s
            Ep_Day=ET_ignored*(1-SCF)/1000/86400; % potential evaporation rate, m/s
        else
            X.frozen_flag=true;
        end
    end
    if Year == 2002 && Month == 12
        X.frozen_flag = true;
    elseif Year == 2003 && any(Month==[1;2;12])%;3
        X.frozen_flag = true;
    elseif Year == 2004 && any(Month==[1;2;3;12])%;11
        X.frozen_flag = true;
    elseif Year == 2005 && any(Month==[1;2;12])%;3;11
        X.frozen_flag = true;
    elseif Year == 2006 && any(Month==[1;2;12])%;3
        X.frozen_flag = true;
    elseif Year == 2007 && Month==1
        X.frozen_flag = true;
    end


    if X5ip2.CalTT>=3600*X5ip2.climate_ii_hour%=
        X5ip2.climate_i_hour=X5ip2.climate_i_hour+1;
        X5ip2.climate_ii_hour=X5ip2.climate_ii_hour+1;
        X5ip2.climate_ii_hour=ceil(X5ip2.CalTT/3600);
        X5ip2.hour_flag=true;
    end

    if X5ip2.hour_flag%CalT<3600*climate_i_hour&&hour_flag %assuming climate is constant within one hour
        climate_hour_clock=mod(X5ip2.climate_ii_hour+climate_i_start_clock_orignal,24);
        if climate_hour_clock==0
            climate_hour_clock=24;
        end
        Ta=(Tmax+Tmin)/2+(Tmax-Tmin)/2*cos(pi()/12*(climate_hour_clock-13));%Tavg+(Tmax+Tmin)/2*cos(2*pi()*((climate_hour_clock-15)/24));
        hra=(hrmax+hrmin)/2+(hrmax-hrmin)/2*cos(pi()/12*(climate_hour_clock-1));%Daily_hr+(hrmax+hrmin)/2*cos(2*pi()*((climate_hour_clock-6)/24));
        Ta0=(Tmax+Tmin)/2+(Tmax-Tmin)/2*cos(pi()/12*(climate_hour_clock-13));%Tavg+(Tmax+Tmin)/2*cos(2*pi()*((climate_hour_clock-15)/24));
        hr0=(hrmax+hrmin)/2+(hrmax-hrmin)/2*cos(pi()/12*(climate_hour_clock-1));
        if Daily_hr ==888
            hra = 0.4;
        end
        %hra=0.56;
        [Rnet,ETp,AeroTerm,RadTerm]=Rnetfun("Hourly", Month, Day, Year, X5ip2.degree, X5ip2.minute, X5ip2.South,...
            Sunshinehour, X5ip2.elevation, Ta, Ta, hra, u, X5ip2.zu,...
            X5ip2.LZ, X5ip2.LM, climate_hour_clock, X5ip2.root_flag, X.theta, X.T, i);

        cva=610.78*exp(17.27*Ta/(Ta+237.3))*X.M/X.R/X.rho/(Ta+273.15)*hra;     % concentration of water vapor at atmosphere
        es=610.78*exp(17.27*Ta./(Ta+237.3));                             % pa
        rhov=0.622*es*hra/287.04/(Ta+273.15);                            % T in k   kg/m3
        % Cp calculation heat capacity of air
        vapor_pressure=0.61078*exp(17.27*Ta/(Ta+237.3))*1000*hra;
        air_pressure=101.325*exp(-X5ip2.elevation/8200)*1000;
        air_rho=air_pressure/287.04/(Ta+273.15)*(1-0.378*vapor_pressure/air_pressure);
        Cp=1200;%1004.7*(0.522*vapor_pressure/air_pressure+1)*air_rho;
        if ~X5ip2.root_flag
            if climate_hour_clock<7||climate_hour_clock>17 %Hydrus Page 17
                Tp=Tp_Day*0.24;
                Ep=Ep_Day*0.24;
            else
                Tp=Tp_Day*2.75*sin(2*pi()*climate_hour_clock/24-pi()/2);% potential transpiration rate, m/s
                Ep=Ep_Day*2.75*sin(2*pi()*climate_hour_clock/24-pi()/2);% potential evaporation rate, m/s
            end
        else
            SCF=1-exp(-0.4*LAI);
            Ep=ETp*(1-SCF)/1000/86400;
            Tp=ETp*(SCF)/1000/86400;
        end

        qprec=qprec_Day*(1+cos(2*pi()*climate_hour_clock/24-pi())); % hourly prec based on dayly.
        %qprec=qprec_Day;

        if Ta~=(Tmax+Tmin)/2+(Tmax-Tmin)/2*cos(pi()/12*(climate_hour_clock-13))
            error("different Ta values");
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X5ip2.hour_flag=false;
    end%0405

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