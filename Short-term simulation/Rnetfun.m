function [Rn,ET0,AeroTerm,RadTerm]=Rnetfun(Hourly,Month,Day,Year,degree,minute,South,...
    Sunshinehour,elevation,Tmax,Tmin,hra,u,zu,...
    LZ,LM,climate_hour_clock, root_flag, theta, T,i)
%global Rnet Ta hra
% Rnet (W/(m2)) returns net radiation
% Syntax:RNET(Hourly,Month,Year,degree,minute,South,Sunshinehour,elevation,Tmax,Tmin,hr,LZ,LM)
%             Hourly, calcualtion interval. "Hourly", "Daily","None"
%             Month,Year are the input date
%             Degree,minute are latitude
%             South or North mean Southern or Northern Hemisphere
%             Sunshinehour is the daily sunshinehour
%             Elevation
%             Tmax,Tmin are daily or hourly maximum T and minmium T
%             hr is daily or hourly relative humidity
%             u is wind speed
%             zu is height of wind speed
%             LZ is the longitude of zone where site located
%             LM is the longitude of site
% If hydrus is selected, climate data is hourly.
%%Example: [rnet,et]=Rnetfun("Daily",7,6,2013,50,48,"North",9.25,100,21.5,12.3,0.735,10,10,45,43+11/60)
% Note: "Hourly" and "None" are based on hourly climate data and assume that,
% daily is based on daily data
% temperature during one hour is constant. "Hourly" and "None" Final unit is w/m2
% "Daily" for ET calculation KJ/m2/day and mm/day

%SBc=4.903*10^(-9)*10^6;  %J/K4/m2/day


x=i;
x=x/x+1;

if strcmp(Hourly,"Hourly")||strcmp(Hourly,"None")
    SBc=4.903*10^(-9)*10^6/24;  %J/K4/m2/hour
    J=floor(275*Month/9-30+Day)-2;
    if Month<3
        J=J+2;
    end
    if mod(Year,4)==0&&Month>2
        J=J+1;
    end
else
    SBc=4.903*10^(-9)*10^6;  %J/K4/m2/day
    J=floor(275*Month/9-30+Day)-2;
    if Month<3
        J=J+2;
    end
    if mod(Year,4)==0&&Month>2
        J=J+1;
    end
end

%hourly Ra J/m2/hour

Gsc=0.082.*10^6; %J/m2/min

dr=1+0.033.*cos(2.*pi()./365.*J);%J is the day number of the year

delta=0.409.*sin(2.*pi()./365.*J-1.39);

Sc=0.1645.*sin(4.*pi().*(J-81)/364)-0.1255.*cos(2.*pi().*(J-81)./364)-0.025.*sin(2.*pi().*(J-81)./364);

if strcmp(South, "South")

    phai=-(degree+minute./60).*(pi()/180);% N(positive)/S(negative)

else

    phai=(degree+minute./60).*(pi()/180);

end

ws=acos(-tan(phai).*tan(delta));

N=24/pi().*ws; %theoritical sunshine hour

% for hourly calculation
%climate_hour_clock1=10;
t=climate_hour_clock+0.5;
if t==24.5
    t=0.5;
end

w=pi()./12.*((t+0.06667.*(LZ-LM)+Sc)-12);%t is the midpoint of clcok time (from 0.5-23.5), LZ  longtitude of centr of the loca ltime zone(120 degree for usa), LM is longitude of measurementsite

w1=w-pi().*1/24;

w2=w+pi().*1/24;

if strcmp(Hourly,"Hourly")


    Ra=12*(60)./pi().*Gsc.*dr.*((w2-w1).*sin(phai).*sin(delta)+cos(phai).*cos(delta).*(sin(w2)-sin(w1)));

    if w<-ws||w>ws
        Ra=0;
    end
end
%%%% for hourly calculation end

Gsc_daily=0.082; %MJ/m2/d
Ra_daily=24.*(60)./pi().*Gsc_daily.*dr.*(ws.*sin(phai).*sin(delta)+cos(phai).*cos(delta).*sin(ws));%MJ/m2/d

if strcmp(Hourly, "Daily")
    Gsc=0.082; %MJ/m2/d
    Ra=24.*(60)./pi().*Gsc.*dr.*(ws.*sin(phai).*sin(delta)+cos(phai).*cos(delta).*sin(ws));
end
if strcmp(Hourly,"None")
    Gsc=1366.7;
    Ra=1/pi().*Gsc.*dr.*(ws.*sin(phai).*sin(delta)+cos(phai).*cos(delta).*(sin(ws)));
end

Rs=(0.25+0.5.*Sunshinehour./N).*Ra; %rad

Rs0=(0.75+2.*10^(-5).*elevation).*Ra;

Rs_daily=(0.25+0.5.*Sunshinehour./N).*Ra_daily;
Rs0_daily=(0.75+2.*10^(-5).*elevation).*Ra_daily;

sin_e=sin(phai)*sin(delta)+cos(phai)*cos(delta)*cos(2*pi()/24*(t+0.5-12));
for i=0:23
    sum_sin_e(i+1,1)=max(sin(phai).*sin(delta)+cos(phai).*cos(delta).*cos(2.*pi()./24.*(i+1-12)),0);
end
Sum_sin_e=sum(sum_sin_e);
Rs_hydrus=max(Rs_daily*sin_e/Sum_sin_e,0);


if theta(1,x-1)<0.1
    coe=0.25;
end
if theta(1,x-1)>=0.25
    coe=0.1;
end
if theta(1,x-1)>=0.1&&theta(1,x-1)<0.25
    coe=0.35-theta(1,x-1);
end

Rns=(1-coe).*Rs; % or select coe as 0.23.
Rns_hydrus=(1-coe).*Rs_hydrus*1000000/24;%J/m2/hour

%%%Rns from hydrus   all in MJ/m2/day


Tt=(0.25+0.5.*Sunshinehour./N);
Rns_hydrus_hourly=(1-coe).*Ra_daily*max(sin_e,0)*Tt; %Rns_hydrus_daily here is based on hourly time period


%%%%%%%%%%
%es=0.6108.*exp(17.27.*(Tmax+Tmin)./2./((Tmax+Tmin)./2+237.3));

esmin=0.6108.*exp(17.27.*(Tmin)./((Tmin)+237.3));

esmax=0.6108.*exp(17.27.*(Tmax)./((Tmax)+237.3));

%ea=es.*hra; %hr~[0,1]

es=(esmin+esmax)/2;

ea=hra*100/50*esmin*esmax/(esmin+esmax);

if strcmp(Hourly,"Hourly")||strcmp(Hourly,"Daily")
    if strcmp(Hourly,"Daily")
        SBc=4.9.*10^(-9);
    end
    if Rs/Rs0<=1
        Rnl=SBc.*((Tmax+273.15).^4+(Tmin+273.15).^4)./2.*(0.34-0.14.*ea.^(0.5)).*(1.35.*Rs./Rs0-0.35);
    else
        Rnl=SBc.*((Tmax+273.15).^4+(Tmin+273.15).^4)./2.*(0.34-0.14.*ea.^(0.5)).*(1.35-0.35);
    end
    if Ra==0
        Rnl=SBc.*((Tmax+273.15).^4+(Tmin+273.15).^4)./2.*(0.34-0.14.*ea.^(0.5)).*(1.35.*0.5-0.35);
    end
else
    SBc=4.9.*10^(-9).*10^6/24./3600;
    ccc=1-Sunshinehour./N;
    eea=1.24.*(ea/((Tmax+Tmin)/2+273.15))^(1/7);
    ees=min(0.9+0.18.*theta(1,x-1),1);
    Rnl=(SBc.*((Tmax+Tmin)./2+273.15).^4.*((1-0.84.*ccc).*eea+0.84.*ccc)-SBc.*(T(1,x-1)+273.15).^4.*ees);
end

if strcmp(Hourly,"Hourly")
    Rn=Rns-Rnl;
    Rn=Rn./3600;%W/m2
    if Ra==0
        Rn_hydrus=(Rns_hydrus-Rnl)/3600;
    else
        Rn_hydrus=(Rns_hydrus+Rnl)/3600;
    end
    %Rn_hydrus=(Rns_hydrus+Rnl)/3600; %this commond could further reduce temperature
    if root_flag
        Rn=Rns_hydrus_hourly*1000000/86400;
    else
        Rn=Rn_hydrus;
    end
    %%%%%%%%%%%%%%%%%%Rn=Rn_hydrus;

else
    if strcmp(Hourly,"Daily")
        Rn=Rns-Rnl;%MJ/m2/day
        %Rn=Rn*10^6/24/3600; %J/m2/s
    end
end
if strcmp(Hourly, "None")
    Rn=Rns+Rnl;% J/m2/s
end


%-----------------About ET, daily Rn should employed, that is "Daily" option  %ET:
%if strcmp(Hourly, "Daily")
T2m=(Tmax+Tmin)./2;
es=0.611.*exp(17.27.*T2m./(T2m+237.3));
esmin=0.6108.*exp(17.27.*Tmin./(Tmin+237.3));
esmax=0.6108.*exp(17.27.*Tmax./(Tmax+237.3));
es=(esmin+esmax)/2;
%ea=es.*hra; %hr~[0,1]
ea=hra*100/50*esmin*esmax/(esmin+esmax);
u2=u.*4.87./(log(67.8.*zu-5.42));
%triangle=4098.*(0.6108.*exp(17.27.*T2m./(T2m+237.3)))./((T2m+237.3).^2);% kpa/c
triangle=2049*esmax/(Tmax+237.2)^2+2049*esmin/(Tmin+237.2)^2; % kpa/c
P=101.3.*((293-0.0065.*elevation)./293).^5.253; %Kpa
lov=(2.501-2.361*10^(-3)*T2m);%latentheat of vaporisation
%psyc=0.665.*10.^(-3).*P; % kpa/c
psyc=0.0016286*P/lov;
G=0;

canopy_height=9;
LAI=2.9;
dis=canopy_height*0.667;
%zom=0.22*(canopy_height-dis);
zom=0.123*canopy_height;
zoh=0.1*zom;
ra=log((10-dis)/zom)*log((10-dis)/zoh)/0.41/0.41/u;
rc=200/LAI;
cp=1.103/1000;%Mj/kg； specific heat
rhoa=3.486*P/(1.01*(T2m+273));%kg/m3


%-----------------bare soil------------------%
if ~root_flag
    dis=0;
    zom=0.001;
    zoh=0.001;
    ra=log((10-dis)/zom)*log((10-dis)/zoh)/0.41/0.41/u;
    rc=0;
    %ET0=1/lov*(triangle.*(Rn-G)+cp*(es-ea)/ra*rhoa*86400)/(triangle+psyc.*(1+rc/ra));%/1000/86400  mm/day
end


%-----------------------%


AerDynRes=log((10-dis)/(zom))*log((10-dis)/(zoh))/0.41^2;

raa=AerDynRes/u;
AeroTCff=0.622*3.486*86400./AerDynRes/1.01;
PAtm=101.3*((293.-0.0065*elevation)/293.)^5.253;
lambda=2.501-0.002361*T2m;
gamma=0.0016286*PAtm/lambda;
gamma1=gamma;
if raa>0
    gamma1=gamma*(1+rc/raa);
end
Dlt=2049.*esmax/(Tmax+237.3)^2+2049.*esmin/(Tmin+237.3)^2;
dl_dl=Dlt/(Dlt+gamma1);
gm_dl=gamma/(Dlt+gamma1);
EaMean=(esmax+esmin)/2;
AeroTerm=gm_dl*AeroTCff/(T2m+273.15)*u*(EaMean-ea);
RadTerm=dl_dl*Rn/lambda;
%RadTerm=1/lov*(triangle.*(Rn-G))/(triangle+psyc.*(1+rc/ra));
ET0=AeroTerm+RadTerm;

if root_flag&&strcmp(Hourly,"Hourly")
    AerDynRes=log((10-dis)/(zom))*log((10-dis)/(zoh))/0.41^2;

    raa=AerDynRes/u;
    AeroTCff=0.622*3.486*86400./AerDynRes/1.01;
    PAtm=101.3*((293.-0.0065*elevation)/293.)^5.253;
    lambda=2.501-0.002361*T2m;
    gamma=0.0016286*PAtm/lambda;
    gamma1=gamma;
    if raa>0
        gamma1=gamma*(1+rc/raa);
    end
    Dlt=2049.*esmax/(Tmax+237.3)^2+2049.*esmin/(Tmin+237.3)^2;
    dl_dl=Dlt/(Dlt+gamma1);
    gm_dl=gamma/(Dlt+gamma1);
    EaMean=(esmax+esmin)/2;
    AeroTerm=gm_dl*AeroTCff/(T2m+273.15)*u*(EaMean-ea);
    %RadTerm=dl_dl*max(0,Rns_hydrus/3600)/lambda; %Rnet here is MJ/m2/Day
    RadTerm=dl_dl*max(0,Rns_hydrus_hourly)/lambda;

    %RadTerm=1/lov*(triangle.*(Rn-G))/(triangle+psyc.*(1+rc/ra));
    ET0=AeroTerm+RadTerm;
end
