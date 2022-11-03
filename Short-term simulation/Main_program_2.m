% Note: dispersion is considered only ql is positive, as well as qls
clear all
global out_X out_X1 out_X4 out_X3 out_X2 surface_flag

Initial_condition=readmatrix("Magali","Range",'A:E',"Sheet","Initial_condition");
Rain_record=readmatrix("Magali","Range",'A:C',"Sheet","Rain_record");
Daily_climate_record=readmatrix("Magali","Range",'A:L',"Sheet","Daily_climate_record");
IVA_array=readmatrix("Magali","Range",'A:B',"Sheet","Rainfall_isotope");
IVA_array=readmatrix("Magali","Range",'A:C',"Sheet","Rainfall_isotope");

X5ip1=struct('Rain_record',Rain_record,'Daily_climate_record',Daily_climate_record);

pTest = 6;
poi = [100;150];
hcria = -1e10;
bottom_condition = 1;
iso_spe = 1; % 1 is H and 0 is O
itm = 0;

%----- Time information
tMax =4800*900;

%------ Space information
L = 2;
dz =0.01;
z = dz/2:dz:L;
n = length(z);
sii = 0;

X=struct('iso_spe',iso_spe,'itm',itm,'pTest',pTest,'n',n,'poi',poi);
[X2]=soil_parameter_function(X);
[X4]=isotope_parameter_function(X);

theta = zeros(n,1);
theta_ = zeros(n,1);
h =zeros(n, 1);
T =zeros(n, 1);
cil = zeros(n, 1);
cim = zeros(n, 1);
dh_cria = zeros(n,1);

z_theta = zeros(n,length(Rain_record));
z_h =zeros(n,length(Rain_record));
z_T =zeros(n,length(Rain_record));
z_cil = zeros(n,length(Rain_record));
z_cim = zeros(n,length(Rain_record));


%--------------------initial condition
theta0=Initial_condition(:,1);
T0=Initial_condition(:,2);
if iso_spe %H
    cil0=Initial_condition(:,4);
else %O
    cil0=Initial_condition(:,5);
end

cim0=0;


%--------------------initial_parameters
theta(:,1)=theta0;
S = (theta(:,1)-X2.theta_res)./(X2.theta_sat-X2.theta_res);
h(:,1) = S.^(-1./X2.lamda).*X2.he;
T(:,1) = T0;
cil(:,1) = X4.ref*(cil0+1)*1000*0.02/0.018;
cim(:,1) = cim0;
y0 = [h(:,1);T(:,1);cim(:,1);cil(:,1)];

%-------climate parameters

climate_i_day=1; %0402
climate_i_hour=1; %0402
climate_ii_hour=1; %0402 relate to CalTT
climate_i_start_clock=0;%0405
climate_i_start_clock_orignal=climate_i_start_clock;
day_flag=true;
hour_flag=true;

%---------------------------GEOGRAPHIC INFORMATION------------------------%0402
degree=35;                                                                 % laitude
minute=28;                                                                 % laitude
South='North';                                                             % South or North mean Southern or Northern Hemisphere
elevation=1220;                                                            % m
zu=10;                                                                     % height of weed speed measurement, m
zT=2;                                                                      % height of temperature measurement, m
LZ=255;%110;  %255                                                         % LZ is the longitude of zone where site located
LM=251.53;%107+88/60; %253                                                 % LM is the longitude of site
%
R1=1;
R11=1;
R111=1;

%----------------------------ROOTUPTAKE PARAMETERS------------------------%
root_flag=1;
variable_Pt=1;                                                             % 0 means potential transpiration rate is constant; 1 means potrntial trnaspiration rate is infulenced by soil water content and root distribution
canopy_height=2.5;
LAI=2.5;
if root_flag==0
    LAI=0;
    canopy_height=0;
end
dis=canopy_height*0.67;
zom=0.123*canopy_height;
zoh=0.1*zom;


dydt=zeros(4*n,1);

frozen_flag=false;
g=9.8;                                                                     % gravitational acceleration, m/s2
M=0.018;                                                                   % molar weight of water, kg/mol
R=8.314;
rho = 1000;
out_day = 1;

X=struct('itm',itm,'frozen_flag',frozen_flag,'theta',theta,'T',T,'g',g,'M',M,'R',R,'rho',rho,'LAI',LAI);

t_test=0;
dt_initial =100;
dt = dt_initial;
tic
i = 2;
flag_05=1;
total_qevap = 0;
total_rain = 0;

AbsTol = zeros(1,4*n);
AbsTol(1,1:n) = 1e-5;
AbsTol(1,n+1:2*n) = 1e-5;
AbsTol(1,2*n+1:3*n) = 1e-5;
AbsTol(1,3*n+1:4*n)= 1e-5;

options = odeset('AbsTol',AbsTol,'RelTol',1e-5, 'InitialStep', dt_initial, 'MaxStep',dt_initial);%'AbsTol',[1e-10 1e-10 1e-10 1e-10],

sat_flag = 0;
while 1
    surface_flag = 1;

    CalTT = t_test;
    tspan = [t_test, t_test+dt];

    if i==2
        X5ip2=struct('CalTT',CalTT,'climate_i_day',climate_i_day,'climate_i_hour',climate_i_hour,...
            'climate_ii_hour',climate_ii_hour,'climate_i_start_clock',climate_i_start_clock,...
            'day_flag',day_flag,'hour_flag',hour_flag,...
            'degree',degree,'minute',minute,'South',South,'elevation',elevation,'zu',zu,'zT',zT,...
            'LZ',LZ,'LM',LM,...
            'canopy_height',canopy_height,'LAI',LAI,'root_flag',root_flag,'dis',dis,'zom',zom,'zoh',zoh);
    else
        X5ip2.CalTT=CalTT;

    end

    [X5op,X5ip2_update]=climate_function(X5ip2,X5ip1,X,i);
    X5ip2=X5ip2_update;

    X = struct('h',h,'theta',theta,'T',T,'cil',cil,...
        'cim',cim,'dydt',dydt,'i',i,'Ep',X5op.Ep,'Tp',X5op.Tp,'ETp',X5op.ETp,'frozen_flag',X5op.frozen_flag,'hcria',hcria,'n',n,'dz',dz,...
        'Ta',X5op.Ta,'hra',X5op.hra,'u',X5op.u,'Rnet',X5op.Rnet,'root_flag',root_flag,'he',X2.he,'qprec',X5op.qprec,'qprec_Day',X5op.qprec_Day,...
        'climate_i_day', X5op.climate_i_day,'bottom_condition',bottom_condition,'itm',itm,'IVA_array',IVA_array,...
        'poi',poi,'pTest',pTest,'iso_spe',iso_spe,...
        'g',g,'M',M,'R',R,'rho',rho,'CalTT', CalTT,'variable_Pt', variable_Pt,'LAI',LAI);


    [tt,y] = ode113(@(t,y) coupledequations(t,y,X,X2,dydt), tspan, y0,options);%

    if isnan(y(1))
        error('Nan results')
    end
    %%-----post_checking

    [num_4, num_5, num_8, dh_cria, b, dh, h_, theta_, cil__]=post_check(y, X2, h, theta, cil, n,dz,dt,sat_flag);
    if isnan(h_(1))
        error('Nan results')
    end


    if length(num_5) >1 && any(dh>dh_cria)
        dt = dt/2;
        continue
    end

    %%-----post_checking end

    h(:,1)= y(end,1:n)';
    T(:,1)  =  y(end,n+1:2*n)';
    cil(:,1)=y(end,3*n+1:4*n)';
    cim(:,1)=y(end,2*n+1:3*n)';
    theta(:,1)=(h(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;

    sat_index = find( h > X2.he);
    if ~isempty(sat_index)
        theta(sat_index,1) = X2.theta_sat(sat_index,1);
    end


    t_test=t_test+dt;
    i=i+1;

    y0 = y(end,:)';
    if X.climate_i_day<=640 %Magali
        cil(:,1)= X4.ref*(cil0+1)*1000*0.02/0.018;
        y0(4*n-n+1:4*n,1) = X4.ref*(cil0+1)*1000*0.02/0.018;
    end

    %---------------output----------------------%
    if t_test>out_day*900 || t_test==out_day*900
        z_theta(:,out_day) = theta;
        z_h(:,out_day) =h;
        z_T(:,out_day) =T;
        z_cil(:,out_day) = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
        z_cim(:,out_day) = (cim./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
        out_day=out_day+1;
        if any(abs(z_cil)>200)
            error('wrong delta values')
        end
    end
    % total evaporation output

    z_ql(1:n+1,i) = out_X1.ql;
    z_q(1:n+1,i) =  out_X1.q;
    z_rain(i,1) = out_X.qprec;
    z_evap(i,1) = out_X1.qevap;
    z_sink(1:n,i) = out_X1.sink;
    z_Tp(1,i) = out_X.Tp;
    z_Ep(1,i) = out_X.Ep;
    z_Ta(1,i) = X.Ta;
    z_hra(1,i) = X.hra;




    zz_cil(1:n,i) = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    zz_theta(1:n,i) = theta;
    z_cils(1,i) = (out_X1.cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;


    z_dt(i,1) = dt;

    %---------------Time Display-------------------%
    disptime=num2str(t_test/tMax*100);
    fprintf(repmat('\b',1,sii));
    sii=fprintf(disptime);
    out_X=X;
    dt = dt_initial;
    if  any(theta>0.98*X2.theta_sat)%any(theta == X2.theta_sat)
        sat_flag = 1;
        dt = 100;%50;%100;%400;10
    else
        sat_flag = 0;
    end

    options = odeset(options,'InitialStep', dt, 'MaxStep',dt);

    Delta=(cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;

    if CalTT>tMax
        break
    end
   
end

toc


Measurements=readmatrix("P1--FIgures.xlsx","Range",'T3:W4802',"Sheet","Magali");
Measurements_H=readmatrix("P1--FIgures.xlsx","Range",'AE32:AI43',"Sheet","Magali");
Measurements(tMax/900+1:end,:) = [];
simulations_H = [mean(z_cil(1:10,640));mean(z_cil(1:10,724));mean(z_cil(1:10,1112));mean(z_cil(1:10,1484));mean(z_cil(1:10,1860));mean(z_cil(1:10,2244));...
    mean(z_cil(1:10,2640));mean(z_cil(1:10,3016));mean(z_cil(1:10,3428));mean(z_cil(1:10,3796));mean(z_cil(1:10,4180));mean(z_cil(1:10,4564))];%z_cil(10,2640);z_cil(10,3016);z_cil(10,3428);]
simulations_H(:,2) = [mean(z_cil(1:25,640));mean(z_cil(1:25,724));mean(z_cil(1:25,1112));mean(z_cil(1:25,1484));mean(z_cil(1:25,1860));mean(z_cil(1:25,2244));...
    mean(z_cil(1:25,2640));mean(z_cil(1:25,3016));mean(z_cil(1:25,3428));mean(z_cil(1:25,3796));mean(z_cil(1:25,4180));mean(z_cil(1:25,4564))];
simulations_H(:,3) = [z_cil(50,640);z_cil(50,724);z_cil(50,1112);z_cil(50,1484);z_cil(50,1860);z_cil(50,2244);...
    z_cil(50,2640);z_cil(50,3016);z_cil(50,3428);z_cil(50,3796);z_cil(50,4180);z_cil(50,4564)];
simulations_H(:,4) = [z_cil(80,640);z_cil(80,724);z_cil(80,1112);z_cil(80,1484);z_cil(80,1860);z_cil(80,2244);...
    z_cil(80,2640);z_cil(80,3016);z_cil(80,3428);z_cil(80,3796);z_cil(80,4180);z_cil(80,4564)];
simulations_H(:,5) = [z_cil(150,640);z_cil(150,724);z_cil(150,1112);z_cil(150,1484);z_cil(150,1860);z_cil(150,2244);...
    z_cil(150,2640);z_cil(150,3016);z_cil(150,3428);z_cil(150,3796);z_cil(150,4180);z_cil(150,4564)];

iso_period = [640;724;1112;1484;1860;2244;2640;3016;3428;3796;4180;4564];


simulations = [(z_theta(25,1:tMax/900))',z_theta(75,1:tMax/900)',z_theta(125,1:tMax/900)',z_theta(175,1:tMax/900)'];
simulations(tMax/900+1:end,:) = [];
subplot(4,1,1)
plot(Measurements(:,1),'.')
hold on
plot(simulations(:,1),'r')
hold off
subplot(4,1,2)
plot(Measurements(:,2),'.')
hold on
plot(simulations(:,2),'r')
hold off
subplot(4,1,3)
plot(Measurements(:,3),'.')
hold on
plot(simulations(:,3),'r')
hold off
subplot(4,1,4)
plot(Measurements(:,4),'.')
hold on
plot(simulations(:,4),'r')
hold off

Nash1 = (sum((Measurements(:,1)-mean(Measurements(:,1))).^2) - sum((Measurements(:,1)-(simulations(:,1))).^2))/sum((Measurements(:,1)-mean(Measurements(:,1))).^2);
Nash2 = (sum((Measurements(:,2)-mean(Measurements(:,2))).^2) - sum((Measurements(:,2)-(simulations(:,2))).^2))/sum((Measurements(:,2)-mean(Measurements(:,2))).^2);
Nash3 = (sum((Measurements(:,3)-mean(Measurements(:,3))).^2) - sum((Measurements(:,3)-(simulations(:,3))).^2))/sum((Measurements(:,3)-mean(Measurements(:,3))).^2);
Nash4 = (sum((Measurements(:,4)-mean(Measurements(:,4))).^2) - sum((Measurements(:,4)-(simulations(:,4))).^2))/sum((Measurements(:,4)-mean(Measurements(:,4))).^2);

for ii = 1:4
    R2(ii,1) = ((sum((Measurements(:,ii)-mean(Measurements(:,ii))).*(simulations(:,ii)-mean(simulations(:,ii)))))/...
        ((sum((Measurements(:,ii)-mean(Measurements(:,ii))).^2).*sum((simulations(:,ii)-mean(simulations(:,ii))).^2))^0.5))^2;
end

Measurement_total = [Measurements(:,1);Measurements(:,2);Measurements(:,3);Measurements(:,4)];
simulations_total = [simulations(:,1);simulations(:,2);simulations(:,3);simulations(:,4)];
R2(ii+1,1) = ((sum((Measurement_total-mean(Measurement_total)).*(simulations_total-mean(simulations_total))))/...
    ((sum((Measurement_total-mean(Measurement_total)).^2).*sum((simulations_total-mean(simulations_total)).^2))^0.5))^2;
Nash5 = (sum((Measurement_total-mean(Measurement_total)).^2) - sum((Measurement_total-(simulations_total)).^2))/sum((Measurement_total-mean(Measurement_total)).^2);

disp('water content Nash and R2')
Nash =[Nash1;Nash2;Nash3;Nash4;Nash5];
disp([Nash,R2])

Nash1H = (sum((Measurements_H(:,1)-mean(Measurements_H(:,1))).^2) - sum((Measurements_H(:,1)-(simulations_H(:,1))).^2))/sum((Measurements_H(:,1)-mean(Measurements_H(:,1))).^2);
Nash2H = (sum((Measurements_H(:,2)-mean(Measurements_H(:,2))).^2) - sum((Measurements_H(:,2)-(simulations_H(:,2))).^2))/sum((Measurements_H(:,2)-mean(Measurements_H(:,2))).^2);
Nash3H = (sum((Measurements_H(:,3)-mean(Measurements_H(:,3))).^2) - sum((Measurements_H(:,3)-(simulations_H(:,3))).^2))/sum((Measurements_H(:,3)-mean(Measurements_H(:,3))).^2);
Nash4H = (sum((Measurements_H(:,4)-mean(Measurements_H(:,4))).^2) - sum((Measurements_H(:,4)-(simulations_H(:,4))).^2))/sum((Measurements_H(:,4)-mean(Measurements_H(:,4))).^2);
Nash5H = (sum((Measurements_H(:,5)-mean(Measurements_H(:,5))).^2) - sum((Measurements_H(:,5)-(simulations_H(:,5))).^2))/sum((Measurements_H(:,5)-mean(Measurements_H(:,5))).^2);
Measurement_totalH = [Measurements_H(:,1);Measurements_H(:,2);Measurements_H(:,3);Measurements_H(:,4);Measurements_H(:,5)];
simulations_totalH = [simulations_H(:,1);simulations_H(:,2);simulations_H(:,3);simulations_H(:,4);simulations_H(:,5)];
NashH = (sum((Measurement_totalH-mean(Measurement_totalH)).^2) - sum((Measurement_totalH-(simulations_totalH)).^2))/sum((Measurement_totalH-mean(Measurement_totalH)).^2);
Nash_H = [Nash1H;Nash2H;Nash3H;Nash4H;Nash5H;NashH];

for ii = 1:5
    R2H(ii,1) = (sum((Measurements_H(:,ii)-mean(Measurements_H(:,ii))).*(simulations_H(:,ii)-mean(simulations_H(:,ii))))/...
        ((sum((Measurements_H(:,ii)-mean(Measurements_H(:,ii))).^2).*sum((simulations_H(:,ii)-mean(simulations_H(:,ii))).^2))^0.5))^2;
end
R2H(ii+1,1) = (sum((Measurement_totalH-mean(Measurement_totalH)).*(simulations_totalH-mean(simulations_totalH)))/...
    ((sum((Measurement_totalH-mean(Measurement_totalH)).^2).*sum((simulations_totalH-mean(simulations_totalH)).^2))^0.5))^2;
disp('Isotope Nash and R2')
disp([Nash_H,R2H])
figure
for i = 1:5
    subplot(5,1,i)
    plot(simulations_H(:,i),'red')
    hold on
    plot(Measurements_H(:,i),'Blue')
end
