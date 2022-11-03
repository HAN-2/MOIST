clear all
global out_X out_X1 out_X4 out_X3 out_X2 surface_flag

Initial_condition=readmatrix("stumpp","Range",'A:F',"Sheet","Initial_condition");
Rain_record=readmatrix("stumpp","Range",'A:C',"Sheet","Rain_record");
Daily_climate_record=readmatrix("stumpp","Range",'A:M',"Sheet","Daily_climate_record");
IVA_array=readmatrix("stumpp","Range",'A:C',"Sheet","Rainfall_isotope");
scalingfactor = Initial_condition (:,end);
mea = readmatrix("stumpp","Range",'B:E',"Sheet","result");

X5ip1=struct('Rain_record',Rain_record,'Daily_climate_record',Daily_climate_record);

pTest = 6;
poi = [30;90];
hcria = -1000;
bottom_condition = 1;
iso_spe = 0; % 1 is H and 0 is O
itm = 0;

%----- Time information
tMax = 1735*86400;

%------ Space information
L =1.51;
dz =0.01;
z = dz/2:dz:L;
n = length(z);
sii = 0;
bc = 0;

X=struct('iso_spe',iso_spe,'itm',itm,'pTest',pTest,'n',n,'poi',poi,'bc',bc,'scalingfactor',scalingfactor,'initial_scaling',scalingfactor);
[X2]=soil_parameter_function(X);
[X4]=isotope_parameter_function(X);
scalingfactor(1:end) = 1;

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
cil0=Initial_condition(:,4);
cim0=0;


%--------------------initial_parameters
theta(:,1)=theta0;

S = (theta(:,1)-X2.theta_res)./(X2.theta_sat-X2.theta_res);
h(:,1) = [-1.505:dz:-dz/2]';
if ~bc
    theta = X2.theta_res+(X2.theta_sat-X2.theta_res)./(1+(abs(X2.alp.*h).^X2.N)).^X2.M.*scalingfactor;
    theta(theta>=X2.theta_sat) = X2.theta_sat(theta>=X2.theta_sat);
    theta1 = theta;
end

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
degree=47;                                                                 % laitude
minute=30;                                                                 % laitude
South='North';                                                             % South or North mean Southern or Northern Hemisphere
elevation=700;                                                            % m
zu=10;                                                                     % height of weed speed measurement, m
zT=2;                                                                      % height of temperature measurement, m
LZ=345;%110;  %255                                                         % LZ is the longitude of zone where site located
LM=346+6/60;%107+88/60; %253                                                 % LM is the longitude of site
%
R1=1;
R11=1;
R111=1;

%----------------------------ROOTUPTAKE PARAMETERS------------------------%
root_flag=1;
variable_Pt=1;                                                             % 0 means potential transpiration rate is constant; 1 means potrntial trnaspiration rate is infulenced by soil water content and root distribution
canopy_height = 0;
LAI = 0;
dis = 0;
zoh = 0;
zom = 0;

dydt=zeros(4*n,1);

frozen_flag=false;
g=9.8;                                                                     % gravitational acceleration, m/s2
M=0.018;                                                                   % molar weight of water, kg/mol
R=8.314;
rho = 1000;
out_day = 1;

X=struct('itm',itm,'frozen_flag',frozen_flag,'theta',theta,'T',T,'g',g,'M',M,'R',R,'rho',rho,'variable_Pt',variable_Pt);

t_test=0;
dt_initial = 86400;
dt = dt_initial;
tic
i = 2;
flag_05=1;
total_qevap = 0;
total_rain = 0;

AbsTol = zeros(1,4*n);
AbsTol(1,1:n) = 1e-2;
AbsTol(1,n+1:2*n) = 1e-2;
AbsTol(1,2*n+1:3*n) = 1e-2;
AbsTol(1,3*n+1:4*n)= 1e-2;

options = odeset('AbsTol',AbsTol,'RelTol',1e-5,'InitialStep', dt_initial, 'MaxStep',dt_initial);
if any(theta>=X2.theta_sat)
    sat_flag = 1;
else
    sat_flag = 0;
end
while 1
    surface_flag = 1;

    CalTT = t_test;
    tspan = [t_test,t_test+dt];


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
        'g',g,'M',M,'R',R,'rho',rho,'CalTT', CalTT,'treeheight',X5ip2.treeheight,'rootlength',X5ip2.rootlength,'variable_Pt',variable_Pt,'bc',bc,'scalingfactor',scalingfactor,'initial_scaling',Initial_condition (:,end) );

    [tt,y] = ode23tb(@(t,y) coupledequations(t,y,X,X2,dydt), tspan, y0);%,options1

    if isnan(y(1))
        error('Nan results')
    end


    %%-----post_checking

    [num_4, num_5, num_8, num_10, num_9, dh_cria, b, dh, h_, theta_, cil__,redo,ql]=post_check(y, X2, h, theta, cil, n,dz,dt,sat_flag,poi,bc,scalingfactor);
    if isnan(h_(1))
        error('Nan results')
    end
    if length(num_5) >1 && any(dh>dh_cria)
        dt = dt/2;
        continue
    end


    %%-----post_checking end

    h(:,1)= y(end,1:n)';

    T(:,1)= y(end,n+1:2*n)';
    cil(:,1)=y(end,3*n+1:4*n)';
    cim(:,1)=y(end,2*n+1:3*n)';

    if ~bc
        theta = X2.theta_res+ scalingfactor.*(X2.theta_sat- X2.theta_res)./ (1+(abs(X2.alp.*h)).^X2.N).^X2.M;
        theta(theta>X2.theta_sat) = X2.theta_sat(theta>X2.theta_sat);
    end
    t_test=t_test+dt;
    i=i+1;
    y0 = [h;T;cim;cil];
    %---------------output----------------------%
    if t_test>out_day*86400 || t_test==out_day*86400
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
    z_qls(i,1) = out_X1.qls;
    z_rain(i,1) = out_X.qprec;
    z_qi(1:n+1,i) = out_X1.qi;
    z_dt(i,1) = dt;

    zz_theta(1:n,i) = theta;
    zz_h(1:n,i) = h;
    zz_sink(1:n,i) = out_X1.sink;
    zz_cil(1:n,i) = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    z_cils(1,i) = (out_X1.cils./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;

    %---------------Time Display-------------------%
    disptime=num2str(t_test/tMax*100);
    fprintf(repmat('\b',1,sii));
    sii=fprintf(disptime);
    out_X=X;
    dt = dt_initial*1;

    if any(theta>0.99*X2.theta_sat(n,1))%any(h>=X2.he)%any(theta == X2.theta_sat)
        sat_flag = 1;
        dt = 86400/2;
    else
        sat_flag = 0;
        dt = dt_initial;
    end

    options = odeset(options,'InitialStep', dt, 'MaxStep',dt);
    Delta=(cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    if CalTT>tMax
        break
    end
end
plot(z_cil(end,:))
hold on
scatter(mea(:,1),mea(:,2))
plot(mea(:,3),mea(:,4))