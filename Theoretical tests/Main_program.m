clear
global out_X out_X1 out_X4 out_X3 out_X2 surface_flag

%--------------------------Test 1 to 6 --------------------------------%
pTest =6;
poi = 0;
hcria = -20000;
bottom_condition = 0;
iso_spe = 1; % 1 is H and 0 is O
itm = 1;

%----- Time information
tMax =250*86400;

%------ Space information
L =1;
dz =0.1;
z = dz/2:dz:L;
n = length(z);
sii = 0;

X=struct('iso_spe',iso_spe,'itm',itm,'pTest',pTest,'n',n,'poi',poi);

[X4]=isotope_parameter_function_test(X);

theta = zeros(n,1);
theta_ = zeros(n,1);
h =zeros(n, 1);
T =zeros(n, 1);
cil = zeros(n, 1);
cim = zeros(n, 1);
dh_cria = zeros(n,1);

z_theta = zeros(n,250);
z_h =zeros(n,250);
z_T =zeros(n,250);
z_cil = zeros(n,250);
z_cim = zeros(n,250);


%--------------------initial condition
theta0=0.3499;
T0=30;
cil0=-0.065;
cim0=0;
% soil parameters%
theta_sat =0.35;
theta_res = 0.01;
lamda = 0.22;
yeeta = 9.14;
ksat = 1.23e-7;
he = -0.193;
X2.theta_sat =0.35;
X2.theta_res = 0.01;
X2.lamda = 0.22;
X2.yeeta = 9.14;
X2.ksat = 1.23e-7;
X2.he = -0.193;


%--------------------initial_parameters
theta(:,1)=theta0;
theta(1:end,1)=0.3499;%theta0;
% theta(1,1)= 0.43;
S = (theta-theta_res)./(theta_sat-theta_res);
h(:,1) = S.^(-1./lamda).*he;


T(:,1) = T0;
cil(:,1) = X4.ref*(cil0+1)*1000*0.02/0.018;

y0 = [h(:,1);T(:,1);cil(:,1)];

dydt=zeros(3*n,1);

frozen_flag=false;
g=9.8;                                                                     % gravitational acceleration, m/s2
M=0.018;                                                                   % molar weight of water, kg/mol
R=8.314;
rho = 1000;
out_day = 1;

X=struct('itm',itm,'frozen_flag',frozen_flag,'theta',theta,'T',T,'g',g,'M',M,'R',R,'rho',rho);

t_test=0;
dt_initial =100;%0.22;
dt = dt_initial;
tic
i = 2;
flag_05=1;
total_qevap = 0;

AbsTol = zeros(1,3*n);
AbsTol(1,1:n) = 1e-5;
AbsTol(1,n+1:2*n) = 1e-5;
AbsTol(1,2*n+1:3*n) = 1e-5;

options = odeset('AbsTol',AbsTol,'RelTol',1e-5, 'MaxStep',dt_initial, 'InitialStep', dt_initial);%'AbsTol',[1e-10 1e-10 1e-10 1e-10],

Ta = 30;
hra = 0.2;
Rnet = 0;
Ep = 2e-7;
while 1
    surface_flag = 1;

    CalTT = t_test;
    tspan = [t_test, t_test+dt];

    X = struct('h',h,'theta',theta,'T',T,'cil',cil,...
        'dydt',dydt,'i',i,'Ep',Ep,'frozen_flag',frozen_flag,'hcria',hcria,'n',n,'dz',dz,...
        'Ta',Ta,'hra',hra,'Rnet',Rnet,'root_flag',0,'he',he,'qprec',0,...
        'bottom_condition',bottom_condition,'itm',itm,...
        'poi',poi,'pTest',pTest,'iso_spe',iso_spe,...
        'g',g,'M',M,'R',R,'rho',rho,'CalTT', CalTT,'u',0);

    err = 0;
    [tt,y] = ode113(@(t,y) coupledequations_test(t,y,X,X2,dydt,err), tspan, y0,options);%

    %%-----post_checking

    [num_4, num_5, num_8, dh_cria, b, dh, h_, theta_, cil_]=post_check_test(y, X2, h, theta, cil, n,dz,dt,tspan);

    if  length(num_5) > 1 || abs(max(b))>1e-5%|| length(num_8) > 1 %(num_4(2)<5))    (length(num_4)>2&&num_4(2)<5)    any(dh>dh_cria) ||
        dt = dt/2;
        continue
    end

    %%-----post_checking end

    h(:,1) = h_;
    theta(:,1) = theta_;
    T(:,1)= y(end,n+1:2*n)';
    cil(:,1)=y(end,2*n+1:3*n)';
    theta(:,1)=(h(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
    sat_index = find( h > X2.he);
    if ~isempty(sat_index)
        theta(sat_index,1) = X2.theta_sat(1,1);
    end

    t_test=t_test+dt;
    i=i+1;
    y0 = y(end,:);

    %---------------output----------------------%
    if t_test>out_day*86400 || t_test==out_day*86400
        z_theta(:,out_day) = theta;
        z_h(:,out_day) =h;
        z_T(:,out_day) =T;
        z_cil(:,out_day) = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;

        out_day=out_day+1;
        if any(abs(z_cil)>300)
            error('wrong delta values')
        end
    end
    % total evaporation output
    total_qevap = out_X1.qevap*dt+total_qevap;
    z_qevap(i,1) = out_X1.qevap;
    z_dt(i,1) = dt;
    z_maxb(i,1) = max(b);
    z_temperature(1:n,i) = T;
    %---------------Time Display-------------------%
    disptime=num2str(t_test/tMax*100);
    fprintf(repmat('\b',1,sii));
    sii=fprintf(disptime);
    out_X=X;
    dt = dt_initial;
    if  any(theta>0.98*X2.theta_sat)&&X.qprec ~= 0%any(theta == X2.theta_sat)
        sat_flag = 1;
        dt = 50;%50;%100;%400;10
    else
        sat_flag = 0;
        dt = 500;
    end
    options = odeset(options,'InitialStep', dt,'MaxStep',dt);
    Delta=(cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
    z_delta(1:n,i) = Delta;
    if any(abs(abs(Delta)-80)>2)
        %error('bad balance')
    end
    if CalTT>tMax
        break
    end
end
toc

