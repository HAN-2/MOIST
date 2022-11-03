function [num_4, num_5, num_8, num_10, num_9,dh_cria, b, dh, h_, theta_,cil_,redo,ql]=post_check(y,X2,h,theta,cil, n,dz,dt,sat_flag,poi,bc,scalingfactor)
global out_X out_X2 out_X1 out_X4
redo = 0;
h_(:,1)= y(end,1:n)';
T_ = y(end,n+1:2*n)';
cil_ = y(end,3*n+1:end);
D_cil = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;
theta_(:,1)=(h_(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
if ~bc
    theta_ = X2.theta_res+ (X2.theta_sat- X2.theta_res)./ (1+(abs(X2.alp.*h_)).^X2.N).^X2.M.*scalingfactor;
    theta_(theta_>X2.theta_sat) = X2.theta_sat(theta_>X2.theta_sat);
end
index = find(h_>=X2.he);
if ~isempty(index)
    theta_(index,1) = out_X2.theta_sat(index,1);
end


if sat_flag
    index1 = find(theta>=X2.theta_sat);
    index1 = find(h>=X2.he);
    h_(:,1)= y(end,1:n)';
    if ~isreal(h_)
        error('complex h_21')
    end
    theta_(:,1)=(h_(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
    if ~bc
        theta_ = X2.theta_res+ (X2.theta_sat- X2.theta_res)./ (1+(abs(X2.alp.*h_)).^X2.N).^X2.M.*scalingfactor;
    end
    %index = find(theta_>=X2.theta_sat);
    index = find(h_>=X2.he);
    if ~isempty(index)
        theta_(index,1) = out_X2.theta_sat(index,1);
    end
    T_ = y(end,n+1:2*n)';
    cil_ = y(end,3*n+1:end);

    if ~isempty(index1)
        theta_(index1,1) = theta(index1,1) + dt*out_X1.dydt(index1,1);
    end
    index2 = find(theta_>=X2.theta_sat);
    if ~isempty(index2)
        theta_(index2,1) = out_X2.theta_sat(index2,1);
    end

    h_(:,1) = ((theta_-out_X2.theta_res)./(X2.theta_sat-out_X2.theta_res)).^(-1./X2.lamda).*X2.he;
    if ~bc
        S = (theta_-X2.theta_res)./(X2.theta_sat-X2.theta_res);
        h_ = -abs(((S./scalingfactor).^(-1./X2.M)-1).^(1./X2.N)./X2.alp);
    end

    if ~isreal(h_)
        redo = 1;
    else
        redo = 0;
    end
end



dh_cria = zeros(n,1);
dh = abs(h_- h);
dh_cria( abs(h_)>=10,1) = 10;
dh_cria(abs(h_)<10,1) = 0.1;%0.4;
%dh_cria = 0.1;
if any(h_>= -0.1)
    dh_cria = 0.01;%0.001;
end





%-------------------------------------------------------%

%situation 4:for h oscillation (a sharp peak should not appear)
h__ = roundn(h_,-5);
num_41 = h__(1:n-2);
num_42 = h__(2:n-1);
num_43 = h__(3:n);
num_4 = find(num_42>num_43&num_42>num_41|(num_42<num_43&num_42<num_41));
if isempty(num_4) %|| length(num_4) == 1
    num_4(1:50,1)=1000;
end

nn = length(num_4);
num_5 = num_4(2:nn)-num_4(1:nn-1);
num_5 = find(num_5 == 1);

%situation 7:for cil oscillation (a sharp peak should not appear)
cil__ = roundn(cil_,-5);
num_71 = cil__(1:n-2);
num_72 = cil__(2:n-1);
num_73 = cil__(3:n);
num_7 = find(num_72>num_73&num_72>num_71|(num_72<num_73&num_72<num_71));
if isempty(num_7) %|| length(num_4) == 1
    num_7(1:5,1)=1000;
end
nn = length(num_7);
num_8 = num_7(2:nn)-num_7(1:nn-1);
num_8 = find(num_8 == 1);

%situation 6: Mass balance
cvsat_ = 610.78 .* exp(17.27 .* T_ ./ (T_ + 237.3)) .* out_X.M ./ out_X.R ./ out_X.rho ./ (T_ + 273.15);
hr_ = exp(out_X.M .* h_ .* out_X.g ./ out_X.R ./ (T_+273.15));
left = (theta_-theta)./dt+((out_X2.theta_sat-theta_).*cvsat_.*hr_-(out_X2.theta_sat-theta).*out_X1.cvsat.*out_X1.hr)./dt;
if ~bc

    left = out_X1.C.*(h_-h)./dt + ((out_X2.theta_sat-theta_).*cvsat_.*hr_-(out_X2.theta_sat-theta).*out_X1.cvsat.*out_X1.hr)./dt;
end
right = -(out_X1.q(2:end)-out_X1.q(1:end-1))./dz- out_X1.sink;
b=left-right;


%----------------------
%situation 8:for ql oscillation (a sharp peak should not appear)
S = scalingfactor.*1./ (1+(abs(X2.alp.*h_)).^X2.N).^X2.M;
k = X2.ksat.*(S./scalingfactor).^X2.l.*(1-(1-(S./scalingfactor).^(1./X2.M)).^X2.M).^2;
w = 0.5;
kk = w.*out_X1.k(2:n,1)+(1-w).*out_X1.k(1:n-1,1);
ql(2:n,1)=-(kk)./1.*(h_(2:n,1)-h_(1:n-1,1))./out_X.dz+(kk)./1*0.996;
if h_(end)>=0
    ql(n+1,1) = ql(n,1);
else
    ql(n+1,1)=0;
end
num_81 = ql(1:end-2);%roundn(ql(1:end-2),-8);
num_82 = ql(2:end-1);%roundn(ql(2:end-1),-8);
num_83 = ql(3:end);%roundn(ql(3:end),-8);
num_9 = find(num_82>num_83&num_82>num_81|(num_82<num_83&num_82<num_81));
num_9 = num_9 + 1;
nn = length(num_9);
num_10 = num_9(2:nn)-num_9(1:nn-1);
num_10 = find(num_10 == 1);


end