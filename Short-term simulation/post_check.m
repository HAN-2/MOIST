function [num_4, num_5, num_8, dh_cria, b, dh, h_, theta_,cil_]=post_check(y,X2,h,theta,cil, n,dz,dt,sat_flag)
global out_X out_X2 out_X1 out_X4
h_(:,1)= y(end,1:n)';
T_ = y(end,n+1:2*n)';
cil_ = y(end,3*n+1:end);
D_cil = (cil./1000.*0.018./0.02-out_X4.ref)/out_X4.ref*1000;

if ~isreal(h_)
    error('h- is not real')
end
theta_(:,1)=(h_(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
index = find(h_>=X2.he);
if ~isempty(index)
    theta_(index,1) = out_X2.theta_sat(index,1);
end

if sat_flag
    index1 = find(theta>=X2.theta_sat);
    h_(:,1)= y(end,1:n)';
    if ~isreal(h_)
        error('complex h_21')
    end
    theta_(:,1)=(h_(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;

    index = find(h_>=X2.he);
    if ~isempty(index)
        theta_(index,1) = out_X2.theta_sat(index,1);
    end
    T_ = y(end,n+1:2*n)';
    cil_ = y(end,3*n+1:end);

    if ~isempty(index1)
        theta_(index1,1) = theta(index1,1) + dt*(out_X1.dydt(index1));
    end
    index2 = find(theta_>=X2.theta_sat);
    if ~isempty(index2)
        theta_(index2,1) = out_X2.theta_sat(index2,1);
    end

    h_(:,1) = ((theta_-out_X2.theta_res)./(X2.theta_sat-out_X2.theta_res)).^(-1./X2.lamda).*X2.he;
end

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

dh = abs(h_- h);
index_1 = abs(h_)>=10;
index_2 = abs(h_)<10;
dh_cria(index_1,1) = 10;
dh_cria(index_2,1) = 0.4;
%dh_cria(index)= Inf;

%situation 5:for dh oscillation (a sharp peak should not appear)


%situation 6: Mass balance
cvsat_ = 610.78 .* exp(17.27 .* T_ ./ (T_ + 237.3)) .* out_X.M ./ out_X.R ./ out_X.rho ./ (T_ + 273.15);
hr_ = exp(out_X.M .* h_ .* out_X.g ./ out_X.R ./ (T_+273.15));
left = (theta_-theta)./dt+((out_X2.theta_sat-theta_).*cvsat_.*hr_-(out_X2.theta_sat-theta).*out_X1.cvsat.*out_X1.hr)./dt;
right = -(out_X1.q(2:end)-out_X1.q(1:end-1))./dz;
b=left-right;

cil_ = y(end,3*n+1:4*n)';

right_iso = -(out_X1.qi(2:end)-out_X1.qi(1:end-1))./dz;
left_iso = (theta_.*cil_-theta.*cil)./dt+((out_X2.theta_sat-theta_).*cvsat_.*hr_.*cil_.*out_X1.alpha-(out_X2.theta_sat-theta).*out_X1.cvsat.*out_X1.hr.*out_X.cil.*out_X1.alpha)./dt;
b_iso = right_iso-left_iso;

end