function [num_4, num_5, num_8, dh_cria, b, dh, h_, theta_,cil_]=post_check_test(y,X2,h,theta,cil, n,dz,dt,tspan)
global out_X out_X2 out_X1 out_X4
num_8 = 0;
cil_ = 0;
index = find(theta>=X2.theta_sat);
h_(:,1)= y(end,1:n)';
theta_(:,1)=(h_(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
if ~isempty(index)
    theta_(index,1) = theta(index,1) + dt*(out_X1.dydt(index));
end
index = find(theta_>=X2.theta_sat);
theta_(index,1) = out_X2.theta_sat;
h_(:,1) = ((theta_-out_X2.theta_res)./(X2.theta_sat-out_X2.theta_res)).^(-1/X2.lamda).*X2.he;
if ~isempty(index)
    for i = index(1):index(end)
        h_(i,1) = h_(i-1,1) + dz;
    end
end

T_ = y(end,n+1:2*n)';
h__ = roundn(h_,-4);
num_41 = h__(1:n-2);
num_42 = h__(2:n-1);
num_43 = h__(3:n);
num_4 = find(num_42>num_43&num_42>num_41|(num_42<num_43&num_42<num_41));
if isempty(num_4) %|| length(num_4) == 1
    num_4(1:5,1)=1000;
end

nn = length(num_4);
num_5 = num_4(2:nn)-num_4(1:nn-1);
num_5 = find(num_5 == 1);

dh = abs(h_- h);
index_1 = abs(h_)>=10;
index_2 = abs(h_)<10;
dh_cria(index_1,1) = 10;
dh_cria(index_2,1) = 0.4;

%situation 6: Mass balance
cvsat_ = 610.78 .* exp(17.27 .* T_ ./ (T_ + 237.3)) .* out_X.M ./ out_X.R ./ out_X.rho ./ (T_ + 273.15);
hr_ = exp(out_X.M .* h_ .* out_X.g ./ out_X.R ./ (T_+273.15));
left = (theta_-theta)./dt+((out_X2.theta_sat-theta_).*cvsat_.*hr_-(out_X2.theta_sat-theta).*out_X1.cvsat.*out_X1.hr)./dt;
right = -(out_X1.q(2:end)-out_X1.q(1:end-1))./dz;
b=left-right;

end