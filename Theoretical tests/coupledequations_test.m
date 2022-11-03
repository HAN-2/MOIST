function dydt = coupledequations_test(t,y,X,X2,dydt,err)
global out_X out_X1 out_X2 out_X4 out_X3

y1 = y(1:X.n); % head
y2 = y(X.n+1:X.n*2); % T
y3=  y(X.n*2+1:X.n*3); % cim

X.h = y1;
X.T = y2;
%X.cim = y3;
X.cil = y3;
S=(X.h(:,1)./X2.he).^(-X2.lamda);
X.theta(:,1)=S.*(X2.theta_sat-X2.theta_res)+X2.theta_res;
sat_index = find(X.h(:,1) >= X2.he(:,1));
X.theta(sat_index, 1) = X2.theta_sat;

out_X.h = y1;
out_X.T = y2;
out_X.cil = y3;
%out_X.cil = y4;
out_X.theta = X.theta;

i=X.i/X.i+1;

[X1,X3,X4] = parameter_function_test(X,i);

out_X=X;
out_X1=X1;
out_X2 =X2;
out_X4=X4;
out_X3=X3;

[q,ql,qv,qT,qi] = flux_function_test(i,X,X1,X3,X4,X2);

sink = 0;%rootuptake(X);
out_X.sink = sink;

for ii = 1:X.n

    if X1.C(ii) == 0 || X1.C(ii) < 1e-4 %|| abs(X1.theta(ii)-X2.theta_sat)<1e-6
        X1.C(ii) = 1e-7;
    end
end

a = X1.C-X1.cvsat.*X1.hr.*X1.C+ (X2.theta_sat-X.theta(:,i-1)).*X1.cvsat.*X1.dhrdh;
b = (X2.theta_sat-X.theta(:,i-1)).*X1.cvsat.* X1.dhrdT + (X2.theta_sat-X.theta(:,i-1)).*X1.hr .* X1.dcvsatdT;

aa = ((X2.theta_sat-X.theta(:,i-1)).*X1.cvsat.*X1.dhrdh-X1.cvsat.*X1.hr.*X1.C).*X1.rho.*X1.lamdaE;
bb =  X1.Csoil + ((X2.theta_sat-X.theta(:,i-1)).*X1.cvsat.*X1.dhrdT+(X2.theta_sat-X.theta(:,i-1)).*X1.hr.*X1.dcvsatdT).*X1.rho.*X1.lamdaE;




dydt(X.n+1:2*X.n,1) = (-aa ./ a .* ((q(2:end,1) - q(1:end-1,1)) ./ X.dz-sink) + (qT(2:end,1) - qT(1:end-1,1)) ./ X.dz) ./...
    (b./a.* aa - bb);

dydt(1:X.n,1) =-1 ./ a .* ((q(2:end,1) - q(1:end-1, 1)) ./ X.dz-sink) - b./a.*dydt(X.n+1:2*X.n,1);


dydt(2*X.n+1:3*X.n,1) = (-(qi(2:end,1) - qi(1:end-1,1)) ./ X.dz - ...
    (X.cil.*X1.C-X.cil.*X1.cvsat.*X1.hr.*X1.alpha.*X1.C+(X2.theta_sat-X.theta).*X.cil.*X1.cvsat.*X1.alpha.*X1.dhrdh).*dydt(1:X.n,1)-...
    (X.cil.*X1.hr.*X1.alpha.*X1.dcvsatdT+X.cil.*X1.cvsat.*X1.alpha.*X1.dhrdT+X.cil.*X1.cvsat.*X1.hr.*X1.dalphadT).*dydt(X.n+1:2*X.n,1).*(X2.theta_sat-X.theta))./...
    (X.theta(:,1)+(X2.theta_sat-X.theta).*X1.cvsat.*X1.hr.*X1.alpha);

dydt(2*X.n+1:3*X.n,1) = (-(qi(2:end,1) - qi(1:end-1,1)) ./ X.dz - ...
    (X.cil.*X1.C+X.cil.*X1.cvsat.*X1.hr.*X1.alpha.*X1.C+(X2.theta_sat-X.theta).*X.cil.*X1.cvsat.*X1.alpha.*X1.dhrdh).*dydt(1:X.n,1)-...
    (X.cil.*X1.hr.*X1.alpha.*X1.dcvsatdT+X.cil.*X1.cvsat.*X1.alpha.*X1.dhrdT+X.cil.*X1.cvsat.*X1.hr.*X1.dalphadT).*dydt(X.n+1:2*X.n,1).*(X2.theta_sat-X.theta))./...
    (X.theta(:,1)+(X2.theta_sat-X.theta).*X1.cvsat.*X1.hr.*X1.alpha);

out_X1.q = q;
out_X1.ql = ql;
out_X1.qv = qv;
out_X1.qT = qT;
out_X1.qi = qi;
out_X1.aa = aa;
out_X1.bb = bb;
out_X1.a = a;
out_X1.b = b;
out_X1.dydt = dydt;
end


