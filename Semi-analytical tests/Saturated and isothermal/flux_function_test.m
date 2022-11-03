function[q,ql,qv,qT,qi]=flux_function_test(i,X,X1,X3,X4,X2)
global out_X1 out_X4

Tcoe = 0;
qTs=0;%X1.G;
i=i/i+1;
cils = cils_function_test(i,X,X1,X3,X4);
%cils = 1.6;
cil0 = cils;
civs = cil0 * X1.cvsats * X1.hrs * X1.alphas;
civ = X.cil(:,i-1).*X1.alpha.*X1.cvsat.*X1.hr;
cvs=(X1.hrs(1)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1,1)-X1.Ts));
qis= cil0 * X1.qls - X1.Dil0 * (X.cil(1,i-1) - cil0) / (X.dz/2) - ...
    X1.Div0 * (civ(1,1) - civs) / (X.dz/2);
if X1.qls<0
    ws = 0;%0.75;%0.75;
else
    ws = 0;%0.25;
end

qis = ((1-ws)*X.cil(1) + ws.*cil0)*X1.qls -X1.Dil0(1) * (X.cil(1,i-1) - cil0) / (X.dz/2)+...
    ((1-ws)*X.cil(1) + ws.*cil0)*X1.qvs.*X4.alphadiff.*X1.alphas-...
    X1.Div0(1).*cvs.*X1.hrs.*X1.alphas.*(X.cil(1,i-1) - cil0) / (X.dz/2)-...
    X1.Div0(1).*cvs.*X1.hrs.*((1-ws)*X.cil(1) + ws.*cil0).*(X1.alpha(1,i-1) - X1.alphas) / (X.dz/2);

% frozen period
if X.frozen_flag
    qis = 0;
end
% lower boundary
if X.bottom_condition == 0
    qln=0;
    qvn=0;
    qTn=0;
    qin=0;
else
    if X.bottom_condition == 1
        qln=X1.k(X.n,1);
        qvn=0;
        qTn=0;
        qin = (qln+qvn)*X.cil(X.n,i-1);
    end
end
n=X.n;


phaie = X2.ksat.*X2.he./(1-X2.lamda.*X2.yeeta);
phai = phaie.*(X.h./X2.he).^(1-X2.lamda.*X2.yeeta);

ql(1,1)=X1.qls;
ql(2:n,1)=-(X1.k(2:n,1)+X1.k(1:n-1,1))./2.*(X.h(2:n,i-1)-X.h(1:n-1,i-1))./X.dz+(X1.k(2:n,1)+X1.k(1:n-1,1))./2;


ql(n+1,1)=0;%qln;
ql(n+1,1) = -out_X1.qevap;


qv(1,1) = X1.qvs;
qv(2:n,1)=-X3.Dv(1:n-1)./X.dz.*(X1.cvsat(2:n).*X1.hr(2:n)-X1.cvsat(1:n-1).*X1.hr(1:n-1));
qv(n+1,1) = 0;


sat_index = find(X.theta==X3.theta_sat);
if ~isempty(sat_index)
    for ii = sat_index(1):sat_index(end)
        ql(ii,1) = ql(ii+1,1);
        qv(ii,1) = qv(ii+1,1);
    end
end


q=ql+qv;

qTs = -(X1.kH(1,1)).*...
    (X.T(1,i-1)-X1.Ts)./(X.dz./2)+X1.CW.*(X.T(1,i-1)+Tcoe).*ql(1,1)+X1.CV.*qv(1,1).*(X.T(1,i-1)+Tcoe)+X1.rho.*X1.lamdaE(1,1).*qv(1,1);
qT(1,1) = qTs;
qT(2:n,1)=(-(X1.kH(1:n-1,1)+X1.kH(2:n,1))./2.*...
    (X.T(2:n,i-1)-X.T(1:n-1,i-1))./X.dz+X1.CW.*(X.T(1:n-1,i-1)+Tcoe).*ql(2:n)+X1.CV.*qv(2:n).*(X.T(1:n-1,i-1)+Tcoe)+X1.rho.*X1.lamdaE(1:n-1).*qv(2:n));
qT(2:n,1)=(-(X1.kH(2:n,1)).*...
    (X.T(2:n,i-1)-X.T(1:n-1,i-1))./X.dz+X1.CW.*(X.T(1:n-1,i-1)+Tcoe).*ql(2:n)+X1.CV.*qv(2:n).*(X.T(1:n-1,i-1)+Tcoe)+X1.rho.*X1.lamdaE(1:n-1).*qv(2:n));
qT(n+1,1) = qT(n,1);%qTn;
qT(n+1,1) = qTs;%qTn;


X1.Dil = X1.Dil; %+ %0.28*abs(ql(2:end));

qi(1,1) = qis;
qi(2:n,1)=(X.cil(1:n-1,1) + X.cil(2:n,1))./2.* ql(2:n,1) - ...
    (X1.Dil(1:n-1,1) + X1.Dil(2:n,1))./2.*(X.cil(2:n, 1) - X.cil(1:n-1,1))./ X.dz -...
    (X1.Div(1:n-1,1) + X1.Div(2:n,1))./2.*(civ(2:n, 1) - civ(1:n-1,1))./X.dz;

w(1:n) = 0;
qi(2:n,1) = ((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*ql(2:n,1) -X1.Dil(2:n) .* (X.cil(2:n,1) - X.cil(1:n-1,1)) ./ (X.dz)+...
    ((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*qv(2:n).*X4.alphadiff.*X1.alpha(2:n)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*X1.alpha(2:n).*(X.cil(2:n,1) - X.cil(1:n-1)) ./ (X.dz)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*(X1.alpha(2:n) - X1.alpha(1:n-1)) ./ (X.dz);



qi(n+1,1) = 0;
cilb = out_X4.ref*(1)*1000*0.02/0.018;
wb=0;
qi(n+1,1) = (wb*X.cil(n-1)+(1-wb)*cilb).*ql(n+1,1)- ...
    (X1.Dil(n,1)).*(cilb - X.cil(n,1))./ (X.dz/2);

end


