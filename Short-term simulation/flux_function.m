function[q,ql,qv,qT,qi]=flux_function(i,X,X1,X3,X4,X2)
global out_X1

Tcoe = 0;
qTs=X1.G;
i=i/i+1;
n=X.n;
ql(1,1)=X1.qls;
ql(2:n,1)=-(X1.k(2:n,1)+X1.k(1:n-1,1))./2.*(X.h(2:n,i-1)-X.h(1:n-1,i-1))./X.dz+(X1.k(2:n,1)+X1.k(1:n-1,1))./2;

ql(n+1,1)=X1.k(n);

out_X1.ql = ql;

[cils,ws,cil_E,cil_rain,qis] = cils_function(i,X,X1,X3,X4);

cil0 = cils;
civs = cil0 * X1.cvsats * X1.hrs * X1.alphas;
cvs=(X1.hrs(1)*X3.cvsat(1,1)-X3.s(1,1)*X3.hr(1)*(X3.T(1,1)-X1.Ts));

civs = cil0 *cvs* X1.alphas;
civ = X.cil(:,1).*X1.alpha.*X1.cvsat.*X1.hr;

dis_coe =  0.056*0.8;
if ql(1) > 0
    X1.Dil0 = X1.Dil0 + dis_coe * abs(ql(1));
end

qis = ((1-ws)*X.cil(1) + ws.*cil0)*X1.qls -X1.Dil0(1) * (X.cil(1,i-1) - cil0) / (X.dz/2)+...
    ((1-ws)*X.cil(1) + ws.*cil0)*X1.qvs.*X4.alphadiff.*X1.alphas-...
    X1.Div0(1).*cvs.*X1.hrs.*X1.alphas.*(X.cil(1,i-1) - cil0) / (X.dz/2)-...
    X1.Div0(1).*cvs.*X1.hrs.*((1-ws)*X.cil(1) + ws.*cil0).*(X1.alpha(1,i-1) - X1.alphas) / (X.dz/2);

if abs(X.qprec * cil_rain-X1.qevap.*cil_E-qis)>1e-10
    qis = X.qprec * cil_rain-X1.qevap.*cil_E;
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

phaie = X2.ksat.*X2.he./(1-X2.lamda.*X2.yeeta);
phai = phaie.*(X.h./X2.he).^(1-X2.lamda.*X2.yeeta);



qv(1,1) = X1.qvs;
qv(2:n,1)=-X3.Dv(1:n-1)./X.dz.*(X1.cvsat(2:n).*X1.hr(2:n)-X1.cvsat(1:n-1).*X1.hr(1:n-1));
qv(n+1,1) = qv(n,1);

% lysimeter boundary
if X.h(n,1)>=X2.he(n,1) || abs(X1.C(end)) == 0 || X.theta(n,1)>0.98*X2.theta_sat(n,1) %&& ql(n) >0
    %     if  ql(n,1) > 0
    ql(n+1,1) = X2.ksat(end);%ql(n,1);
    ql(n+1,1) = -(X1.k(n,1)+X1.k(n,1))./2.*(X2.he(end,1)-X.h(n,i-1))./(X.dz/2)+(X1.k(n,1)+X1.k(n,1))./2;
    qv(n+1,1) = 0;

    if ql(n-1)<0
        ql(n+1,1) = 0;
        qv(n+1,1) = 0;
    end
else
    ql(n+1,1) = 0;
    qv(n+1,1) = 0;
end


if ql(n+1)<0
    ql(n+1) = 0;
end
out_X1.ql = ql;
%%% lysimeter boundary end



% frozen period
if X3.frozen_flag
    qis = 0;
    ql(1,1)=0;
    qv(1,1)=0;
end


q=ql+qv;


qTs = -(X1.kH(1,1)).*...
    (X.T(1,i-1)-X1.Ts)./(X.dz./2)+X1.CW.*(X.T(1,i-1)+Tcoe).*ql(1,1)+X1.CV.*qv(1,1).*(X.T(1,i-1)+Tcoe)+X1.rho.*X1.lamdaE(1,1).*qv(1,1);
qT(1,1) = qTs;
qT(2:n,1)=(-(X1.kH(1:n-1,1)+X1.kH(2:n,1))./2.*...
    (X.T(2:n,i-1)-X.T(1:n-1,i-1))./X.dz+X1.CW.*(X.T(1:n-1,i-1)+Tcoe).*ql(2:n)+X1.CV.*qv(2:n).*(X.T(1:n-1,i-1)+Tcoe)+X1.rho.*X1.lamdaE(1:n-1).*qv(2:n));
qT(n+1,1) = qT(n,1);%qTn;

for ii = 1:n
    if ql(ii+1)>0
        X1.Dil(ii) = X1.Dil(ii) + dis_coe*abs(ql(ii+1));
    end
end

w(ql(2:end)>0,1) = 0;

w(ql(2:end)<=0,1) = 1;

qi(1,1) = qis;

qi(2:n,1)=(X.cil(1:n-1,1) + X.cil(2:n,1))./2.* ql(2:n,1) - ...
    (X1.Dil(1:n-1,1) + X1.Dil(2:n,1))./2.*(X.cil(2:n, 1) - X.cil(1:n-1,1))./ X.dz -...
    (X1.Div(1:n-1,1) + X1.Div(2:n,1))./2.*(civ(2:n, 1) - civ(1:n-1,1))./X.dz;


qi(2:n,1) = ((1-w(1:n-1)).*X.cil(1:n-1) + w(1:n-1).*X.cil(2:n)).*ql(2:n,1) -X1.Dil(2:n) .* (X.cil(2:n,1) - X.cil(1:n-1,1)) ./ (X.dz)+...
    ((1-w(1:n-1)).*X.cil(1:n-1) + w(1:n-1).*X.cil(2:n)).*qv(2:n).*X4.alphadiff.*X1.alpha(2:n)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*X1.alpha(2:n).*(X.cil(2:n,1) - X.cil(1:n-1)) ./ (X.dz)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*((1-w(1:n-1)).*X.cil(1:n-1) + w(1:n-1).*X.cil(2:n)).*(X1.alpha(2:n) - X1.alpha(1:n-1)) ./ (X.dz);


qi(2:n,1) = ((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*ql(2:n,1) -X1.Dil(2:n) .* (X.cil(2:n,1) - X.cil(1:n-1,1)) ./ (X.dz)+...
    ((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*qv(2:n).*X4.alphadiff.*X1.alpha(2:n)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*X1.alpha(2:n).*(X.cil(2:n,1) - X.cil(1:n-1)) ./ (X.dz)-...
    X1.Div(2:n).*X1.cvsat(2:n).*X1.hr(2:n).*((1-1).*X.cil(1:n-1) + 1.*X.cil(2:n)).*(X1.alpha(2:n) - X1.alpha(1:n-1)) ./ (X.dz); % magali

qi(n+1,1) = qin;
qi(n+1,1) = ql(n+1)*X.cil(n);

end


