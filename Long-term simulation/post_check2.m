function [zzh,theta_] = post_check2(h_,n,X2,scalingfactor,bc) 

    h__ = roundn(h_,-8);
    num_41 = h__(1:n-2);
    num_42 = h__(2:n-1);
    num_43 = h__(3:n);
    num_4 = find(num_42>num_43&num_42>num_41|(num_42<num_43&num_42<num_41));


    num_4 = num_4+1;

if length(num_4)>=1
    h_(num_4) = NaN;
    zzh = fillmissing(h_,'pchip');
else
    zzh = h_;
end

    theta_(:,1)=(zzh(:,1)./X2.he).^(-X2.lamda).*(X2.theta_sat-X2.theta_res)+X2.theta_res;
   if ~bc
       theta_ = X2.theta_res+ (X2.theta_sat- X2.theta_res)./ (1+(abs(X2.alp.*zzh)).^X2.N).^X2.M.*scalingfactor;
       theta_(theta_>X2.theta_sat) = X2.theta_sat(theta_>X2.theta_sat);
   end

end