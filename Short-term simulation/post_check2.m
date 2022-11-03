function zzh = post_check2(h,n)

h__ = roundn(h,-5);
num_41 = h__(1:n-2);
num_42 = h__(2:n-1);
num_43 = h__(3:n);
num_4 = find(num_42>num_43&num_42>num_41|(num_42<num_43&num_42<num_41));
if isempty(num_4) %|| length(num_4) == 1
    num_4(1:50,1)=1000;
end

num_4 = num_4+1;
zx = [num_4(1)-1,num_4(end)];
zy= [h(num_4(1)-1),h(num_4(end))];
zxx = num_4(1)-1:1:num_4(end);
zyy = spline(zx,zy,zxx);
zzh = h;
zzh(zxx)=zyy;
end