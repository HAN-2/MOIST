function sink = rootuptake(X)
persistent distribution_F10_tem distribution_Z_tem
% h_1~h_4 are absolute value
if X.i-1 ==1
    distribution_Z_tem = 0;
    distribution_F10_tem = 0;
end
%
if X.root_flag
    dz = X.dz;
    n = X.n;
    h = X.h;
    ct = X.i-1;
    poi = X.poi;


    sink_alpha=zeros(n,1);

    Pt=X.Tp;

    p2h=5;
    p2l=9;
    r2h=5/1000/86400;
    r2l=1/1000/86400;
    h_1(1:n,1)=0.02901;%0.08; %sk0.1;
    h_2(1:n,1)=0.065;%0.15; %sk0.25;
    if Pt>r2l&&Pt<r2h
        h_3(1:n,1)=p2h+(p2l-p2h)/(r2h-r2l)*(r2h-Pt);
    end
    if Pt<=r2l
        h_3(1:n,1)=p2l;
    end
    if Pt>=r2h
        h_3(1:n,1)=p2h;
    end
    h_4(1:n,1)=160;%1.65;%100; %sk300; %1.65

    distribution_F10=0.101477806341917;%%%%%%%0.336;%0.434;%sk %0.5477;
    distribution_Z_max=0.9;






    distribution_Z = X.rootlength;

    Temarray=find(abs(h(1:poi(1)))<h_1(1:poi(1)));
    sink_alpha(Temarray,1)=0;

    Temarray=find(abs(h(1:poi(1)))>=h_1(1:poi(1))&abs(h(1:poi(1)))<h_2(1:poi(1)));
    sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_2(Temarray)-h_1(Temarray))-h_1(Temarray)./(h_2(Temarray)-h_1(Temarray));

    Temarray=find(abs(h(1:poi(1)))>=h_2(1:poi(1))&abs(h(1:poi(1)))<h_3(1:poi(1)));
    sink_alpha(Temarray,1)=1;

    Temarray=find(abs(h(1:poi(1)))>=h_3(1:poi(1))&abs(h(1:poi(1)))<h_4(1:poi(1)));
    sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_3(Temarray)-h_4(Temarray))-h_4(Temarray)./(h_3(Temarray)-h_4(Temarray));

    Temarray=find(abs(h(1:poi(1)))>=h_4(1:poi(1)));
    sink_alpha(Temarray,1)=0;

    %%(poi(1):n)
    Temarray=find(abs(h(poi(1)+1:n))<h_1(poi(1)+1:n));
    Temarray=Temarray+poi(1);
    sink_alpha(Temarray,1)=0;

    Temarray=find(abs(h(poi(1)+1:n))>=h_1(poi(1)+1:n)&abs(h(poi(1)+1:n))<h_2(poi(1)+1:n));
    Temarray=Temarray+poi(1);
    sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_2(Temarray)-h_1(Temarray))-h_1(Temarray)./(h_2(Temarray)-h_1(Temarray));


    Temarray=find(abs(h(poi(1)+1:n))>=h_2(poi(1)+1:n)&abs(h(poi(1)+1:n))<h_3(poi(1)+1:n));
    Temarray=Temarray+poi(1);
    sink_alpha(Temarray,1)=1;

    Temarray=find(abs(h(poi(1)+1:n))>=h_3(poi(1)+1:n)&abs(h(poi(1)+1:n))<h_4(poi(1)+1:n));
    Temarray=Temarray+poi(1);
    sink_alpha(Temarray,1)=abs(h(Temarray))*1./(h_3(Temarray)-h_4(Temarray))-h_4(Temarray)./(h_3(Temarray)-h_4(Temarray));

    Temarray=find(abs(h(poi(1)+1:n))>=h_4(poi(1)+1:n));
    Temarray=Temarray+poi(1);

    %%(poi(2):n)
    Temarray=find(abs(h(poi(2)+1:n))<h_1(poi(2)+1:n));
    Temarray=Temarray+poi(2);
    sink_alpha(Temarray,1)=0;

    Temarray=find(abs(h(poi(2)+1:n))>=h_1(poi(2)+1:n)&abs(h(poi(2)+1:n))<h_2(poi(2)+1:n));
    Temarray=Temarray+poi(2);
    sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_2(Temarray)-h_1(Temarray))-h_1(Temarray)./(h_2(Temarray)-h_1(Temarray));


    Temarray=find(abs(h(poi(2)+1:n))>=h_2(poi(2)+1:n)&abs(h(poi(2)+1:n))<h_3(poi(2)+1:n));
    Temarray=Temarray+poi(2);
    sink_alpha(Temarray,1)=1;

    Temarray=find(abs(h(poi(2)+1:n))>=h_3(poi(2)+1:n)&abs(h(poi(2)+1:n))<h_4(poi(2)+1:n));
    Temarray=Temarray+poi(2);
    sink_alpha(Temarray,1)=abs(h(Temarray))*1./(h_3(Temarray)-h_4(Temarray))-h_4(Temarray)./(h_3(Temarray)-h_4(Temarray));

    Temarray=find(abs(h(poi(2)+1:n))>=h_4(poi(2)+1:n));
    Temarray=Temarray+poi(2);
    sink_alpha(Temarray,1)=0;

    %------solution for distribution_F10

    if distribution_Z>0 %|| ct==1
        %distribution_Z=1.3;%dz=0.01;
        dz = X.dz;
        y1=@(x1) x1-((log((1+exp(-24.66*x1^1.59/distribution_Z.*dz.*0))./(1+exp(-24.66*x1^1.59/distribution_Z.*distribution_Z*0.1)))+0.5.*(exp(-24.66*x1^1.59/distribution_Z.*dz.*0)-exp(-24.66*x1^1.59/distribution_Z.*distribution_Z*0.1)))./...
            (log(2./(1+exp(-24.66*x1^1.59/distribution_Z.*distribution_Z)))+0.5.*(1-exp(-24.66*x1^1.59/distribution_Z.*distribution_Z))));
        if distribution_Z >0
            distribution_F10=fzero(y1,0.5);
            distribution_F10_tem=distribution_F10;
        else
            distribution_F10=0;
            distribution_F10_tem=distribution_F10;
        end
    else
        distribution_F10=distribution_F10_tem;
    end
    if distribution_Z>0
        sink_alpha(ceil(distribution_Z/dz)+1:end,1)=0;
        distribution_b=24.66*distribution_F10^1.59/distribution_Z;



        %%%root distribution function F
        i=zeros(n,1);
        i(1:ceil(distribution_Z/dz),1)=1:ceil(distribution_Z/dz);% during luck lake validation, this line is deleted, why?
        distribution_F=(log((1+exp(-distribution_b.*dz.*(i-1)))./(1+exp(-distribution_b.*dz.*(i))))+0.5.*(exp(-distribution_b.*dz.*(i-1))-exp(-distribution_b.*dz.*(i))))./...
            (log(2./(1+exp(-distribution_b.*distribution_Z)))+0.5.*(1-exp(-distribution_b.*distribution_Z)));%Li
        distribution_F(ceil(distribution_Z/dz)+1:end,1)=0;
        sink_max=distribution_F.*Pt./dz;

        distribution_F(n)=0;
        sink_max(n)=distribution_F(n).*Pt./dz;
    else
        sink_alpha(1:end,1) = 0;
        sink_max = 0;
        distribution_F(1:n,1) = 0;
    end
    %%%%%%%%

    sink=sink_alpha.*sink_max;%.*dz;

    sink_uncomp=sink_alpha.*distribution_F.*Pt./dz;

    %------------------------------variable Pt-------------------------------%
    if X.variable_Pt==1

        sum_alphaF=0;

        sum_alphaF=sum(sink_alpha.*distribution_F.^0.5);

        if sum_alphaF~=0
            %for i=1:distribution_Z*100

            sink_beta=sink_alpha.*distribution_F.^0.5./sum_alphaF;

            sink_max=sink_beta.*Pt./dz;

            %end
        else
            sink_max(:)=0;
        end

        if isnan(sink_max)
            sink_max(:)=0;
        end

        sink=sink_alpha.*sink_max;%;
        %sink(1)=0;
    end

    distribution_Z_tem = distribution_Z;
else
    sink=zeros(X.n,1);
end

end