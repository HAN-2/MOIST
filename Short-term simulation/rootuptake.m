function sink = rootuptake(X)
persistent distribution_F10_tem

if X.root_flag
    dz = X.dz;
    n = X.n;
    h = X.h;
    ct = X.i-1;
    poi = X.poi;
    if isempty(poi) || poi(1) == 0
        poi = 1;
    end

    sink_alpha=zeros(n,1);

    Pt=X.Tp;

    p2h=0.12;
    p2l=0.6;
    r2h=5/1000/86400;
    r2l=1/1000/86400;
    h_1(1:n,1)=0.03;%0.08; %sk0.1;
    h_2(1:n,1)=0.205;%0.15; %sk0.25;
    if Pt>r2l&&Pt<r2h
        h_3(1:n,1)=p2h+(p2l-p2h)/(r2h-r2l)*(r2h-Pt);
    end
    if Pt<=r2l
        h_3(1:n,1)=p2l;
    end
    if Pt>=r2h
        h_3(1:n,1)=p2h;
    end

    h_4(1:n,1)=20;


    distribution_F10=0.101477806341917;%%%%%%%0.336;%0.434;%sk %0.5477;
    distribution_Z=2;


    if isempty(poi)|| poi(1) ==1
        Temarray=find(abs(h(1:n))<h_1(1:n));
        sink_alpha(Temarray,1)=0;

        Temarray=find(abs(h(1:n))>=h_1(1:n)&abs(h(1:n))<h_2(1:n));
        sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_2(Temarray)-h_1(Temarray))-h_1(Temarray)./(h_2(Temarray)-h_1(Temarray));

        Temarray=find(abs(h(1:n))>=h_2(1:n)&abs(h(1:n))<h_3(1:n));
        sink_alpha(Temarray,1)=1;

        Temarray=find(abs(h(1:n))>=h_3(1:n)&abs(h(1:n))<h_4(1:n));
        sink_alpha(Temarray,1)=abs(h(Temarray)).*1./(h_3(Temarray)-h_4(Temarray))-h_4(Temarray)./(h_3(Temarray)-h_4(Temarray));

        Temarray=find(abs(h(1:n))>=h_4(1:n));
        sink_alpha(Temarray,1)=0;
    end
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
    sink_alpha(Temarray,1)=0;


    %------------------------------variable Pt-------------------------------%
    lamta=0.5;
    if X.variable_Pt==1
        distribution_z_cm=distribution_Z*100;
        dz_cm=dz*100;
        p1=-1.053e-9;
        p2=7.065e-7;
        p3=-0.000171;
        p4=0.01763;
        p5=-0.6921;
        p6=14.27;

        for i=0:n-1
            distribution_F(i+1,1)=(p1/6*((dz_cm*(i+1))^6-(dz_cm*(i))^6)+p2/5*((dz_cm*(i+1))^5-(dz_cm*(i))^5)+p3/4*((dz_cm*(i+1))^4-(dz_cm*(i))^4)+p4/3*((dz_cm*(i+1))^3-(dz_cm*(i))^3)+p5/2*((dz_cm*(i+1))^2-(dz_cm*(i))^2)+p6*((dz_cm*(i+1))^1-(dz_cm*(i))^1))/...
                (p1/6*distribution_z_cm^6+p2/5*distribution_z_cm^5+p3/4*distribution_z_cm^4+p4/3*distribution_z_cm^3+p5/2*distribution_z_cm^2+p6*distribution_z_cm);
        end
        if abs(sum(distribution_F)-1)>0.1
            error("check root distribution function, which sum up is not 1")
        end
        sum_alphaF=0;

        sum_alphaF=sum(sink_alpha.*distribution_F.^lamta);


        if sum_alphaF~=0
            %for i=1:distribution_Z*100

            sink_beta=sink_alpha.*distribution_F.^lamta./sum_alphaF;

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

else
    sink=zeros(X.n,1);
end
end