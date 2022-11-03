function [rbh,rbw]=AeroRes(T1,Ta,u,Cp,root_flag,r)

Treeheight=9;
Treeheight=2.5;
zm=0.001 ;                 % roughness parameter for momentum [m]
zh=zm   ;                  % roughness parameter for heat transport [m]
dl=0  ;                  % displacement level for heat transport [m]
if root_flag
    dl=2/3*Treeheight;
    zm=0.123*Treeheight;
    zh=0.0123*Treeheight;
end
rK=0.41;
TempHeight=0.00001;

WindHeight=1000;
WindHeight=200;

Wind=u;
Ca=Cp;
g=9.8;
THeight=TempHeight/100;   % conversions to m
WHeight=WindHeight/100;
TKelvA=Ta+273.15;
TKelvS=T1+273.15;

if Wind>0
    if abs(TKelvS-TKelvA) < 0.01 % The atmosphere is neutral
        r_v=1./Wind/rK/rK*(log((THeight-dl)/zh))*(log((WHeight-dl)/zm));
        rbh=r_v;
        rbw=r_v+r;
        return
    end
    psim=0;	% The atmosphere is not neutral
    psih=0;
    for i=1:6
        uu=Wind*rK/(log((WHeight-dl+zm)/zm)+psim);
        r_v=1/uu/rK*(log((THeight-dl+zh)/zh)+psih);
        if i==1
            rvMin= 0.1*r_v;
            rvMax=10*r_v;
        else
            if r_v>rvMax
                r_v=rvMax;
                break
            else
                if r_v<rvMin
                    r_v=rvMin;
                    break
                end
            end
            Mo=-Ca*TKelvA*uu*uu*uu/rK/g/(Ca*(TKelvS-TKelvA)/r_v);  % Monin-Obukhov scaling length
            zeta=(THeight-dl)/Mo;	    % unitless height
            if zeta<0 		    % unstable
                xx=(max(0,1-16*zeta))^0.25;
                psih=-2*log((1+xx*xx)/2);
                psim=-2*log((1+xx)/2)-log((1+xx*xx)/2)+2*atan(xx)-pi()/2;
            else
                if zeta>0 	  %stable
                    if zeta<1
                        psih=5*zeta;
                        psim=psih;
                    else
                        if zeta>1
                            psih=5;
                            psim=psih;
                        end
                    end
                end
            end
        end

    end
else       % no wind

    Diff0=2.12e-5;
    DiffT=Diff0*(TKelvA/273.15)^2;
    r_v=THeight/DiffT;

end

rbh=r_v;
rbw=r_v+r;

end
