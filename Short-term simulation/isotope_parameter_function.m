function [X4]=isotope_parameter_function(X)
Dv0=2.12e-5;
pTest = X.pTest;
itm = X.itm;
if ~X.iso_spe %O
    iso_initial=-0.008;
    alphadiff = 0.9723;%0.9691
    if pTest~=6
        alphadiff = 1;
    end

    Dca=0.9723;                                                            % Div coefficient a (1.0251 for D, 1.0285 for O)0.9755 for D
    ai=0.9669;                                                             % ai=0.9833 for D 0.9669 for 18O
    aca=1137;                                                              % alpha coefficient a (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    acb=-0.4156;                                                           % alpha coefficient b (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    acc=-0.0020667;                                                        % alpha coefficient c (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    ref=2005.2e-6;
    Dil_coe=1.026;
    if pTest==1||pTest==3
        atmos=iso_initial;
    else
        atmos=-0.015;
    end                                                                    %atmosphere isotopic concentration;
else
    iso_initial=-0.065;%-0.008
    alphadiff = 0.9755;%%0.9839;
    if pTest~=6
        alphadiff = 1;
    end

    Dil_coe=1.013;
    Dca=0.9455;                                                                % Div coefficient a (1.0251 for D, 1.0285 for O)0.9755 for D
    ai=0.9833;                                                                 % ai=0.9833 for D 0.9669 for 18O
    aca=24844;                                                                 % alpha coefficient a (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    acb=-76.284;                                                               % alpha coefficient b (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    acc=0.052612;                                                              % alpha coefficient c (aca, acb, acc, 1137,-0.4156,-0.0020667 for 18O);
    ref=155.76e-6;
    if pTest==1||pTest==3
        atmos=iso_initial;
    else
        atmos=-0.112;                                                      %atmosphere isotopic concentration;
    end
end


X4 = struct('Dv0',Dv0,'Dil_coe',Dil_coe,'alphadiff', alphadiff,'aca',aca,'acb',acb,'acc',acc,'ref',ref,...
    'pTest',pTest,'itm',itm);

end