function [String] = getParamsChabassier (numString, winding)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 0. Preamble
%
% getParamsChabassier() provides the pysical parameters
% of a piano string under two cases: wrapped and unwrapped
% approximations. Choose the number of the string 1-10
% sorted in ascending order of fundamental frequency.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(winding,'wrapped')
switch numString

case 1 %%%Dd1

        String.Params.f0 = 38.98;     %%% Fundamental frequency (Hz)
        String.Params.L = 1.965;      %%% Length (m)
        String.Params.d = 1.492;      %%% Diameter(mm)
        String.Params.rho = 43195;    %%% kg/m3 Material density (kg/m^3)
        String.Params.T = 1773;       %%% Longitudinal tension (N)
        String.Params.x0 = 0.236;     %%% Point of collision hammer-string

case 2 %%% C2

        String.Params.f0 = 65.56; %%% (Hz)
        String.Params.L = 1.602; %%% (m)
        String.Params.d = 1.051; %%% diameter(mm)
        String.Params.rho = 23919; %%% kg/m3
        String.Params.T = 915; %%% (N)
        String.Params.x0 = 0.192; %%% (m)


case 3 %%% Ab2

        String.Params.f0 = 104.06; %%% (Hz)
        String.Params.L = 1.591; %%% (m)
        String.Params.d = 1.095; %%% diameter(mm)
        String.Params.rho = 7589; %%% kg/m3
        String.Params.T = 783; %%% (N)
        String.Params.x0 = 0.191; %%% (m)



case 4 %%%F3

        String.Params.f0 = 175.01; %%% (Hz)
        String.Params.L = 0.96; %%% (m)
        String.Params.d = 1.049; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3 %%
        String.Params.T = 766; %%% (N)
        String.Params.x0 = 0.115; %%% (m)

case 5 %%% A3

        String.Params.f0 = 220.49; %%% (Hz)
        String.Params.L = 0.773; %%% (m)
        String.Params.d = 1.027; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 755; %%% (N)
        String.Params.x0 = 0.093; %%% (m)

case 6 %%%Cd5

        String.Params.f0 = 555.6; %%% (Hz)
        String.Params.L = 0.326; %%% (m)
        String.Params.d = 0.928; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 695; %%% (N)
        String.Params.x0 = 0.039; %%% (m)

case 7 %%%Ab5

        String.Params.f0 = 832.54; %%% (Hz)
        String.Params.L = 0.223; %%% (m)
        String.Params.d = 0.904; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 696; %%% (N)
        String.Params.x0 = 0.027; %%% (m)

case 8 %%% G6

        String.Params.f0 = 1571.35; %%% (Hz)
        String.Params.L = 0.124; %%% (m)
        String.Params.d = 0.86; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 689; %%% (N)
        String.Params.x0 = 0.015; %%% (m)

case 9 %%%C7

        String.Params.f0 = 2098.32; %%% (Hz)
        String.Params.L = 0.095; %%% (m)
        String.Params.d = 0.831; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 670; %%% (N)
        String.Params.x0 = 0.011; %%% (m)




case 10 %%%Fd7

        String.Params.f0 = 2965.15; %%% (Hz)
        String.Params.L = 0.069; %%% (m)
        String.Params.d = 0.787; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 631; %%% (N)
        String.Params.x0 = 0.008; %%% (m)

end %%%End of switch
end %%% End of IF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%UNWRAPPED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(winding,'unwrapped')
switch numString

case 1  %%%Dd1

        String.Params.f0 = 38.98; %%% (Hz)
        String.Params.L = 3.945; %%% (m)
        String.Params.d = 1.237; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 892; %%% (N)
        String.Params.x0 = 0.473; %%% (m)

case 2 %%%C2

        String.Params.f0 = 65.55; %%% (Hz)
        String.Params.L = 2.417; %%% (m)
        String.Params.d = 1.056; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 690; %%% (N)
        String.Params.x0 = 0.29; %%% (m)

case 3 %%%Ab2

        String.Params.f0 = 104.06; %%% (Hz)
        String.Params.L = 1.564; %%% (m)
        String.Params.d = 1.16; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 879; %%% (N)
        String.Params.x0 = 0.188; %%% (m)

case 4 %%%F3

        String.Params.f0 = 175.01; %%% (Hz)
        String.Params.L = 0.96; %%% (m)
        String.Params.d = 1.146; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 914; %%% (N)
        String.Params.x0 = 0.115; %%% (m)

case 5 %%%A3

        String.Params.f0 = 220.49; %%% (Hz)
        String.Params.L = 0.773; %%% (m)
        String.Params.d = 1.102; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 870; %%% (N)
        String.Params.x0 = 0.093; %%% (m)

case 6 %%%Cd5

        String.Params.f0 = 555.60; %%% (Hz)
        String.Params.L = 0.326; %%% (m)
        String.Params.d = 0.921; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 684; %%% (N)
        String.Params.x0 = 0.041; %%% (m)

case 7 %%%Ab5

        String.Params.f0 = 832.54; %%% (Hz)
        String.Params.L = 0.223; %%% (m)
        String.Params.d = 0.868; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 637; %%% (N)
        String.Params.x0 = 0.027; %%% (m)

case 8 %%%G6

        String.Params.f0 = 1571.35; %%% (Hz)
        String.Params.L = 0.124; %%% (m)
        String.Params.d = 0.794; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 587; %%% (N)
        String.Params.x0 = 0.015; %%% (m)

case 9 %%%C7

        String.Params.f0 = 2098.32; %%% (Hz)
        String.Params.L = 0.095; %%% (m)
        String.Params.d = 0.757; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 555; %%% (N)
        String.Params.x0 = 0.011; %%% (m)

case 10 %%%Fd7

        String.Params.f0 = 2965.15; %%% (Hz)
        String.Params.L = 0.069; %%% (m)
        String.Params.d = 0.707; %%% diameter(mm)
        String.Params.rho = 7850; %%% kg/m3
        String.Params.T = 510; %%% (N)
        String.Params.x0 = 0.008; %%% (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%% ENd switch

end %%% ENd of IF


end
%%% End of function
