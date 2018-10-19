%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 0. Preamble. FDTD_EB_ThetaScheme_NoLoss
%
% Implementation of FDTD scheme with free parameter
% (\theta) for numerical dispersion optimisation.
% The optimal value for theta is obtained matching
% number of numerical and analytical modes
% available within SR/2.
%
% Euler-Bernoulli beam (string) without loss under
% Simply Supported (SS) Boundary Conditions excited
% by setting Initial Conditions (IC): triangular,
% windowed randomised spectrum, and raised cosine.
%
% getParamsChabassier() provides with the pysical
% parameters of a piano string under two cases:
% 'wrapped' and 'unwrapped' approximations.
% Choose the number of the string 1-10 sorted
% in ascending order of fundamental frequency.
%
%
% Gonzalo Villegas Curulla, Edinburgh October 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Custom parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% START EDITS HERE %%%%%%%%%%%%%%%%%%%%%%
% I/O
SR = 44100;                   %%% Sampling rate (Hz)
k = 1/SR;                     %%% Time step (s)
Tf = 0.5;                     %%% Duration simulation (s)
Nf = floor(SR*Tf);            %%% Duration simulation (samples)
theta = 0.84;                 %%% Manual entered theta (try theta = 0.64)

% String
[String] = getParamsChabassier(1, 'wrapped'); % 1-10 / 'wrapped','unwrapped'

% Provide initial conditions (IC)
IC_flag = 2; %%% (1) triangular, (2) randomised, (3) raised cosine (requires: wid)
ctr = String.Params.x0;       %%% In-point of excitation (m)
u0 = 0.002;                   %%% Amplitude of excitation (m)
wid = 0.04;                   %%% Width of the raised cosine (m) (for IC)

rp = 0.02;                    %%% Read-out point
win_dur = 0.02;               %%% Fade-in-out (s)


flag_visualisation = 0;       %%% True for drawnow visualisation of the beam
refresh_rate = 60;            %%% Plot every # samples

%%%%%%%%%% END EDITS %%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Derived parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_st = String.Params.rho;                  %%% Material density (kg/m^3)
r_st = String.Params.d/2000;                 %%% Radius string (m)
L_st = String.Params.L;                      %%% Length string (m)
T_st = String.Params.T;                      %%% Tension string (N)
E_st = 2.02e11;                              %%% Young's modulus (Pa)

A_st = pi*(r_st^2);                          %%% Cross-section area (m^2)
I_st = 0.25*pi*(r_st^4);                     %%% Radius of gyration or area impulse

alpha = T_st/rho_st/A_st;
beta = E_st*I_st/rho_st/A_st;

h2 = (alpha*k^2 + sqrt(alpha^2*k^4+(32*theta-16)*beta*k^2))/(4*theta-2);
hmin = sqrt(h2);            %%% Minium grid spacing for stability

N = floor(L_st/hmin);       %%% Number of segments
h = L_st/N;                 %%% Readjust h

m = (1:200);                %%% Modes indices
modal_freqs = m/2/L_st.*sqrt(alpha + beta*(m.^2)*(pi^2)/(L_st^2)); %%% (Hz)

% Optimised theta and hmin_opt
AA = (E_st*I_st*pi^2)/(rho_st*A_st*L_st^4);
BB = T_st / (rho_st*A_st*L_st^2);
CC = - SR^2;
m2 = (-BB + sqrt(BB^2 -4*AA*CC))/(2*AA);
N_max_modes = ceil(sqrt(m2));
N = N_max_modes;
h_opt = L_st / N_max_modes;
theta_opt = 0.5 + ( alpha*(k^2)*(h_opt^2) + 4*beta*(k^2) )/(2*h_opt^4);
%%% w2(end) < w2_max (!)
w2 = alpha*((1:N_max_modes)*pi/L_st).^2 + beta * ((1:N_max_modes)*pi/L_st).^4;
w_analytical = sqrt(w2);
modal_freqs = w_analytical' / (2*pi);


% Error check on THETA
if lt(theta,0.5)
  error('Theta needs to be greater than or equal to 0.5.');
end
if gt(theta,1)
  warning('Theta should be comprised between 0.5 and 1.')
end

% Prompt user with optimal theta if different
MSSG = "Your THETA value is: %.3f. Optimal theta is: %.3f.";
MSSG = compose(MSSG,theta, theta_opt);
MSSG = char(MSSG);
disp(MSSG);
dth = abs(theta-theta_opt);   %%% Difference for threshold
if dth > 1e-4
  disp('It is recommended that you change theta to optimal theta.');
  changeT = input('Do you want to change it? Y/N: ','s');
  if changeT == 'Y' || changeT == 'y'
        theta = theta_opt;
  end
end


% Read-out
rp_int  = floor(rp/h) +1;           %%% First order Interpolator
rp_frac = rp/h - rp_int +1;
lo      = floor((rp/L_st)*N);       %%% Truncation

% Excitation Initial Conditions

xax = h*[1:N-1]';       %%% Axis along beam in physical units
%%% (1) TRIANGULAR
if IC_flag == 1
      tri = u0 + u0*min(xax/ctr-1,0) + u0*min((L_st-xax)/(L_st-ctr)-1,0);
      IC = tri;
end

%%% (2) WINDOWED RANDOMNESS
if IC_flag == 2
  IC = 2*rand(N-1,1) - 1;
  [b,a] = butter(4,0.4); %%% Cutoff freq 0.4*SR --apply Butterworth
  IC = filter(b,a,IC);
  IC = IC.*abs(2*u0);    %%% Amplitude
  IC = IC.*(1-cos(2*pi*(1:N-1)'/(N-1)))*0.5; %%% Windowed IC
end

%%% (3) RAISED COSINE
if IC_flag == 3
      ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
      rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid))*u0;
      IC = rc;
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Initialisations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


u2 = IC;
u1 = u2;
u = zeros(N-1,1);
out = zeros(Nf,1);

% Matrix update

Dxxxx = toeplitz([6/(h^4) -4/(h^4) 1/(h^4) zeros(1,N-1-3)]); %%% (N-1)x(N-1)
Dxxxx(1,1) = 5/(h^4); Dxxxx(end,end) = 5/(h^4); %% (BC)
Dxxxx = sparse(Dxxxx);

Dxx = toeplitz([-2/(h^2) 1/(h^2) zeros(1,N-1-2)]);          %%% (N-1)x(N-1)
Dxx = sparse(Dxx);

v = ones(N,1);
Dxplus = 1/h*spdiags([-v, v],-1:0,N,N-1);


% System matrix-form update
A = eye(N-1) + 0.5*(1-theta)*(h^2)*Dxx;
B = alpha*(k^2)*Dxx;
C = beta*(k^2)*Dxxxx;
update = sparse(A\(B-C) + 2*eye(N-1));


% Energy coefficients
Ecoef1 = 0.5*h*rho_st*pi*(r_st^2)/k^2;
Ecoef2 = -0.25*(1-theta)*(h^3)*1/(k^2)*rho_st*pi*r_st^2;
Ecoef3 = 0.5*T_st*h;
Ecoef4 = 0.5*E_st*I_st*h;

% Energy vectors
Ek    = zeros(Nf,1);
Etot  = zeros(Nf,1);
Ev1   = zeros(Nf,1);
Ev2   = zeros(Nf,1);
Ev3   = zeros(Nf,1);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for n = 1:Nf %Nf

      u = update*u1 - u2;

      if flag_visualisation & ~mod(n,refresh_rate) %%% if PLOT
        figure(1);
        plot(xax, u, 'k');
        axis([0 L_st -1.5*u0 1.5*u0]);
        drawnow
      end


      % out(n) = u(lo); % Truncation read-out
      out(n) = (1-rp_frac)*u(rp_int) + rp_frac*u(rp_int+1); % 1st order interp

      % Energy
      Ek(n)  = Ecoef1*(u-u1)' * (u-u1);
      Ev1(n) = Ecoef2*(Dxplus*(u-u1))'*(Dxplus*(u-u1));
      Ev2(n) = Ecoef3* (Dxplus*u)' * (Dxplus*u1);
      Ev3(n) = Ecoef4* (Dxx*u)'*(Dxx*u1);


      %  Shift state
      u2 = u1; u1 = u;

end
toc

Etot = Ek+Ev1+Ev2+Ev3;
Enorm = 1 - Etot/Etot(1);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Outs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot setups
t = (0:k:Tf-k)';              %%% Time vector (s)
lineWidth   = 1.1;
font_title  = 14;
font_axis   = 12;
font_legend = 12;
font_label  = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 1. Spectra.

f1 = figure('name','MyResults','units','normalized','position',[0,0,1,1]);
f1.Color = [0.6 0.6 0.6];

OFFT = 20*log10(abs(fft(out)));
df = SR/length(OFFT);
fax = [0:df:SR-df]';

ax1 = axes(f1,'Units','Normalized');
% ax1.Position = [0.075 0.6 0.35 0.35];
ax1.OuterPosition = [0.0 0.5 0.5 0.5];

semilogx(ax1,fax,OFFT-max(OFFT),'linewidth',lineWidth);
if String.Params.f0<2500
  ax1.XLim = [0.9*String.Params.f0 2500];
else
  ax1.XLim = [String.Params.f0 7500];
end

line(ax1,[modal_freqs(:) modal_freqs(:)], get(gca, 'ylim'),'color','r');
title(ax1,'Spectra','FontSize',font_title);
xlabel(ax1,'Magnitude (dbFS)','FontSize', font_axis);
ylabel(ax1,'Frequency (Hz)','FontSize', font_axis);

xticks(ax1,[10 20 50 100 200 500 1e3 2e3 5e3 1e4]);
xticklabels(ax1,{'10','20','50','100','200','500','1k','2k','5k','10k'});
legend(ax1,{'\bf FDTD','\bf Cont. Modes'}, 'location','northwest');
ax1.Legend.FontSize = font_legend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2. Waveform

fade_N = floor(win_dur*SR);
fade_half_win = 0.5 - 0.5 * cos(pi*(fade_N-1:-1:0)/fade_N)';
out(end-fade_N+1:end) = out(end-fade_N+1:end).*fade_half_win;
out(1:fade_N) = out(1:fade_N).*flipud(fade_half_win);
out = out/max(abs(out));

ax2 = axes(f1,'Units','Normalized');
% ax2.Position = [0.56 0.6 0.38 0.35];
ax2.OuterPosition = [0.5 0.5 0.5 0.5];
title(ax2,'Waveform','FontSize',font_title);

plot(ax2,[1:Nf]'*k,out);
ax2.FontSize = font_axis;
xlim(ax2,[0 Tf]);
xlabel(ax2,'Time (s)','fontsize',font_label);
ylabel(ax2,'Displacement (m)','fontsize',font_label);
title(ax2,'Waveform','fontsize',font_title);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 3. Energy analysis.

ax3 = axes(f1,'Units','Normalized');
ax3.OuterPosition = [0 0 0.5 0.5];
title(ax3,'Ek Ev','FontSize',font_title);

plot(ax3,(1:Nf),Ek,'linewidth',lineWidth); hold on
plot(ax3,(1:Nf),Ev2+Ev3,'linewidth',lineWidth);
plot(ax3,(1:Nf),Etot,':','color','k','linewidth',lineWidth);
ax3.FontSize = font_axis;

xlim(ax3,[0 300]);
ylabel(ax3,'Energy (J)','fontsize',font_label);
legend(ax3,{'Kinetic ','Potential','Total'},'fontsize',font_legend);
title(ax3,['Theta = ',num2str(theta),' (\theta_{opt} = ',num2str(theta_opt),')'],'fontsize',font_title );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 4. Normalised Energy.
ax4 = axes(f1,'Units','Normalized');
ax4.OuterPosition = [0.5 0. 0.5 0.5];

scatter(ax4,1:150,Enorm(1:150),14.3, 'filled');
xlabel(ax4,'Samples','FontSize',font_label);
ylabel(ax4,'\delta H Enorm: 1 - H_{tot} / H_0 ','FontSize',font_label);
title(ax4,'E conservation','FontSize',font_title);
