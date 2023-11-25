clear;
clc;
close all;

% copyright: Xiaolong Zhang (UIUC), xzhng157@illinois.edu

% Test V6.1 voltage
% Test V6.2 current with LC and Resistive load
% Test V6.3 include filter inductor resistance into model
% Test V6.4 motor current prediction
% Test V6.5 motor current prediction
%           read in ia/ib/ic/vab in seperate Excel files
% Test V6.6 add motor impedance nonlinear model (2022.4)
% Test V6.8 use new filter parameters and new analytical impedance model (2022.11.15)
%           add R_c to reference voltage (modulation index) setting function
%           change capacitor bank to Y connection
%           add R_c to the calculation of Veq
% Test V6.9 add capacitor's ESL (including cable inductance)
% Test V6.10 Sawtooth carrier - use phase angle representation to calculate harmonics
%            add deadtime PWM spectrum #9,#10,#11,#12,#13 and verified
%            dead time loss is included in M calculation
% Test V6.11 Use only one phase in each sequence circuit calculation and then inverse
%            transform back to all phases
%            Phase shift in Vs1 due to deadtime is included in M calculation by vetorial sum
% Test V6.12 Modify modulation type 14 with switch and diode voltage drop
%                               distortion, changed deadtime expression, verified with test
% Test V6.13 Improved linear circuit impedance and harmonic calculation
% Test V6.14 Compare frequency-invariant method (ita0 = vsa./Zeq, line 1670-1689) vs.
%            proposed method (it1 = vs1./Zeq1; it2 = vs2./Zeq2; ita = it1 +
%            it2), constrasted plot in Figure 45, prediction error in
%            Figure 46.
% Test V6.15 Simplified harmonic synthesis using complex-number representation: Cmn_a = Amn_ph.*exp(j*Theta_a)... (line 1296-1233)

tStart = tic;

%% machine and system parameters
modulation_type = 14;
% #1, naturally sampled reference, single-edge carrier PWM
% #2, naturally sampled reference, double-edge carrier PWM
% #3, regularly sampled reference, single-edge carrier PWM
% #4, regularly sampled reference, symmetrical double-edge carrier PWM
% #5, regularly sampled reference, asymmetrical double-edge carrier PWM
% #6, naturally sampled reference, double-edge carrier SVM
% #7, regularly sampled reference, symmetrical double-edge carrier SVM
% #8, regularly sampled reference, asymmetrical double-edge carrier SVM
% #9, naturally sampled reference, double-edge carrier PWM, w/ deadtime
% #10, regularly sampled reference, symmetrical double-edge carrier PWM, w/deadtime
% #11, regularly sampled reference, aymmetrical double-edge carrier PWM, w/deadtime
% #12, regularly sampled reference, symmetrical double-edge carrier SVM, w/deadtime
% #13, regularly sampled reference, aymmetrical double-edge carrier SVM, w/deadtime

machine_model = 2;
% 1, linear model; 2, nonlinear model

%constants
C120 = exp(j*pi*2/3);
C240 = exp(-1i*pi*2/3);
s = tf ('s') ;

%drive parameters
n_pu = 0.4
fn = 1000;
fo = fn*n_pu;
Fsw_idiq_ctrl = 8e3;
p = round(Fsw_idiq_ctrl/fo);
fc = p*fo; %8*fo;
% M =  0.605*2/sqrt(3) %  for n =1;
p = fc/fo % carrier to fundamental ratio
Ht = 20*p; % 10*p % number of harmonics considered, t = total
Mt = ceil(Ht/p) + 10; %2*ceil(Ht/p)+1;
Nt = 10*p; %50*p %2*Ht; 
Kt_0 =  80; %100;% PWM spectrum
Kt_1 = 50; % deadtime spectrum
Vdc = 60*(1-0.04) %
Vnph = Vdc/2/sqrt(2); %phase to neutral, rms
Vnll = sqrt(3)*Vdc/2/sqrt(2); % line-line, rms

order = zeros(Ht,1);
for h = 1:1:Ht
    order(h) = h;
end

Frequency = order.*fo;
w = 2*pi.*Frequency; % frequency vector

%machine parameters
Phase_Resistance= 2/3*0.185; %2/3*0.183 (previous measurement)  % thingap nominal: 125.2e-3; 
Synchronous_Inductance= 13.05e-6;     
PM_Flux_Linkage= 0.007879*0.97%0.98; %*0.965; rated: V/Hz = 0.0476, phi = 0.00758, Vll/kRPM = 22.0; normal temperature: V/Hz = 0.04948, phi = 0.007879
%Rotor_Inertia=0.6627*0.2*1.0e5;
%Viscous_Damping=8.09e-3*0.2;
Pole_Pairs=16;

% rated and per unit values
In=1; % rated current, Amp-peak 
Tn=3/2*Pole_Pairs*PM_Flux_Linkage*In; % rated torque, Nm
Nn=3750; % rated speed, rpm
fn=Nn/60*Pole_Pairs; % rated electrical frequency, Hz
Vn=2*pi*fn*PM_Flux_Linkage; % rated phase voltage, Volt-peak
Zn=Vn/In; % impedance base, Ohm
Vdcn = Vn/0.577;% rated DC bus voltage corresponding to rated phase voltage
                % SVPWM: Vmax,alpha = sqrt(3)/3*Vdc = 0.577*Vdc
                % SPWM: Vmax,alpha = 1/2*Vdc = 0.50*Vdc

% filter parameters
R_filter = 1000e3;        % wye connection
Lf_0 = (305.36e-6 + 3e-6);       % series in line, cable inductance 5 uH
L_filter2 =1e-6;
kC_series = 1;
Cf_0 = 5.9e-5/kC_series;       % wye connection
L_esl = 1e-6;
R_Lf_0 = 0.125 +0.1; % cable resistance 0.2 Ohm
R_Cf_0 = 0.004;
R_d_0 = 0.21 + 0.05;%0.02;
R_c = R_Cf_0*kC_series + R_d_0 ; % 


Psi1 = PM_Flux_Linkage;   
Rs_0 = Phase_Resistance;
Ls_0 = Synchronous_Inductance;


%% frequency-dependent impedance functions


Rf = R_filter;
Xs_0 = (Synchronous_Inductance + L_filter2).*w;
Rr_0 = 0;
       
% inductor function

ro_cu = 1.8*10^-8;
mu0 = 0.00000125663706212;
airgap_factor_Lf = 1;

leakage_factor_Lf_fit = 0.0539;
d_Lf_fit = 1.0e-3;
zt_Lf_fit = 2.3461;

Lfm_0 = (1-leakage_factor_Lf_fit)*Lf_0;
Lfl_0 = leakage_factor_Lf_fit*Lf_0;

[Lf,R_Lf] = Ls_Rs_analytical(R_Lf_0,Lf_0,leakage_factor_Lf_fit,w,d_Lf_fit,ro_cu,mu0,airgap_factor_Lf,zt_Lf_fit);


%machine function

ro_PM = 0.0000018;
airgap_factor_Ls = 1;
airgap_factor_Lrl = 1;

leakage_factor_Ls_fit = 0.2569;
d_Ls_fit = 1.0008e-04;
zt_Ls_fit = 3.1278;
Rr_0_fit = 0.1244;
Lrl_0_fit = 1.2056e-04;
d_Lrl_fit = 0.0687;

Lsm_0 = (1-leakage_factor_Ls_fit)*Ls_0;
Lsl_0 = leakage_factor_Ls_fit*Ls_0;

[Ls,Rs] = Ls_Rs_analytical(Rs_0,Ls_0,leakage_factor_Ls_fit,w,d_Ls_fit,ro_cu,mu0,airgap_factor_Ls,zt_Ls_fit);
Lsl = Ls - Lsm_0;
[Lrl,Rr] = Lrl_Rr_analytical(Rr_0_fit,Lrl_0_fit,w,d_Lrl_fit,ro_PM,mu0,airgap_factor_Lrl);

Ls_complex = Lsl + Lsm_0.*(Rr + j.*w.*Lrl)./(Rr + j.*w.*(Lsm_0 + Lrl));

Leq_machine_fit = real(Ls_complex) + L_filter2;
Req_machine_fit = (Rs-w.*imag(Ls_complex));
Xeq_machine_fit = w.*Leq_machine_fit;

% positive and negative sequence impedance
w1 = abs(2*pi*fo.*(order - 1));
w2 = abs(2*pi*fo.*(order + 1));
w1(1) = 2*pi*fo*1e-3;
slip1 = w1./w;
slip2 = w2./w;

[Lrl_1,Rr_1] = Lrl_Rr_analytical(Rr_0_fit,Lrl_0_fit,w1,d_Lrl_fit,ro_PM,mu0,airgap_factor_Lrl);
Ls_complex_1 = Lsl + Lsm_0.*(Rr_1 + j.*w1.*Lrl_1)./(Rr_1 + j.*w1.*(Lsm_0 + Lrl_1));
Leq_machine_fit_1 = real(Ls_complex_1) + L_filter2;
Req_machine_fit_1 = (Rs-w.*imag(Ls_complex_1));
Xeq_machine_fit_1 = w.*Leq_machine_fit_1;

[Lrl_2,Rr_2] = Lrl_Rr_analytical(Rr_0_fit,Lrl_0_fit,w2,d_Lrl_fit,ro_PM,mu0,airgap_factor_Lrl);
Ls_complex_2 = Lsl + Lsm_0.*(Rr_2 + j.*w2.*Lrl_2)./(Rr_2 + j.*w2.*(Lsm_0 + Lrl_2));
Leq_machine_fit_2 = real(Ls_complex_2) + L_filter2;
Req_machine_fit_2 = (Rs-w.*imag(Ls_complex_2));
Xeq_machine_fit_2 = w.*Leq_machine_fit_2;
%% Per unit impedance and modulation index calculation
Xs_pu=2*pi*fn*(Synchronous_Inductance+L_filter2)/Zn;
Rs_pu=Phase_Resistance/Zn;
Vdc_pu = Vdc/Vdcn;

Xr_pu = R_filter/Zn;
Xl_pu = 2*pi*fn*Lf_0/Zn;
Xc_pu = (1/(2*pi*fn*Cf_0))/Zn; % line to neutral
Xlr_pu = R_Lf_0/Zn;
Xcr_pu = R_c/Zn;
Xcl_pu = 2*pi*fn*L_esl/Zn;

fLC = 1/sqrt(Lf_0*Cf_0)/2/pi;

wm0=n_pu*2*pi*Nn/60;
ia0=0; %200;
ib0=0; %-100;

%{
filename = ['IntegraVision.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  'IntegraVision';
opts.SelectedVariableNames = [3:3]; 
opts.DataRange = '20:59';
Test_Harmonic = readmatrix(filename,opts);

opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '18:18';
Test_THD = readmatrix(filename,opts)
%}

% reference machine current 
Id0 = 3.0 + Xc_pu/(Xc_pu^2/n_pu^2 + Xcr_pu^2);
Iq0 = 0.91 - n_pu*Xcr_pu/(Xc_pu^2/n_pu^2 + Xcr_pu^2);%Test_Harmonic(1)*sqrt(2);
% Id0 =  2.6 + n_pu^2/Xc_pu;
% Iq0 = 0.91;%Test_Harmonic(1)*sqrt(2);
Eq0 = n_pu;
It1 = Id0 + j*Iq0
Vt1 = j*Eq0 + It1*(j*Xs_pu*n_pu+Rs_pu);
Ic = Vt1/(-j*Xc_pu/n_pu + j*Xcl_pu*n_pu + Xcr_pu);
Vs1 = Vt1 + (j*Xl_pu*n_pu + Xlr_pu)*(It1+Ic);
Ud0 = real(Vs1);
Uq0 = imag(Vs1);
Is1 = Ic + It1
PF_angle = angle(Vs1)-angle(Is1);
PF = cos(PF_angle)
It1_rms = abs(It1)*In/sqrt(2)
%{
% reference inverter current 
It1_rms = 0;
Id0 = -n_pu^2/Xc_pu;
Iq0 = 0.2*sqrt(2);
Is1 = Id0 + j*Iq0;
Eq0 = n_pu;

while abs(It1_rms - Test_Harmonic(1))>0.01
    if It1_rms > Test_Harmonic(1)
        Iq0 = Iq0 - 0.01;
    else
        Iq0 = Iq0 + 0.01;
    end

    Is1 = Id0 + j*Iq0
    Vm = j*Eq0*(-j*Xc_pu/n_pu)/(j*Xs_pu*n_pu+Rs_pu-j*Xc_pu/n_pu);
    Zm = (-j*Xc_pu/n_pu)*(Rs_pu+j*Xs_pu*n_pu)/(Rs_pu+j*Xs_pu*n_pu-j*Xc_pu/n_pu);
    Vs1 = Vm + (Zm+j*Xl_pu*n_pu+Xlr_pu)*Is1;
    It1 = (Vs1 - Is1*(j*Xl_pu*n_pu+Xlr_pu) - j*Eq0)/(j*Xs_pu*n_pu+Rs_pu)
    Ic = (Vs1 - Is1*(j*Xl_pu*n_pu+Xlr_pu))/(-j*Xc_pu/n_pu + Xcr_pu)
    It1_rms = abs(It1)/sqrt(2);
end
%}
Ud0 = real(Vs1);
Uq0 = imag(Vs1);
PF_angle = angle(Vs1)-angle(Is1)
PF = cos(PF_angle)
Vs1_ll_rms = abs(Vs1)*Vn*sqrt(3)/sqrt(2);
Vt1_ll_rms = abs(Vt1)*Vn*sqrt(3)/sqrt(2);
Eq1_ll_rms = abs(Eq0)*Vn*sqrt(3)/sqrt(2);
% dead time
Td = 0;%2e-6;
TdTs = Td*fc
delta = 1;
% voltage loss due to DT
    m=0;
    n=1;
    q = m + n/p;
    Omega_mn =2*pi*q;
    Kt = 600;   
    M = abs(Vs1)*Vn*2/Vdc;
            sum3=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Gmnp = (1-(-1)^(k-n))*(exp(-1/2*j*Omega_mn)+(-1)^k);
                            sum3 = sum3 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Gmnp;
                        end
                    end                    
% Fsw_idiq_ctrl = fc;
Vdt_pu = 2*TdTs;
Vdiode = 0.7;
Vswitch = 0.8;
Vdrop_pu = (Vdiode+Vswitch)/2*2/Vdc + 2*2e-6*fc;

k_comp = cos(2*pi/p*TdTs)%*sin(wm0*Pole_Pairs/(2*Fsw_idiq_ctrl))/(wm0*Pole_Pairs/(2*Fsw_idiq_ctrl));
sum3
M_raw = abs(Vs1)*Vn*2/Vdc/k_comp*exp(-1/2*j*Omega_mn*delta*TdTs) - sum3 + 4/pi*Vdrop_pu*exp(-j*PF_angle);
M = abs(M_raw)
%M0 = abs(Vs1)/Vdc_pu/k_comp 
M0 = M*sqrt(3)/2


PF_angle = PF_angle - 1*asin(4/pi*(Vdt_pu + Vdrop_pu)*Vdc/2/6/Lf_0/4/fo/(abs(Is1)*In)) %- 0.2870/5*2*pi*0.7 
ThetaO = angle(Vs1) + angle(M_raw) %+ pi % 2.0788 for n =1;
ThetaC = pi; %pi - 3.5e-7/1e-4*2*pi% pi;


%% vector and matrix creation
Amn_ph = zeros(Mt+1,2*Nt+1);
Bmn_ph = zeros(Mt+1,2*Nt+1);
Vm_mn_ph = zeros(Mt+1,2*Nt+1);

Hac_ph = zeros(Ht,1);
Amn_ll = zeros(Mt+1,2*Nt+1);
Bmn_ll = zeros(Mt+1,2*Nt+1);
Hac_ll = zeros(Ht,1);

%phase to neutral matrices
Theta_a = zeros(Mt+1,2*Nt+1);
Theta_b = zeros(Mt+1,2*Nt+1);
Theta_c = zeros(Mt+1,2*Nt+1);

Cmn_a = zeros(Mt+1,2*Nt+1);
Dmn_a = zeros(Mt+1,2*Nt+1);
Cmn_b = zeros(Mt+1,2*Nt+1);
Dmn_b = zeros(Mt+1,2*Nt+1);
Cmn_c = zeros(Mt+1,2*Nt+1);
Dmn_c = zeros(Mt+1,2*Nt+1);

Vac_a = zeros(Ht+1,1); % cos component
Vac_b = zeros(Ht+1,1);
Vac_c = zeros(Ht+1,1);

Hac_a = zeros(Ht,1); % magnitude component
Kac_a = zeros(Ht,1); % phase angle component
Hac_b = zeros(Ht,1);
Kac_b = zeros(Ht,1);
Hac_c = zeros(Ht,1);
Kac_c = zeros(Ht,1);


%line to line matrices
Theta_ab = zeros(Mt+1,2*Nt+1);
Theta_bc = zeros(Mt+1,2*Nt+1);
Theta_ca = zeros(Mt+1,2*Nt+1);
Cmn_ab = zeros(Mt+1,2*Nt+1);
Dmn_ab = zeros(Mt+1,2*Nt+1);
Cmn_bc = zeros(Mt+1,2*Nt+1);
Dmn_bc = zeros(Mt+1,2*Nt+1);
Cmn_ca = zeros(Mt+1,2*Nt+1);
Dmn_ca = zeros(Mt+1,2*Nt+1);

Cac_ab = zeros(2*Ht+1,1); % cos component
Dac_ab = zeros(2*Ht+1,1); % sin component
Cac_bc = zeros(2*Ht+1,1);
Dac_bc = zeros(2*Ht+1,1);
Cac_ca = zeros(2*Ht+1,1);
Dac_ca = zeros(2*Ht+1,1);

Hac_ab = zeros(Ht,1); % magnitude component
Kac_ab = zeros(Ht,1); % phase angle component
Hac_bc = zeros(Ht,1);
Kac_bc = zeros(Ht,1);
Hac_ca = zeros(Ht,1);
Kac_ca = zeros(Ht,1);

WTHD0_ph = 0;
WTHD0_ll = 0;

WTHD0_a = 0;
WTHD0_b = 0;
WTHD0_c = 0;
WTHD0_ab = 0;
WTHD0_bc = 0;
WTHD0_ca = 0;

THD_isa = 0;
THD_isb = 0;
THD_isc = 0;
THD_ita = 0;
THD_itb = 0;
THD_itc = 0;


%% harmonics calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #1, naturally sampled reference, single-edge carrier PWM

if (modulation_type==1)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if m==0
                if n==1
                    Amn_ph(m+1,n+Nt+1) = M;
                end
            elseif n==0
                Amn_ph(m+1,n+Nt+1) = 2/pi/m*(cos(m*pi)-besselj(n,m*pi*M));
                Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) +(-pi/2);
                Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) +(-pi/2);
                Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) +(-pi/2);
            else
                Amn_ph(m+1,n+Nt+1) = 2/pi/m*besselj(n,m*pi*M);
                Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) +( pi/2 - n*pi/2);
                Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) +( pi/2 - n*pi/2);
                Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) +( pi/2 - n*pi/2);
            end
        end
    end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #2, naturally sampled reference, double-edge carrier PWM
elseif (modulation_type==2)
% phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if m==0
                if n==1
                    Amn_ph(m+1,n+Nt+1) = M;
                end
            else
                Amn_ph(m+1,n+Nt+1) = 4/pi/m*besselj(n,m*pi/2*M)*sin((m+n)*pi/2);
            end
        end
    end


 
  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #3, regularly sampled reference, single-edge carrier PWM

elseif (modulation_type==3)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if m==0
                if n > 0
                    Amn_ph(m+1,n+Nt+1) = 2/pi/q*besselj(n,q*pi*M);
                    Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + pi/2 - n*pi/2;
                    Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + pi/2 - n*pi/2;
                    Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + pi/2 - n*pi/2;
                else
                    Amn_ph(m+1,n+Nt+1) = 0;
                end
            elseif n==0
                Amn_ph(m+1,n+Nt+1) = 2/pi/m*(cos(m*pi)-besselj(n,m*pi*M));
                Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1)-pi/2;
                Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1)-pi/2;
                Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1)-pi/2;
            else
                Amn_ph(m+1,n+Nt+1) = 2/pi/q*besselj(n,q*pi*M);
                Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + pi/2 - n*pi/2;
                Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + pi/2 - n*pi/2;
                Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + pi/2 - n*pi/2;
            end
        end
    end


    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #4, regularly sampled reference, symmetrical double-edge carrier PWM

elseif (modulation_type==4)
    % phase leg Voltage Spectrum


    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if  m == 0
                if n > 0
                     Amn_ph(m+1,n+Nt+1) = 4/pi/q*besselj(n,q*pi/2*M)*sin((q+n)*pi/2);
                else
                    Amn_ph(m+1,n+Nt+1) = 0;
                end
            else
                 Amn_ph(m+1,n+Nt+1) = 4/pi/q*besselj(n,q*pi/2*M)*sin((q+n)*pi/2);
            end
        end
    end
    


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #5, regularly sampled reference, asymmetrical double-edge carrier PWM

elseif (modulation_type==5)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if  m == 0
                if n > 0
                    Amn_ph(m+1,n+Nt+1) = 4/pi/q*besselj(n,q*pi/2*M)*sin((m+n)*pi/2);
                else
                    Amn_ph(m+1,n+Nt+1) = 0;
                end
            else
                Amn_ph(m+1,n+Nt+1) = 4/pi/q*besselj(n,q*pi/2*M)*sin((m+n)*pi/2);
            end
        end
    end

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #6, naturally sampled reference, double-edge carrier SVM

elseif (modulation_type==6)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            if  and(m==0,n==0)
                ;
            elseif and(m==0,n==1)
                Amn_ph(m+1,n+Nt+1) = M;
            elseif and(m==0,~mod(n,3))
                if n>0
                    Amn_ph(m+1,n+Nt+1) = 3*sqrt(3)*M/pi/(n^2-1)*sin(n*pi/6)*sin(n*pi/2);
                end
            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,m*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,m*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:10
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,m*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,m*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,m*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,m*sqrt(3)*pi/4*M));
                    end
                end
                Amn_ph(m+1,n+Nt+1) = 8/m/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);
            elseif and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,m*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,m*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(m*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,m*3*pi/4*M)-besselj(0,m*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:10
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,m*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,m*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,m*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,m*sqrt(3)*pi/4*M));
                    end
                end
                Amn_ph(m+1,n+Nt+1) = 8/m/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);
            end
        end
    end
    Amn_ph(1,Nt+1) = 0;

    


    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #7, regularly sampled reference, symmetrical double-edge carrier SVM

elseif (modulation_type==7)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);

            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Amn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);

            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0;   

            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Amn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);

            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Amn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);
            end
 
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Amn_ph(m+1,n+Nt+1) = abs(Amn_ph(m+1,n+Nt+1));
       
        end
    end
    Amn_ph(1,Nt+1) = 0;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #8, regularly sampled reference, asymmetrical double-edge carrier SVM

elseif (modulation_type==8)
  % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);

            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(m*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0;    
            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(m*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
            end
            Amn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);
            
            
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Amn_ph(m+1,n+Nt+1));
            Amn_ph(m+1,n+Nt+1) = abs(Amn_ph(m+1,n+Nt+1));

        end
        
        
    end
    Amn_ph(1,Nt+1) = 0;

   

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#9, naturally sampled reference, double-edge carrier PWM, w/ deadtime
     elseif (modulation_type==9)
     % phase leg Voltage Spectrum

     for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Omega_mn =2*pi*q;
            
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            
            if m==0
                if n==1
                    Vm_mn_ph(m+1,n+Nt+1) = - M - 8*TdTs/pi*(-1)^(1/2*(n-1))/n*exp(-j*n*PF_angle);
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
                elseif (n>1) && (mod(n, 2) == 1)
                    Vm_mn_ph(m+1,n+Nt+1) = - 8*TdTs/pi*(-1)^(1/2*(n-1))/n*exp(-j*n*PF_angle);
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
                end                
            else            
                sum1=0;
                for k = -Kt:1:Kt
                    if k~=n
                        Emnp = ((-1)^(k-n)-1)*(1+(-1)^(m+k));
                        sum1 = sum1 + exp(j*(k-n)*PF_angle)/pi/(k-n)*besselj(k,1/2*pi*m*M)*Emnp;
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = j^(m+n)/(pi*m*j)*exp(-j*m*pi*delta*TdTs)*((cos(pi*m*TdTs))*besselj(n,1/2*pi*m*M)*(1-(-1)^(m+n)) - sin(pi*m*TdTs)*sum1);
            end
            
            Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
              Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
              Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
              Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));

        end
     end
      Amn_ph(1,Nt+1) = 0;

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #10, regularly sampled reference, symmetrical double-edge carrier PWM, w/ deadtime
    elseif (modulation_type==10)
    % phase leg Voltage Spectrum


    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Omega_mn =2*pi*q;
            
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);

            if  m == 0
                if n > 0
                    Vm_mn_ph(m+1,n+Nt+1) = 2/j/Omega_mn*besselj(n,Omega_mn/4*M)*(cos(1/2*Omega_mn*TdTs))*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*(exp(-1/2*j*Omega_mn)-(-1)^n); 
                    sum1=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Gmnp = (1-(-1)^(k-n))*(exp(-1/2*j*Omega_mn)+(-1)^k);
                            sum1 = sum1 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Gmnp;
                        end
                    end
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) + sum1 - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
                else
                    Vm_mn_ph(m+1,n+Nt+1) = 0;
                end
            else
                    Vm_mn_ph(m+1,n+Nt+1) = 2/j/Omega_mn*besselj(n,Omega_mn/4*M)*(cos(1/2*Omega_mn*TdTs))*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*(exp(-1/2*j*Omega_mn)-(-1)^n);
                    sum1=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Gmnp = (1-(-1)^(k-n))*(exp(-1/2*j*Omega_mn)+(-1)^k);
                            sum1 = sum1 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Gmnp;
                        end
                    end
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) + sum1;
            end
           
             Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
             Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
             Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
             Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
        %    Amn_ph(m+1,n+Nt+1) =   (Vm_mn_ph(m+1,n+Nt+1));
           
        end
        
        
    end
    Amn_ph(1,Nt+1) = 0;
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #11, regularly sampled reference, asymmetrical double-edge carrier PWM, w/ deadtime
    elseif (modulation_type==11)
    % phase leg Voltage Spectrum


    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            q = m + n/p;
            Omega_mn =2*pi*q;
            
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);

            if  m == 0
                if n > 0
                    Vm_mn_ph(m+1,n+Nt+1) = 2/j/Omega_mn*besselj(n,Omega_mn/4*M)*cos(1/2*Omega_mn*TdTs)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*((-1)^m-(-1)^n);
                    sum1=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Fmnp = (1-(-1)^(k-n))*((-1)^m+(-1)^k);
                            sum1 = sum1 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Fmnp;
                        end
                    end
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) + sum1;
                else
                    Vm_mn_ph(m+1,n+Nt+1) = 0;
                end
            else
                    Vm_mn_ph(m+1,n+Nt+1) = 2/j/Omega_mn*besselj(n,Omega_mn/4*M)*cos(1/2*Omega_mn*TdTs)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*((-1)^m-(-1)^n);
                    sum1=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Fmnp = (1-(-1)^(k-n))*((-1)^m+(-1)^k);
                            sum1 = sum1 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Fmnp;
                        end
                    end
                    Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) + sum1;
            end
           
             Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
             Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
             Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
             Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
        %    Amn_ph(m+1,n+Nt+1) =   (Vm_mn_ph(m+1,n+Nt+1));
           
        end
        
        
    end
    Amn_ph(1,Nt+1) = 0;
    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #12, regularly sampled reference, symmetrical double-edge carrier SVM, w/ deadtime

elseif (modulation_type==12)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Omega_mn =2*pi*q;
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            
            sum3=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Gmnp = (1-(-1)^(k-n))*(exp(-1/2*j*Omega_mn)+(-1)^k);
                            sum3 = sum3 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Gmnp;
                        end
                    end
            
        
            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*(cos(1/2*Omega_mn*TdTs))*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
                Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0; 
                Vm_mn_ph(m+1,n+Nt+1) = 0;
                
            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*(cos(1/2*Omega_mn*TdTs))*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
                
            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*(cos(1/2*Omega_mn*TdTs))*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
            end
            
            Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
           
        end
    end
    Amn_ph(1,Nt+1) = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #13, regularly sampled reference, asymmetrical double-edge carrier SVM, w/ deadtime

elseif (modulation_type==13)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Omega_mn =2*pi*q;
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            
            sum3=0;
                    for k = -Kt_1:1:Kt_1
                        if k~=n
                            Fmnp = (1-(-1)^(k-n))*((-1)^m+(-1)^k);
                            sum3 = sum3 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*(PF_angle))/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Fmnp;
                        end
                    end
            
            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(m*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*cos(1/2*Omega_mn*TdTs)*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
                Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0; 
                Vm_mn_ph(m+1,n+Nt+1) = 0;
                
            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(m*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*cos(1/2*Omega_mn*TdTs)*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
                
            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((m+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((m+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((m+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2)*cos(1/2*Omega_mn*TdTs)*exp(-1/2*j*Omega_mn*delta*TdTs) + sum3;
                
            end
            
            Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
           
        end
    end
    Amn_ph(1,Nt+1) = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #14, regularly sampled reference, symmetrical double-edge carrier SVM, w/ deadtime - 2

elseif (modulation_type==14)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Omega_mn =2*pi*q;
           % Kt = 3*ceil(abs(q));
 
           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            %{
            sum3=0;
                    for k = -Kt:1:Kt
                        if k~=n
                            Gmnp = (1-(-1)^(k-n))*(exp(-1/2*j*Omega_mn)+(-1)^k);
                            sum3 = sum3 + 2/j/Omega_mn*besselj(k,Omega_mn/4*M)*sin(1/2*Omega_mn*TdTs)*exp(j*(k-n)*PF_angle)/pi/(k-n)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*Gmnp;
                        end
                    end

           sum4 = 2/j/Omega_mn*besselj(n,Omega_mn/4*M)*(cos(1/2*Omega_mn*TdTs)-1)*exp(-1/2*j*Omega_mn*delta*TdTs)*j^(n)*(exp(-1/2*j*Omega_mn)-(-1)^n);
            %}
           %sum3=sum3*exp(-j*Omega_mn/4);
           %sum4=sum4*exp(-j*Omega_mn/4);
           sum3 = 0;
           sum4 = 0;
            
            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) + sum4 + sum3;
               % Vm_mn_ph(m+1,n+Nt+1) = sum4 + sum3;
               
                Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
                
                
            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0; 
                Vm_mn_ph(m+1,n+Nt+1) = 0;
                
            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) + sum4 + sum3;
                % Vm_mn_ph(m+1,n+Nt+1) = sum4 + sum3;

            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) + sum4 + sum3;
                %Vm_mn_ph(m+1,n+Nt+1) = sum4 + sum3;
            end
            
            Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
           
        end
    end
    Amn_ph(1,Nt+1) = 0;

    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #15, regularly sampled reference, symmetrical double-edge carrier SVM, w/ deadtime - 3

elseif (modulation_type==15)
    % phase leg Voltage Spectrum
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            
            q = m + n/p;
            Omega_mn =2*pi*q;
            Kt = 3*ceil(abs(q));
 %           Kt = 4*ceil(abs(q));
            if Kt < Kt_0
                Kt = Kt_0;
            end
            Theta_a(m+1,n+Nt+1) = m*ThetaC + n*ThetaO;
            Theta_b(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO - 2*pi/3);
            Theta_c(m+1,n+Nt+1) = m*ThetaC + n*(ThetaO + 2*pi/3);
            
                sum3=0;
                for k = -Kt:1:Kt
                    if k~=n
                        Emnp = ((-1)^(k-n)-1)*(1+(-1)^(m+k));
                        sum3 = sum3 + exp(j*(k-n)*PF_angle)/pi/(k-n)*besselj(k,1/2*pi*m*M)*Emnp;
                    end
                end

                sum3 = j^(m+n)/(pi*m*j)*exp(-j*m*pi*delta*TdTs)*(- sin(pi*m*TdTs))*sum3;

           sum4 =  j^(m+n)/(pi*m*j)*exp(-j*m*pi*delta*TdTs)*((cos(pi*m*TdTs)-1)*besselj(n,1/2*pi*m*M)*(1-(-1)^(m+n)));
          
          % sum3=sum3*exp(-j*Omega_mn/4);
          % sum4=sum4*exp(-j*Omega_mn/4);
            
            if  and(m==0,n>0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2);

                Vm_mn_ph(m+1,n+Nt+1) = Vm_mn_ph(m+1,n+Nt+1) - 8*TdTs/pi*sin(n*pi/2)/n*exp(-j*n*PF_angle) - 4/pi/n*Vdrop_pu*sin(n*pi/2)*exp(-j*n*PF_angle);
                
                
            elseif and(m==0,n<=0)
                Amn_ph(m+1,n+Nt+1)=0;
                sum1=0;
                sum2=0; 
                Vm_mn_ph(m+1,n+Nt+1) = 0;
                
            elseif  and(m~=0,n~=0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M)) + ...
                                     1/n*sin(q*pi/2)*cos(n*pi/2)*sin(n*pi/6)*(besselj(0,q*3*pi/4*M)-besselj(0,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) + sum4 + sum3;
                % Vm_mn_ph(m+1,n+Nt+1) = sum4 + sum3;

            elseif and(m~=0,n==0)
                Amn_ph(m+1,n+Nt+1) = pi/6*sin((q+n)*pi/2)*(besselj(n,q*3*pi/4*M)+2*cos(n*pi/6)*besselj(n,q*sqrt(3)*pi/4*M));
                sum1=0;
                sum2=0;
                for k = 1:1:Kt
                    if (k~=-n)
                        sum1 = sum1 + 1/(n+k)*sin((q+k)*pi/2)*cos((n+k)*pi/2)*sin((n+k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n+3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                    if (k~=n)
                        sum2 = sum2 + 1/(n-k)*sin((q+k)*pi/2)*cos((n-k)*pi/2)*sin((n-k)*pi/6)*(besselj(k,q*3*pi/4*M)+2*cos((2*n-3*k)*pi/6)*besselj(k,q*sqrt(3)*pi/4*M));
                    end
                end
                Vm_mn_ph(m+1,n+Nt+1) = 8/q/pi^2*(Amn_ph(m+1,n+Nt+1) + sum1 + sum2) + sum4 + sum3;
                %Vm_mn_ph(m+1,n+Nt+1) = sum4 + sum3;
            end
            
            Amn_ph(m+1,n+Nt+1) = abs(Vm_mn_ph(m+1,n+Nt+1));
            Theta_a(m+1,n+Nt+1) = Theta_a(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_b(m+1,n+Nt+1) = Theta_b(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
            Theta_c(m+1,n+Nt+1) = Theta_c(m+1,n+Nt+1) + angle(Vm_mn_ph(m+1,n+Nt+1));
           
        end
    end
    Amn_ph(1,Nt+1) = 0;

    
end

%% three phase and line-line voltage


    Cmn_a = Amn_ph.*exp(j*Theta_a);
    Cmn_b = Amn_ph.*exp(j*Theta_b);
    Cmn_c = Amn_ph.*exp(j*Theta_c);

    
    for m = 0:1:Mt
        for n = (-Nt):1:Nt
            if and((m*p+n)>=(-Ht),(m*p+n)<0)
                Vac_a(-m*p-n+1) = Vac_a(-m*p-n+1) + conj(Cmn_a(m+1,n+Nt+1)); 
                Vac_b(-m*p-n+1) = Vac_b(-m*p-n+1) + conj(Cmn_b(m+1,n+Nt+1));
                Vac_c(-m*p-n+1) = Vac_c(-m*p-n+1) + conj(Cmn_c(m+1,n+Nt+1));
            elseif and((m*p+n)>=0,(m*p+n)<=Ht)
                Vac_a(m*p+n+1) = Vac_a(m*p+n+1) + Cmn_a(m+1,n+Nt+1); 
                Vac_b(m*p+n+1) = Vac_b(m*p+n+1) + Cmn_b(m+1,n+Nt+1);
                Vac_c(m*p+n+1) = Vac_c(m*p+n+1) + Cmn_c(m+1,n+Nt+1);
            end
        end
    end

 
    Hac_a = abs(Vac_a(2:(Ht+1)));
    Kac_a = angle(Vac_a(2:(Ht+1)));
    Hac_b = abs(Vac_b(2:(Ht+1)));
    Kac_b = angle(Vac_b(2:(Ht+1)));
    Hac_c = abs(Vac_c(2:(Ht+1)));
    Kac_c = angle(Vac_c(2:(Ht+1)));

      
    Vac_ab = Vac_a - Vac_b;
    Vac_bc = Vac_b - Vac_c;
    Vac_ca = Vac_c - Vac_a;
    
    Hac_ab = abs(Vac_ab(2:(Ht+1)));
    Kac_ab = angle(Vac_ab(2:(Ht+1)));
    Hac_bc = abs(Vac_bc(2:(Ht+1)));
    Kac_bc = angle(Vac_bc(2:(Ht+1)));
    Hac_ca = abs(Vac_ca(2:(Ht+1)));
    Kac_ca = angle(Vac_ca(2:(Ht+1)));
    
    %THD and spectrums
    for h = 2:1:Ht
        WTHD0_a = WTHD0_a + Hac_a(h)^2/h^2;
        WTHD0_b = WTHD0_b + Hac_b(h)^2/h^2;
        WTHD0_c = WTHD0_c + Hac_c(h)^2/h^2;
    end
    
    WTHD0_a = sqrt(WTHD0_a)
    WTHD0_b = sqrt(WTHD0_b)
    WTHD0_c = sqrt(WTHD0_c)
    
    figure(3);
    % bar(order,Hac_pu)
    bar(order,abs(Hac_a))
    set(gca,'yscale','log')
    xlabel('Harmonics Number');
    ylabel('Harmonic Magnitude(p.u.)'); 
    title('Single Leg Voltage Spectrum - Hac_a');
    xlim([0 Ht]);
    ylim([1e-4 2]);
    grid on;

    dim = [.65 .62 .3 .3];
    txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['WTHD0=' num2str(WTHD0_a*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);


    Hac_ab_pu = Hac_ab./sqrt(3);
    Hac_bc_pu = Hac_bc./sqrt(3);
    Hac_ca_pu = Hac_ca./sqrt(3);
    
    for h = 2:1:Ht
        WTHD0_ab = WTHD0_ab + Hac_ab_pu(h)^2/h^2;
        WTHD0_bc = WTHD0_bc + Hac_bc_pu(h)^2/h^2;
        WTHD0_ca = WTHD0_ca + Hac_ca_pu(h)^2/h^2;
    end
    
    WTHD0_ab = sqrt(WTHD0_ab)
    WTHD0_bc = sqrt(WTHD0_bc)
    WTHD0_ca = sqrt(WTHD0_ca)
    
    figure(4);
    bar(order,Hac_ab_pu)
 %  bar(order.*fo,abs(Hac_ab_pu.*Vn))
    set(gca,'yscale','log')
    xlabel('Harmonics Number');
    ylabel('Harmonic Magnitude(p.u.)'); 
    title('3-phase Voltage Spectrum - Hac_ab_pu');
    xlim([0 Ht]);
    ylim([1e-4 2]);
    grid on;
  
    dim = [.65 .62 .3 .3];
    txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['WTHD0=' num2str(WTHD0_ab*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);
    
    figure(5);
   % bar(order.*fo,Kac_a*180/pi) %abs(Hac_ab_pu.*Vnll)
   % bar(order.*fo,abs(Hac_a.*Vnph))
    bar(Frequency,abs(Hac_ab_pu.*Vnll))
    set(gca,'yscale','log')
    xlabel('Frequency(Hz)');
    ylabel('Harmonic Magnitude(Vrms)'); 
    title('Line-line AB Voltage Spectrum - Vab');
   % xlim([0 Ht].*fo);
    xlim([0 1e5]);
    ylim([1e-2 1e2]); 
    grid on;
  
    dim = [.65 .62 .3 .3];
    txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['WTHD0=' num2str(WTHD0_ab*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);
    
    figure(6);
   % bar(order.*fo,Kac_b*180/pi) %abs(Hac_bc_pu.*Vnll)
   % bar(order.*fo,abs(Hac_b.*Vnph))
    bar(Frequency,abs(Hac_bc_pu.*Vnll))
    set(gca,'yscale','log')
    xlabel('Frequency(Hz)');
    ylabel('Harmonic Magnitude(Vrms)'); 
    title('Line-line BC Voltage Spectrum - Vbc');
   % xlim([0 Ht].*fo);
    xlim([0 1e5]);
    ylim([1e-2 1e2]);
    grid on;
  
    dim = [.65 .62 .3 .3];
    txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['WTHD0=' num2str(WTHD0_ab*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);
    
    figure(7);
   % bar(order.*fo,Kac_c*180/pi) %abs(Hac_ca_pu.*Vnll)
   % bar(order.*fo,abs(Hac_c.*Vnph))
    bar(Frequency,abs(Hac_ca_pu.*Vnll))
    set(gca,'yscale','log')
    xlabel('Frequency(Hz)');
    ylabel('Harmonic Magnitude(Vrms)'); 
    title('Line-line CA Voltage Spectrum - Vca');
   % xlim([0 Ht].*fo);
    xlim([0 1e5]);
    ylim([1e-2 1e2]);
    grid on;
  
    dim = [.65 .62 .3 .3];
    txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['WTHD0=' num2str(WTHD0_ab*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);
   
    
%%
% phase to neutral voltage waveform reproduction
N_sample = 1e5;
tt = 0:1/fo/N_sample:1/fo;
%tt = 0:1e-7:2.5e-3;

ua2 = zeros(1,length(tt));
ub2 = zeros(1,length(tt));
uc2 = zeros(1,length(tt));
%{
for h = 1:1:Ht
    ua2 = ua2 + sqrt(2)*Hac_ab_pu(h)*Vnll*cos(w(h).*tt+Kac_ab(h));
    ub2 = ub2 + sqrt(2)*Hac_bc_pu(h)*Vnll*cos(w(h).*tt+Kac_bc(h));
    uc2 = uc2 + sqrt(2)*Hac_ca_pu(h)*Vnll*cos(w(h).*tt+Kac_ca(h));
end
%}

Ht_res_start = 1;
Ht_res_end = Ht;

for h = Ht_res_start:1:Ht_res_end
    ua2 = ua2 + Hac_a(h)*Vdc/2*cos(w(h).*tt+Kac_a(h));
    ub2 = ub2 + Hac_b(h)*Vdc/2*cos(w(h).*tt+Kac_b(h));
    uc2 = uc2 + Hac_c(h)*Vdc/2*cos(w(h).*tt+Kac_c(h));
end


%{
Vmag_ab = sqrt(2)*abs(Hac_ab_pu.*Vnll);
Vang_ab = Kac_ab;
Vmag_ab_simu = uori_f(2:end,2);
Vang_ab_simu = uori_f(2:end,3)/180*pi;
%}
Vmag_a = abs(Hac_a.*Vdc/2);
Vang_a = Kac_a;
Vmag_b = abs(Hac_b.*Vdc/2);
Vang_b = Kac_b;
Vmag_c = abs(Hac_c.*Vdc/2);
Vang_c = Kac_c;

    

figure(8);
plot(tt,ua2,'-');%,tt,ub2,'-',tt,uc2,'-');
hold on;
legend('ua\_res','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Voltage(V)'); 
title('Line-neutral An Voltage Waveforms - Van, reproduced and simulated');
grid on;


figure(9);
plot(tt,ub2,'-');%,tt,uc2,'-');
hold on;
legend('ub\_res','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Voltage(V)'); 
title('Line-neutral Bn Voltage Waveforms - Vbn, reproduced and simulated');
grid on;

figure(10);
plot(tt,uc2,'-');%,tt,uc2,'-');
hold on;
legend('uc\_res','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Voltage(V)'); 
title('Line-neutral Cn Voltage Waveforms - Vcn, reproduced and simulated');
grid on;

Va = Vac_a(2:(Ht+1))*Vdc/2;
Vb = Vac_b(2:(Ht+1))*Vdc/2;
Vc = Vac_c(2:(Ht+1))*Vdc/2;


%% source current calculation
% linear circuit
Zt_temp = (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c).*(j*Ls_0.*w + Rs_0)./((j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)+(j*Ls_0.*w + Rs_0));
Zt = Rf.*(j*Lf_0.*w + R_Lf_0 + Zt_temp)./(Rf + j*Lf_0.*w + R_Lf_0 + Zt_temp);
% Zt = j*Lf.*w + Zt_temp;

Zt_temp_1 = (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c).*(j.*Xeq_machine_fit_1 + Req_machine_fit_1 )./((j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)+(j.*Xeq_machine_fit_1 + Req_machine_fit_1));
Zt_temp_2 = (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c).*(j.*Xeq_machine_fit_2 + Req_machine_fit_2 )./((j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)+(j.*Xeq_machine_fit_2 + Req_machine_fit_2));
Zt1 = Rf.*(j*Lf.*w + R_Lf + Zt_temp_1)./(Rf + j*Lf.*w + R_Lf + Zt_temp_1);
Zt2 = Rf.*(j*Lf.*w + R_Lf + Zt_temp_2)./(Rf + j*Lf.*w + R_Lf + Zt_temp_2);

vs1 = 1/3*(Va + Vb.*C120 + Vc.*C240);
vs2 = 1/3*(Va + Vb.*C240 + Vc.*C120);
vsa1 = vs1;
vsa2 = vs2;
vsb1 = vs1.*C240;
vsb2 = vs2.*C120;
vsc1 = vs1.*C120;
vsc2 = vs2.*C240;
%{
vsa1 = 1/3*(Va + Vb.*C120 + Vc.*C240);
vsa2 = 1/3*(Va + Vb.*C240 + Vc.*C120);
vsb1 = 1/3*(Va.*C240 + Vb + Vc.*C120);
vsb2 = 1/3*(Va.*C120 + Vb + Vc.*C240);
vsc1 = 1/3*(Va.*C120 + Vb.*C240 + Vc);
vsc2 = 1/3*(Va.*C240 + Vb.*C120 + Vc);
%}
%{
vsa1 = 1/3*(Va_simu + Vb_simu.*C120 + Vc_simu.*C240);
vsa2 = 1/3*(Va_simu + Vb_simu.*C240 + Vc_simu.*C120);
vsb1 = 1/3*(Va_simu.*C240 + Vb_simu + Vc_simu.*C120);
vsb2 = 1/3*(Va_simu.*C120 + Vb_simu + Vc_simu.*C240);
vsc1 = 1/3*(Va_simu.*C120 + Vb_simu.*C240 + Vc_simu);
vsc2 = 1/3*(Va_simu.*C240 + Vb_simu.*C120 + Vc_simu);
%}


Eq1 = 2*pi*fo*Psi1;
Eq1eq = Eq1*(j.*(w(1).*L_esl) - j./(w(1).*Cf_0) + R_c)/((j.*(w(1).*L_esl) - j./(w(1).*Cf_0) + R_c)+j*Ls_0*w(1) + Rs_0);

switch machine_model
    case 1
        %{
        vsa = vsa1 + vsa2;
        vsb = vsb1 + vsb2;
        vsc = vsc1 + vsc2;
        %}
        vsa = 1/3*(Vac_ab(2:(Ht+1)) - Vac_ca(2:(Ht+1))).*Vdc/2;
        vsb = 1/3*(Vac_bc(2:(Ht+1)) - Vac_ab(2:(Ht+1))).*Vdc/2;
        vsc = 1/3*(Vac_ca(2:(Ht+1)) - Vac_bc(2:(Ht+1))).*Vdc/2;
        isa = vsa./Zt;
        isb = vsb./Zt;
        isc = vsc./Zt;
        isa(1) = (vsa(1) - j*Eq1eq)/Zt(1);
        isb(1) = (vsb(1)  - j*Eq1eq*C240)/Zt(1);
        isc(1) = (vsc(1) - j*Eq1eq*C120)/Zt(1);
    case 2
        is1 = vs1./Zt1;
        is2 = vs2./Zt2;
        isa = is1 + is2;
        isb = is1*C240 + is2*C120;
        isc = is1*C120 + is2*C240;
%         isa = vsa1./Zt1 + vsa2./Zt2;
%         isb = vsb1./Zt1 + vsb2./Zt2;
%         isc = vsc1./Zt1 + vsc2./Zt2;
        isa(1) = (vsa1(1) - j*Eq1eq)/Zt1(1) + vsa2(1)./Zt2(1);
        isb(1) = (vsb1(1)  - j*Eq1eq*C240)/Zt1(1) + vsb2(1)./Zt2(1);
        isc(1) = (vsc1(1) - j*Eq1eq*C120)/Zt1(1) + vsc2(1)./Zt2(1);
    otherwise
        ;
end

isa_abs = abs(isa./sqrt(2));
isb_abs = abs(isb./sqrt(2));
isc_abs = abs(isc./sqrt(2));
isa_ang = angle(isa);
isb_ang = angle(isb);
isc_ang = angle(isc);

    for h = 2:1:Ht
        THD_isa = THD_isa + isa_abs(h)^2;
        THD_isb = THD_isb + isb_abs(h)^2;
        THD_isc = THD_isc + isc_abs(h)^2;
    end
    
    THD_isa = sqrt(THD_isa)/isa_abs(1)
    THD_isb = sqrt(THD_isb)/isb_abs(1)
    THD_isc = sqrt(THD_isc)/isc_abs(1)

%% machine current calculation

% linear circuit
Zeq = (j.*w*Ls_0 + Rs_0) + (j*Lf_0.*w + R_Lf_0) + (j.*w*Ls_0 + Rs_0).*(j*Lf_0.*w + R_Lf_0)./(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c);

Zm = (j.*w*Ls_0 + Rs_0) + (j*Lf_0.*w + R_Lf_0).*(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)./((j*Lf_0.*w + R_Lf_0) + (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c));

% nonlinear circuit - static rotor
%Zeq = (j.*Xs + Rs) + (-j*Zcf./w).*(j*Lf.*w)./(-j*Zcf./w + j*Lf.*w);
Zeq1 = (j.*Xeq_machine_fit_1 + Req_machine_fit_1) + (j*Lf.*w + R_Lf) +  (j.*Xeq_machine_fit_1 + Req_machine_fit_1).*(j*Lf.*w + R_Lf)./(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c);
Zeq2 = (j.*Xeq_machine_fit_2 + Req_machine_fit_2) + (j*Lf.*w + R_Lf) +  (j.*Xeq_machine_fit_2 + Req_machine_fit_2).*(j*Lf.*w + R_Lf)./(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c);

Zm1 =  (j.*Xeq_machine_fit_1 + Req_machine_fit_1)  + (j*Lf.*w + R_Lf).*(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)./((j*Lf.*w + R_Lf) + (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c));
Zm2 =  (j.*Xeq_machine_fit_2 + Req_machine_fit_2)  + (j*Lf.*w + R_Lf).*(j.*(w.*L_esl) - j./(w.*Cf_0) + R_c)./((j*Lf.*w + R_Lf) + (j.*(w.*L_esl) - j./(w.*Cf_0) + R_c));

switch machine_model
    case 1
        ita = vsa./Zeq;
        itb = vsb./Zeq;
        itc = vsc./Zeq;
        ita(1) = ita(1) - (j*Eq1)/Zm(1);
        itb(1) = itb(1) - (j*Eq1*C240)/Zm(1);
        itc(1) = itc(1) - (j*Eq1*C120)/Zm(1);       
        eta = ita.*(j*Ls_0.*w);
        etb = itb.*(j*Ls_0.*w);
        etc = itc.*(j*Ls_0.*w);
      
    case 2
        it1 = vs1./Zeq1;
        it2 = vs2./Zeq2;
        it1(1) = it1(1) - (j*Eq1)/Zm1(1);
        ita = it1 + it2;
        itb = it1*C240 + it2*C120;
        itc = it1*C120 + it2*C240;
    otherwise
        ;
end

        
        
ita_abs = abs(ita./sqrt(2));
itb_abs = abs(itb./sqrt(2));
itc_abs = abs(itc./sqrt(2));
ita_ang = angle(ita);
itb_ang = angle(itb);
itc_ang = angle(itc);

    for h = 2:1:Ht
        THD_ita = THD_ita + ita_abs(h)^2;
        THD_itb = THD_itb + itb_abs(h)^2;
        THD_itc = THD_itc + itc_abs(h)^2;
    end
    
    ita_abs1 = ita_abs(1)
    THD_ita = sqrt(THD_ita)/ita_abs(1)
    itb_abs1 = itb_abs(1)
    THD_itb = sqrt(THD_itb)/itb_abs(1)
    itc_abs1 = itc_abs(1)
    THD_itc = sqrt(THD_itc)/itc_abs(1)
    
% linear circuit model

        vsa = 1/3*(Vac_ab(2:(Ht+1)) - Vac_ca(2:(Ht+1))).*Vdc/2;
        vsb = 1/3*(Vac_bc(2:(Ht+1)) - Vac_ab(2:(Ht+1))).*Vdc/2;
        vsc = 1/3*(Vac_ca(2:(Ht+1)) - Vac_bc(2:(Ht+1))).*Vdc/2;

        ita0 = vsa./Zeq;
        itb0 = vsb./Zeq;
        itc0 = vsc./Zeq;
        ita0(1) = ita0(1) - (j*Eq1)/Zm(1);
        itb0(1) = itb0(1) - (j*Eq1*C240)/Zm(1);
        itc0(1) = itc0(1) - (j*Eq1*C120)/Zm(1);      


        ita0_abs = abs(ita0./sqrt(2));
        itb0_abs = abs(itb0./sqrt(2));
        itc0_abs = abs(itc0./sqrt(2));
        ita0_ang = angle(ita0);
        itb0_ang = angle(itb0);
        itc0_ang = angle(itc0);

figure(11);
% bar(order,Hac_pu)
bar(Frequency,ita_abs)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
title('Machine Phase Current Spectrum Ia');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;
  
dim = [.65 .62 .3 .3];
txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['THD=' num2str(THD_ita*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);

figure(12);
% bar(order,Hac_pu)
bar(Frequency,itb_abs)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
title('Machine Phase Current Spectrum Ib');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);

grid on;
  
dim = [.65 .62 .3 .3];
txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['THD=' num2str(THD_itb*100,3) '%']};  %  num2str(M)],'WTHD0:'};
    annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);
    
figure(13);
% bar(order,Hac_pu)
bar(Frequency,itc_abs)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
title('Machine Phase Current Spectrum Ic');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

dim = [.65 .62 .3 .3];
txt = {['f_c//f_o=' num2str(p) ', M=' num2str(M,3)], ['THD=' num2str(THD_itc*100,3) '%']};  %  num2str(M)],'WTHD0:'};
annotation('textbox',dim,'String',txt,'FitBoxToText','on','FontSize',10);

%% post processing - impedance

% THD_current =[ THD_ita;THD_itb;THD_itc;1/3*(THD_ita+THD_itb+THD_itc)];

Transfer_s =  (1/s/Cf_0 + s*Lf_0 + R_Lf_0 + R_c)/(1/s/Cf_0+R_c)*(s*Ls_0 + Rs_0 + (1/s/Cf_0+R_c)*(s*Lf_0+R_Lf_0)/(1/s/Cf_0+s*Lf_0+R_Lf_0+R_c));
figure(14);
bodemag(Transfer_s,'k-')
fig12 = gcr;
setoptions(fig12,'FreqUnits','Hz');
hold;
Impedance(:,1) =order.*fo; 
Impedance(:,2) = abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf_0 + R_Lf_0)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Ls_0 + Rs_0 + (1./(j.*w)/Cf_0 + R_c).*((j.*w).*Lf_0 + R_Lf_0)./(1./(j.*w)/Cf_0 + R_c +(j.*w).*Lf_0 + R_Lf_0))));
Trans =20.*log10(Impedance(:,2));
loglog(order.*fo,Trans,'r.--');
hold;
hold;
Z1 = 20.*log10(abs(Zeq1));
Z2 = 20.*log10(abs(Zeq2));
loglog(Frequency,Z1,'g.--',Frequency,Z2,'b.--');
legend('Synchronous Rotor - Continuous','Synchronous Rotor - Sampled','Positive','Negative');
xlim([1 5e6]);
title('Source voltage to machine current transfer function - frequency response');

f_res = sqrt((Lf_0+Synchronous_Inductance+L_filter2)/(Lf_0*(Synchronous_Inductance+L_filter2)*Cf_0))/2/pi

figure(15);

Trans = 20.*log10(abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit + Req_machine_fit + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf)))));
Z1 = 20.*log10(abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit_1 + Req_machine_fit_1 + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf)))));
Z2 = 20.*log10(abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit_2 + Req_machine_fit_2 + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf)))));

semilogx(Frequency,Trans,'b.--');
hold;
semilogx(Frequency,Z1,'k.--',Frequency,Z2,'r.--');
legend('Stationary Rotor','Positive Sequence','Negative Sequence');
xlim([1e2 1e6]);
xlabel('Frequency(Hz)');
ylabel('Magnitude(dB)'); 
title('Source Voltage to Machine Current Impedance vs Frequency');
grid on

%{
figure(16);

Lth0 = abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit + Req_machine_fit + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf))))./w;
Lth1 = abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit_1 + Req_machine_fit_1 + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf))))./w;
Lth2 = abs(((1./(j.*w)/Cf_0 + R_c + (j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0 + R_c).*((j.*w).*Leq_machine_fit_2 + Req_machine_fit_2 + (1./(j.*w)/Cf_0+ R_c).*((j.*w).*Lf + R_Lf)./(1./(j.*w)/Cf_0+ R_c +(j.*w).*Lf + R_Lf))))./w;

semilogx(Frequency,Lth0,'b.--');
hold;
semilogx(Frequency,Lth1,'k.--',Frequency,Lth2,'r.--');
legend('Stationary Rotor','Positive Sequence','Negative Sequence');
xlim([1e2 1e6]);
xlabel('Frequency(Hz)');
ylabel('Magnitude(dB)'); 
title('Source Voltage to Machine Current Impedance vs Frequency');
grid on
%}

Zratio1 = abs(Zeq1)./abs(Zeq);
Zratio2 = abs(Zeq2)./abs(Zeq);
Zm_ratio1 = abs((Req_machine_fit_1 + j*Xeq_machine_fit_1)./(Rs_0 + j.*Xs_0));
Zm_ratio2 = abs((Req_machine_fit_2 + j*Xeq_machine_fit_2)./(Rs_0 + j.*Xs_0));
figure(16);
plot(Frequency,Zratio1,Frequency,Zratio2,Frequency,Zm_ratio1,Frequency,Zm_ratio2);
set(gca,'xscale','log')
grid on

figure(17)
loglog(Frequency, Impedance(:,2),'r.--', Frequency,abs(Zeq1),'g.--',Frequency,abs(Zeq2),'b.--');

% pm loss
if machine_model == 2
tic
S11=abs(it1).^2.*(-w.*imag(Ls_complex_1)).*slip1;
S22=abs(it2).^2.*(-w.*imag(Ls_complex_1)).*slip2;
Ppm = 1.5*(sum(S11+S22)-S11(1));
toc
end
% 
% Imag_a_simu = iori_f(2:end,2);
% Iang_a_simu = iori_f(2:end,3)/180*pi;
% Imag_b_simu = iori_f(2:end,4);
% Iang_b_simu = iori_f(2:end,5)/180*pi;
% Imag_c_simu = iori_f(2:end,6);
% Iang_c_simu = iori_f(2:end,7)/180*pi;
% 
% Ita_Fundamental_simu=Imag_a_simu(1);
% Ita_Harmonic_simu=sqrt(sum(Imag_a_simu(1:length(Imag_a_simu)-1).^2));
% THD_ita_simu = sqrt(Ita_Harmonic_simu^2-Ita_Fundamental_simu^2)/Ita_Fundamental_simu
% Itb_Fundamental_simu=Imag_b_simu(1);
% Itb_Harmonic_simu=sqrt(sum(Imag_b_simu(1:length(Imag_b_simu)-1).^2));
% THD_itb_simu = sqrt(Itb_Harmonic_simu^2-Itb_Fundamental_simu^2)/Itb_Fundamental_simu
% Itc_Fundamental_simu=Imag_c_simu(1);
% Itc_Harmonic_simu=sqrt(sum(Imag_c_simu(1:length(Imag_c_simu)-1).^2));
% THD_itc_simu = sqrt(Itc_Harmonic_simu^2-Itc_Fundamental_simu^2)/Itc_Fundamental_simu
% 
% THD_current =[THD_ita, THD_ita_simu;THD_itb, THD_itb_simu;THD_itc, THD_itc_simu; ...
%     1/3*(THD_ita+THD_itb+THD_itc), 1/3*(THD_ita_simu+THD_itb_simu+THD_itc_simu)];

% ita_simu = Imag_a_simu.*exp(j.*Iang_a_simu);
% itb_simu = Imag_b_simu.*exp(j.*Iang_b_simu);
% itc_simu = Imag_c_simu.*exp(j.*Iang_c_simu);
% 
% Idelta_a=[ita(1:Rank),ita_simu(1:Rank),ita(1:Rank)-ita_simu(1:Rank)];
% Idelta_b=[itb(1:Rank),itb_simu(1:Rank),itb(1:Rank)-itb_simu(1:Rank)];
% Idelta_c=[itc(1:Rank),itc_simu(1:Rank),itc(1:Rank)-itc_simu(1:Rank)];

%% test results

% use single excel file
data_number = 'tek0000';
filename = ['TestData\' data_number '.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  data_number;
opts.SelectedVariableNames = [1:5]; 
opts.DataRange = '22:1000021';
TestData = readmatrix(filename,opts);

opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '9:9';
Tss = readmatrix(filename,opts);

%{
% using seperate excel files
filename = ['Tek_CH1_Wfm.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  'Tek_CH1_Wfm';
opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '22:1000021';
Test_ia = readmatrix(filename,opts);

filename = ['Tek_CH2_Wfm.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  'Tek_CH2_Wfm';
opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '22:1000021';
Test_ib = readmatrix(filename,opts);

filename = ['Tek_CH3_Wfm.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  'Tek_CH3_Wfm';
opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '22:1000021';
Test_vab  = readmatrix(filename,opts);

filename = ['Tek_CH4_Wfm.xlsx'];
opts = detectImportOptions(filename);
opts.Sheet =  'Tek_CH4_Wfm';
opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '22:1000021';
Test_ic = readmatrix(filename,opts);

opts.SelectedVariableNames = [1:1]; 
Test_time = readmatrix(filename,opts);

opts.SelectedVariableNames = [2:2]; 
opts.DataRange = '9:9';
Tss = readmatrix(filename,opts);
% fft window
% n_pu = 1;
% fn = 500;
MM=10;           % # of fundamental cycles
% f=n_pu*fn;
T=1/fo;          % time duration of one fundamental cycle

N=floor(T/Tss);  % # of sampling points during one cycle
ff = 1/(MM*N*Tss);

Time = Test_time(1:(N*MM));
Current = [Test_ia(1:(N*MM)),Test_ib(1:(N*MM)),Test_ic(1:(N*MM))];
Voltage = Test_vab(1:(N*MM));
%}

%%
% fft window
% n_pu = 1;
% fn = 500;
MM=2;           % # of fundamental cycles
% f=n_pu*fn;
T=1/fo;          % time duration of one fundamental cycle

N=floor(T/Tss);  % # of sampling points during one cycle
ff = 1/(MM*N*Tss);

Time = TestData(1:(N*MM),1);
Current = TestData(1:(N*MM),[2,3,5]);
Voltage = TestData(1:(N*MM),4);


% voltage fft
x = Voltage(1:(N*MM));
% x = Phase_Voltage_2.signals(1).values(T_start:(T_end-1));
%x = LLVoltage_1.signals(1).values(T_start:(T_end-1));
X=fft(x);
C=X(1:1:floor((N*MM)/2))/(N*MM);
FrequencyT2 = (0:length(C)-1)'*(ff/MM); 

MM_Sample = ((1+MM):MM:floor(N*MM/2));
FrequencyT = FrequencyT2(MM_Sample);


Voltage_Mag_A = abs(C(MM_Sample))*2;   
Voltage_Angle_A = angle(C(MM_Sample))*(180/pi);

Voltage_DC = abs(C(1))*2;
Van_Fundamental=Voltage_Mag_A(1)/sqrt(2)
Van_Harmonic=sqrt(sum(Voltage_Mag_A(1:length(Voltage_Mag_A)-1).^2))/sqrt(2)
Van_THD = sqrt(Van_Harmonic^2-Van_Fundamental^2)/Van_Fundamental

figure(25);
bar(FrequencyT,Voltage_Mag_A/sqrt(2))
set(gca,'yscale','log');
%set(gca,'xscale','log','yscale','log')
xlabel('Frequency (Hz)');
ylabel('Voltage (Volt,rms)'); 
title('Voltage FFT');
xlim([-1 2e4]);
ylim([1e-2 1e2]);
grid on;

figure(26);
plot(Time,Voltage,'-');
legend('uab','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Voltage(V)'); 
title('Line-line AB Voltage Waveforms - Vab');
grid on;


% current fft
x = Current(1:N*MM,1);
% x = Phase_Voltage_2.signals(1).values(T_start:(T_end-1));
%x = LLVoltage_1.signals(1).values(T_start:(T_end-1));
X=fft(x);
C=X(1:1:floor((N*MM)/2))/(N*MM);

Current_Mag_A = abs(C(MM_Sample))*2;   
Current_Angle_A = angle(C(MM_Sample))*180/pi;

Ia_Fundamental=Current_Mag_A(1)/sqrt(2)
Ia_Harmonic=sqrt(sum(Current_Mag_A(1:length(Current_Mag_A)-1).^2))/sqrt(2)
Ia_THD = sqrt(Ia_Harmonic^2-Ia_Fundamental^2)/Ia_Fundamental

% phase b
x = Current(1:N*MM,2);
X=fft(x);
C=X(1:1:floor((N*MM)/2))/(N*MM);

Current_Mag_B = abs(C(MM_Sample))*2;   
Current_Angle_B = angle(C(MM_Sample))*180/pi;

Ib_Fundamental=Current_Mag_B(1)/sqrt(2)
Ib_Harmonic=sqrt(sum(Current_Mag_B(1:length(Current_Mag_B)-1).^2))/sqrt(2)
Ib_THD = sqrt(Ib_Harmonic^2-Ib_Fundamental^2)/Ib_Fundamental

% phase c
x = Current(1:N*MM,3);
X=fft(x);
C=X(1:1:floor((N*MM)/2))/(N*MM);

Current_Mag_C = abs(C(MM_Sample))*2;   
Current_Angle_C = angle(C(MM_Sample))*180/pi;

Ic_Fundamental=Current_Mag_C(1)/sqrt(2)
Ic_Harmonic=sqrt(sum(Current_Mag_C(1:length(Current_Mag_C)-1).^2))/sqrt(2)
Ic_THD = sqrt(Ic_Harmonic^2-Ic_Fundamental^2)/Ic_Fundamental



figure(27);
bar(FrequencyT,Current_Mag_A/sqrt(2))
set(gca,'yscale','log');
%set(gca,'xscale','log','yscale','log')
xlabel('Frequency (Hz)');
ylabel('Test Current A(Amp,rms)'); 
title('Current FFT');
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(28);
bar(FrequencyT,Current_Mag_B/sqrt(2))
set(gca,'yscale','log');
%set(gca,'xscale','log','yscale','log')
xlabel('Frequency (Hz)');
ylabel('Test Current B(Amp,rms)'); 
title('Current FFT');
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(29);
bar(FrequencyT,Current_Mag_C/sqrt(2))
set(gca,'yscale','log');
%set(gca,'xscale','log','yscale','log')
xlabel('Frequency (Hz)');
ylabel('Test Current C(Amp,rms)'); 
title('Current FFT');
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(30);
plot(Time(1:N),Current(1:N,1),'-',Time(1:N),Current(1:N,2),'-',Time(1:N),Current(1:N,3),'-');
legend('ia','ib','ic','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Current(A)'); 
title('Test Phase Current Waveforms - iabc');
grid on;



%% compare analytical and test results
%time domain data
ia2 = zeros(1,length(tt));
ib2 = zeros(1,length(tt));
ic2 = zeros(1,length(tt));

%Time_shift =0.001*(0.0109-1.1218-0.007-1.2298+0.1+0.05-0.2);
% match PWM voltage waveform:(0.0109-1.1218-0.007-0.085*0 - 0.0053)*0.001;

Time_shift =-0.001*(0.5791+0.06+0.03+0.22+0.2+0.2523 + 3.7502+0.4731 + 0.02);

for h = 1:1:Ht
    ia2 = ia2 + sqrt(2)*ita_abs(h).*cos(w(h).*(tt+Time_shift)+ita_ang(h));
    ib2 = ib2 + sqrt(2)*itb_abs(h).*cos(w(h).*(tt+Time_shift)+itb_ang(h));
    ic2 = ic2 + sqrt(2)*itc_abs(h).*cos(w(h).*(tt+Time_shift)+itc_ang(h));
end

uab2 = zeros(1,length(tt));
ubc2 = zeros(1,length(tt));
uca2 = zeros(1,length(tt));
for h = Ht_res_start:1:Ht_res_end
    uab2 = uab2 + Hac_a(h)*Vdc/2*cos(w(h).*(tt+Time_shift)+Kac_a(h));
    ubc2 = ubc2 + Hac_b(h)*Vdc/2*cos(w(h).*(tt+Time_shift)+Kac_b(h));
    uca2 = uca2 + Hac_c(h)*Vdc/2*cos(w(h).*(tt+Time_shift)+Kac_c(h));
end



order_min = min(length(isa_abs),(length(Current_Mag_A)-1));
if length(isa_abs)<length(Current_Mag_A)
    Frequency_min = Frequency(1:order_min);
else
    Frequency_min = FrequencyT(2:(order_min+1));
end

% frequency domain data, without power analyzer data, use scope FFT data
Test_Harmonic = (Current_Mag_A(1:order_min) + Current_Mag_B(1:order_min) + Current_Mag_C(1:order_min))/sqrt(2)/3;
%Test_Harmonic = (Current_Mag_B(1:order_min) + Current_Mag_C(1:order_min))/sqrt(2)/2;
Test_THD = 1/3*(Ia_THD + Ib_THD + Ic_THD);

Proposed_Harmonic = (ita_abs(1:order_min) + itb_abs(1:order_min) +itc_abs(1:order_min))/3;
Existing_Harmonic = (ita0_abs(1:order_min) + itb0_abs(1:order_min) +itc0_abs(1:order_min))/3;

CurrentA_compare = [ita_abs(1:order_min),Current_Mag_A(1:order_min)/sqrt(2)];
CurrentB_compare = [itb_abs(1:order_min),Current_Mag_B(1:order_min)/sqrt(2)];
CurrentC_compare = [itc_abs(1:order_min),Current_Mag_C(1:order_min)/sqrt(2)];



figure(31);
plot(1e3*(Time(1:N)+0.05),Current((1:N)+0*N,1),'*-',1e3*(Time(1:N)+0.05),Current((1:N)+0*N,2),'*-',1e3*(Time(1:N)+0.05),Current((1:N)+0*N,3),'*-','LineWidth',1,'MarkerSize',3);
%plot(1e3*(Time(1:N)+0.02),-Current((1:N)+5*N,2),'*-',1e3*(Time(1:N)+0.02),-Current((1:N)+5*N,1),'*-',1e3*(Time(1:N)+0.02),-Current((1:N)+5*N,3),'*-','LineWidth',1,'MarkerSize',1.5);
hold on;
plot(1e3*tt,ia2,'-.',1e3*tt,ib2,'-.',1e3*tt,ic2,'-.','LineWidth',3);
%  lgd = legend('$i_a$ (meas.)','$i_b$ (meas.)','$i_c$ (meas.)','$i_a$ (pred.)','$i_b$ (pred.)','$i_c$ (pred.)','Location', 'northoutside','FontSize',33, 'Interpreter','latex');
%  lgd.Orientation = 'horizontal';
%  lgd.NumColumns = 6;
%  legend('boxoff');
xlabel('Time [ms]','FontSize',40);
ylabel('Current [A]','FontSize',40); 
%title('Machine Phase Current Waveforms Iabc - Analytical');
ylim([-8,8]);
yticks([-8 -4 0 4 8]);
%set(gcf,'Position',[200 200 850 640]);
%set(gcf,'Position',[100 100 800 600]);
%set(gcf,'Position',[300 300 1600 700]);
set(gcf,'Position',[300 300 1600 600]);
%xticks([0 1 2 3 4]);
set(gca,'GridLineStyle','-.');
set(gca,'FontSize',36);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
%title('Machine Phase Current Waveforms Iabc - Test and Prediction');
f = gcf;
filename = ['CurrentWaveform_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(32);
Time2 = Time(1:N) + 5e-2;
Voltage2 = Voltage(1:N);
%Voltage2 = [Voltage2((5000-308):5000);Voltage2(1:(5000-308-1))];
plot(1e3*tt,uca2,'-','Linewidth',2);
hold on;
plot(1e3*Time2,Voltage2,'-','Linewidth',2);
legend('Predicted','Measured','Location', 'southeast','Fontsize',40);
yticks([-60 -40 -20 0 20 40 60])
legend('boxoff');
xlabel('Time [ms]','Fontsize',40);  
ylabel('Voltage [V]','Fontsize',40); 
%title('Line-Line Voltage Waveforms - Vbc, reproduced and measured');
set(gcf,'Position',[300 300 1600 600]);
%xticks([0 1 2 3 4]);
set(gca,'GridLineStyle','-.');
set(gca,'FontSize',36);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% f = gcf;
% filename = ['VoltageWaveform_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
% exportgraphics(f,filename,'Resolution',300);
%print(f,filename,'-dpng','-r600'); 

figure(35);
% bar(order,Hac_pu)
bar(Frequency_min,CurrentA_compare)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
legend('Analytical','Measurement');
title('Machine Phase Current Spectrum Ia');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(36);
% bar(order,Hac_pu)
bar(Frequency_min,CurrentB_compare)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
legend('Analytical','Measurement');
title('Machine Phase Current Spectrum Ib');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(37);
% bar(order,Hac_pu)
bar(Frequency_min,CurrentC_compare)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Magnitude(A)'); 
legend('Analytical','Measurement');
title('Machine Phase Current Spectrum Ic');
% xlim([0 Ht].*fo);
xlim([-1 2e4]);
ylim([1e-5 1e2]);
grid on;

figure(38);
% bar(order,Hac_pu)
bar(Frequency_min,(CurrentA_compare + CurrentB_compare + CurrentC_compare)/3)
set(gca,'yscale','log')
xlabel('Frequency(Hz)');
ylabel('Harmonic Amplitude(A)'); 
legend('Analytical','Measurement');
title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([-1 5e4]);
ylim([1e-3 1e2]);
grid on;


figure(39);
plot(1e3*tt,ib2,'m',1e3*tt,ic2,'b','LineWidth',2);
hold on;
plot(1e3*tt,ia2,'color',[0, 0.75, 0.75],'LineWidth',2);
lgd = legend('i_a','i_b','i_c','Location', 'south','FontSize',30);
lgd.NumColumns = 3;
legend('boxoff');
xlabel('Time [ms]','FontSize',20);
ylabel('Current [A]','FontSize',20); 
%title('Machine Phase Current Waveforms Iabc - Analytical');
ylim([-7,7]);
set(gcf,'Position',[200 200 1050 800]);%[200 200 850 640]);
%xticks([0 1 2 3 4]);
set(gca,'GridLineStyle','-.');
set(gca,'FontSize',20);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

% f = gcf;
% filename = ['AnalyticalCurrentWaveform_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
%exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

order_shown =200;
figure(40);
% bar(order,Hac_pu)
Test_order = (1:order_min)';
Test_analytical_Harmonic = [Proposed_Harmonic(1:order_min), Test_Harmonic(1:order_min)];
bar(Test_order,Test_analytical_Harmonic)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',40);
ylabel('Amplitude [A]','FontSize',40); 

lgd = legend('Predicted','Measured','FontSize',40,'Location', 'northoutside');
lgd.Orientation = 'horizontal';
lgd.NumColumns = 2;
legend('boxoff');

% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
%yticks([0 0.1 0.2])
xticks([0 40 80 120 160 200])
xlim([0 order_shown]);
ylim([0 0.2]);
%set(gcf,'Position',[300 300 800 600]);
%set(gcf,'Position',[300 300 1600 700]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',36);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['CurrentSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(41);
% bar(order,Hac_pu)
Test_order = (1:order_min)';
Test_analytical_Harmonic = [abs(Proposed_Harmonic(1:order_min))-abs(Test_Harmonic(1:order_min))];
bar(Test_order,Test_analytical_Harmonic)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Amplitude (A)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
title('Machine Phase Current Spectrum Difference');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([-0.04 0.02]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));



figure(42);
% bar(order,Hac_pu)
Voltage_order = (1:order_min)';
%Test_analytical_Vll_Harmonic = [(Hac_ab_pu(1:order_min) + Hac_bc_pu(1:order_min) +Hac_ca_pu(1:order_min))*Vnll/3, Voltage_Mag_A(1:order_min)/sqrt(2)];
%bar(Voltage_order,Test_analytical_Vll_Harmonic)
Test_analytical_Vph_Harmonic = [(Hac_a(1:order_min) + Hac_b(1:order_min) +Hac_c(1:order_min))*Vnph/3, Voltage_Mag_A(1:order_min)/sqrt(2)];
bar(Voltage_order,Test_analytical_Vph_Harmonic)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',40);
ylabel('Amplitude [V]','FontSize',40); 
legend('Predicted','Measured','FontSize',40);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-2 30]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',36);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
% f = gcf;
% filename = ['VoltageSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
% exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(43);
%Vll_Harmonic_difference = [abs((Hac_ab_pu(1:order_min) + Hac_bc_pu(1:order_min) +Hac_ca_pu(1:order_min))*Vnll/3)-abs(Voltage_Mag_A(1:order_min)/sqrt(2))];
Vll_Harmonic_difference = [abs((Hac_a(1:order_min) + Hac_b(1:order_min) +Hac_c(1:order_min))*Vnph/3)-abs(Voltage_Mag_A(1:order_min)/sqrt(2))];
bar(Voltage_order,Vll_Harmonic_difference)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (Vrms)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
title('Machine Voltage Spectrum Difference');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([-3 3]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));


figure(44);
% bar(order,Hac_pu)
Vll_Harmonic_difference = [abs((Hac_a(1:order_min) + Hac_b(1:order_min) +Hac_c(1:order_min))*Vnph/3)-abs(Voltage_Mag_A(1:order_min)/sqrt(2))];
Test_analytical_Vll_Harmonic_R = [Vll_Harmonic_difference./Impedance(1:order_min,2)];
bar(Voltage_order,Test_analytical_Vll_Harmonic_R)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (Vrms)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
title('Machine Voltage//Impedance Spectrum Difference');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([-0.08 0.05]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

figure(45);
% bar(order,Hac_pu)
Test_order = (1:order_min)';
Test_analytical_Harmonic = [Test_Harmonic(1:order_min),...
                            Proposed_Harmonic(1:order_min), ...
                            Existing_Harmonic(1:order_min)];
bar(Test_order,Test_analytical_Harmonic)
%set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',40);
ylabel('Amplitude [A]','FontSize',40); 

lgd = legend('Measured','Proposed','Existing','FontSize',40,'Location', 'northoutside');
lgd.Orientation = 'horizontal';
lgd.NumColumns = 3;
legend('boxoff');

% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
%yticks([0 0.1 0.2])
xticks([0 40 80 120 160 200])
xlim([0 order_shown]);
ylim([0 0.2]);
%set(gcf,'Position',[300 300 800 600]);
%set(gcf,'Position',[300 300 1600 700]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',36);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['CurrentSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

%% THD calculation

%band THD
%THD0
THD0_a_analytical = sqrt(sum(ita_abs(2:round(p/2)).^2));
THD0_a_test = sqrt(sum(Current_Mag_A(2:round(p/2)).^2))/sqrt(2);
THD0_b_analytical = sqrt(sum(itb_abs(2:round(p/2)).^2));
THD0_b_test = sqrt(sum(Current_Mag_B(2:round(p/2)).^2))/sqrt(2);
THD0_c_analytical = sqrt(sum(itc_abs(2:round(p/2)).^2));
THD0_c_test = sqrt(sum(Current_Mag_C(2:round(p/2)).^2))/sqrt(2);
THD0_analytical = (THD0_a_analytical+THD0_b_analytical+THD0_c_analytical)/3;
THD0_test = (THD0_a_test+THD0_b_test+THD0_c_test)/3;

%THD1
THD1_a_analytical = sqrt(sum(ita_abs(round(p/2+1):round(1.5*p)).^2));
THD1_a_test = sqrt(sum(Current_Mag_A(round(p/2+1):round(1.5*p)).^2))/sqrt(2);
THD1_b_analytical = sqrt(sum(itb_abs(round(p/2+1):round(1.5*p)).^2));
THD1_b_test = sqrt(sum(Current_Mag_B(round(p/2+1):round(1.5*p)).^2))/sqrt(2);
THD1_c_analytical = sqrt(sum(itc_abs(round(p/2+1):round(1.5*p)).^2));
THD1_c_test = sqrt(sum(Current_Mag_C(round(p/2+1):round(1.5*p)).^2))/sqrt(2);
THD1_analytical = (THD1_a_analytical+THD1_b_analytical+THD1_c_analytical)/3;
THD1_test = (THD1_a_test+THD1_b_test+THD1_c_test)/3;

%THD2
THD2_a_analytical = sqrt(sum(ita_abs(round(1.5*p+1):round(2.5*p)).^2));
THD2_a_test = sqrt(sum(Current_Mag_A(round(1.5*p+1):round(2.5*p)).^2))/sqrt(2);
THD2_b_analytical = sqrt(sum(itb_abs(round(1.5*p+1):round(2.5*p)).^2));
THD2_b_test = sqrt(sum(Current_Mag_B(round(1.5*p+1):round(2.5*p)).^2))/sqrt(2);
THD2_c_analytical = sqrt(sum(itc_abs(round(1.5*p+1):round(2.5*p)).^2));
THD2_c_test = sqrt(sum(Current_Mag_C(round(1.5*p+1):round(2.5*p)).^2))/sqrt(2);
THD2_analytical = (THD2_a_analytical+THD2_b_analytical+THD2_c_analytical)/3;
THD2_test = (THD2_a_test+THD2_b_test+THD2_c_test)/3;

%THD3
THD3_a_analytical = sqrt(sum(ita_abs(round(2.5*p+1):round(3.5*p)).^2));
THD3_a_test = sqrt(sum(Current_Mag_A(round(2.5*p+1):round(3.5*p)).^2))/sqrt(2);
THD3_b_analytical = sqrt(sum(itb_abs(round(2.5*p+1):round(3.5*p)).^2));
THD3_b_test = sqrt(sum(Current_Mag_B(round(2.5*p+1):round(3.5*p)).^2))/sqrt(2);
THD3_c_analytical = sqrt(sum(itc_abs(round(2.5*p+1):round(3.5*p)).^2));
THD3_c_test = sqrt(sum(Current_Mag_C(round(2.5*p+1):round(3.5*p)).^2))/sqrt(2);
THD3_analytical = (THD3_a_analytical+THD3_b_analytical+THD3_c_analytical)/3;
THD3_test = (THD3_a_test+THD3_b_test+THD3_c_test)/3;

%THD4
THD4_a_analytical = sqrt(sum(ita_abs(round(3.5*p+1):round(4.5*p)).^2));
THD4_a_test = sqrt(sum(Current_Mag_A(round(3.5*p+1):round(4.5*p)).^2))/sqrt(2);
THD4_b_analytical = sqrt(sum(itb_abs(round(3.5*p+1):round(4.5*p)).^2));
THD4_b_test = sqrt(sum(Current_Mag_B(round(3.5*p+1):round(4.5*p)).^2))/sqrt(2);
THD4_c_analytical = sqrt(sum(itc_abs(round(3.5*p+1):round(4.5*p)).^2));
THD4_c_test = sqrt(sum(Current_Mag_C(round(3.5*p+1):round(4.5*p)).^2))/sqrt(2);
THD4_analytical = (THD4_a_analytical+THD4_b_analytical+THD4_c_analytical)/3;
THD4_test = (THD4_a_test+THD4_b_test+THD4_c_test)/3;

%THD5
THD5_a_analytical = sqrt(sum(ita_abs(round(4.5*p+1):Ht).^2));
THD5_a_test = sqrt(sum(Current_Mag_A(round(4.5*p+1):Ht).^2))/sqrt(2);
THD5_b_analytical = sqrt(sum(itb_abs(round(4.5*p+1):Ht).^2));
THD5_b_test = sqrt(sum(Current_Mag_B(round(4.5*p+1):Ht).^2))/sqrt(2);
THD5_c_analytical = sqrt(sum(itc_abs(round(4.5*p+1):Ht).^2));
THD5_c_test = sqrt(sum(Current_Mag_C(round(4.5*p+1):Ht).^2))/sqrt(2);
THD5_analytical = (THD5_a_analytical+THD5_b_analytical+THD5_c_analytical)/3;
THD5_test = (THD5_a_test+THD5_b_test+THD5_c_test)/3;



%Total THD
It1_rms = [ita_abs1, Ia_Fundamental;itb_abs1, Ib_Fundamental; itc_abs1, Ic_Fundamental;...
    1/3*(ita_abs1 + itb_abs1 +itc_abs1), 1/3*(Ia_Fundamental + Ib_Fundamental +Ic_Fundamental)];


THD_current =[THD_ita, Ia_THD;THD_itb, Ib_THD;THD_itc, Ic_THD; ...
    1/3*(THD_ita+THD_itb+THD_itc), 1/3*(Ia_THD + Ib_THD + Ic_THD)];

Summary = [It1_rms,THD_current];

Summary(5,1)=Summary(4,1);
Summary(5,2)=Test_Harmonic(1);
Summary(5,3)=Summary(4,3);
Summary(5,4)=Test_THD;
Summary(5,5:16) = [THD0_analytical, THD0_test, ...
    THD1_analytical, THD1_test, ...
    THD2_analytical, THD2_test, ...
    THD3_analytical, THD3_test, ...
    THD4_analytical, THD4_test, ...
    THD5_analytical, THD5_test,];
filename = ['THDSummary_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A_M' num2str(M,'%4.2f') '.xlsx'];
Wavefrom_V_2 = array2table(Summary,...  
    'VariableNames',{'I1_analytical', 'I1_measurement', 'THD_analytical', 'THD_measurement',...
    'THD0_analytical', 'THD0_measurement', 'THD1_analytical', 'THD1_measurement',...
    'THD2_analytical', 'THD2_measurement', 'THD3_analytical', 'THD3_measurement',...
    'THD4_analytical', 'THD4_measurement', 'THD5_analytical', 'THD5_measurement',...
    });
writetable(Wavefrom_V_2,filename,'Sheet',1,'Range','A1');

Operation_Point = [Vs1_ll_rms Vt1_ll_rms Eq1_ll_rms M0 Test_analytical_Vph_Harmonic(1,:)];

Comparison_Table = ...
[
1,              Test_Harmonic(1),  Proposed_Harmonic(1),  Existing_Harmonic(1); 
p-4,            Test_Harmonic(p-4),  Proposed_Harmonic(p-4),  Existing_Harmonic(p-4); 
p-2,            Test_Harmonic(p-2),  Proposed_Harmonic(p-2),  Existing_Harmonic(p-2); 
p-1,            Test_Harmonic(p-1),  Proposed_Harmonic(p-1),  Existing_Harmonic(p-1); 
p+1,            Test_Harmonic(p+1),  Proposed_Harmonic(p+1),  Existing_Harmonic(p+1); 
p+2,            Test_Harmonic(p+2),  Proposed_Harmonic(p+2),  Existing_Harmonic(p+2); 
p+4,            Test_Harmonic(p+4),  Proposed_Harmonic(p+4),  Existing_Harmonic(p+4); 
%2*p-5,          Test_Harmonic(2*p-5),  Proposed_Harmonic(2*p-5),  Existing_Harmonic(2*p-5); 
2*p-1,          Test_Harmonic(2*p-1),  Proposed_Harmonic(2*p-1),  Existing_Harmonic(2*p-1); 
2*p+1,          Test_Harmonic(2*p+1),  Proposed_Harmonic(2*p+1),  Existing_Harmonic(2*p+1);
%2*p+5,          Test_Harmonic(2*p+5),  Proposed_Harmonic(2*p+5),  Existing_Harmonic(2*p+5);  
3*p-4,            Test_Harmonic(3*p-4),  Proposed_Harmonic(3*p-4),  Existing_Harmonic(3*p-4); 
3*p-2,            Test_Harmonic(3*p-2),  Proposed_Harmonic(3*p-2),  Existing_Harmonic(3*p-2); 
3*p+2,            Test_Harmonic(3*p+2),  Proposed_Harmonic(3*p+2),  Existing_Harmonic(3*p+2); 
3*p+4,            Test_Harmonic(3*p+4),  Proposed_Harmonic(3*p+4),  Existing_Harmonic(3*p+4); 
];


Error_Table = ...
[Comparison_Table(:,1),(Comparison_Table(:,3)./Comparison_Table(:,2)-1),(Comparison_Table(:,4)./Comparison_Table(:,2)-1)];

figure(46);
bar((1:1:12),Error_Table(2:end,2:3))

        % Convert x-axis values to ordinals
   a=[cellstr(num2str(Error_Table(2:end,1)))]; 
   new_xticks = [char(a)];
% 'Reflect the changes on the plot
   set(gca,'xticklabel',new_xticks);


   grid on

   xlabel('Harmonic Order','FontSize',32);
ylabel('Percentage Error','FontSize',32); 

lgd = legend('Proposed model','Frequency-invariant model','FontSize',24,'Location', 'southwest');
lgd.Orientation = 'horizontal';
lgd.NumColumns = 1;
legend('boxoff');

set(gcf,'Position',[300 300 800 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

        % Convert y-axis values to percentage values by multiplication
   a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
   pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
   new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
   set(gca,'yticklabel',new_yticks);


f = gcf;
filename = ['HarmonicError_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

Harmonic_rms = (Ia_THD*Ia_Fundamental+Ib_THD*Ib_Fundamental+Ic_THD*Ic_Fundamental)/3
Error_rms = [sqrt(sum((Proposed_Harmonic(2:end) - Test_Harmonic(2:end)).^2))./Harmonic_rms,...
             sqrt(sum((Existing_Harmonic(2:end) - Test_Harmonic(2:end)).^2))./Harmonic_rms]

Harmonic_rms_2 = sqrt(sum(Comparison_Table(2:end,2).^2))
Error_rms_2 = [sqrt(sum((Proposed_Harmonic(Comparison_Table(2:end,1)) - Test_Harmonic(Comparison_Table(2:end,1))).^2))./Harmonic_rms_2,...
             sqrt(sum((Existing_Harmonic(Comparison_Table(2:end,1)) - Test_Harmonic(Comparison_Table(2:end,1))).^2))./Harmonic_rms_2]

tEnd = toc(tStart)

%% test voltage clamping

%{
% voltage fft
x1 = Voltage(1:N*MM);
x2 = x1;
for index = 1:length(x1)
    if x1(index) > 0
       x2(index) = Vdc/2;
    else 
       x2(index) = -Vdc/2;
    end
end

% voltage fft
%x = Voltage(1:(N*MM));
% x = Phase_Voltage_2.signals(1).values(T_start:(T_end-1));
%x = LLVoltage_1.signals(1).values(T_start:(T_end-1));
X=fft(x2);
C=X(1:1:floor((N*MM)/2))/(N*MM);
FrequencyT2 = (0:length(C)-1)'*(ff/MM); 

MM_Sample = ((1+MM):MM:floor(N*MM/2));
FrequencyT = FrequencyT2(MM_Sample);


Voltage_Mag_A_Clamped = abs(C(MM_Sample))*2;   
Voltage_Angle_A_Clamped = angle(C(MM_Sample))*(180/pi);

Voltage_DC_Clamped = abs(C(1))*2;
Van_Fundamental_Clamped=Voltage_Mag_A_Clamped(1)/sqrt(2)
Van_Harmonic_Clamped=sqrt(sum(Voltage_Mag_A_Clamped(1:length(Voltage_Mag_A_Clamped)-1).^2))/sqrt(2)
Van_THD_Clamped = sqrt(Van_Harmonic_Clamped^2-Van_Fundamental_Clamped^2)/Van_Fundamental_Clamped

figure(53);
plot(Time,x2,'-');
hold;
plot(Time,x1,'-');
legend('clamped','original','Location', 'eastoutside');
legend('boxoff');
xlabel('Time(s)');
ylabel('Voltage(V)'); 
title('Line-line AB Voltage Waveforms - Vab');
grid on;

figure(54);
% bar(order,Hac_pu)
Test_clamped_Vll_Harmonic = [Voltage_Mag_A(1:order_min)/sqrt(2), Voltage_Mag_A_Clamped(1:order_min)/sqrt(2)];
bar(Voltage_order,Test_clamped_Vll_Harmonic)
%set(gca,'yscale','log')
set(gca,'yscale','linear')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (Vrms)','FontSize',30); 
legend('Original','Clamped','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 80]);
ylim([0 2]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

figure(55);
% bar(order,Hac_pu)
Test_clamped_Vll_Harmonic = [(Hac_a(1:order_min) + Hac_b(1:order_min) +Hac_c(1:order_min))*Vnph/3, Voltage_Mag_A_Clamped(1:order_min)/sqrt(2)];
bar(Voltage_order,Test_clamped_Vll_Harmonic)
%set(gca,'yscale','log')
set(gca,'yscale','linear')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (Vrms)','FontSize',30); 
legend('Original','Clamped','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 80]);
ylim([0 2]);
set(gcf,'Position',[300 300 1600 600]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
%}

%% sequence plots
%{
%voltage

Sequence_voltage_positive = [vsa1,vsb1, vsc1]./(Vdc/sqrt(3));
Sequence_voltage_negative = [vsa2,vsb2, vsc2]./(Vdc/sqrt(3));

figure(47);
% bar(order,Hac_pu)
bar(Voltage_order,mean(abs(Sequence_voltage_positive(1:order_min,:)),2),0.5)
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (per Unit)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 1e0]);
set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['V1_VoltageSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(48);
% bar(order,Hac_pu)
bar(Voltage_order,mean(abs(Sequence_voltage_negative(1:order_min,:)),2),0.5,"FaceColor","#D95319")
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (per Unit)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 1e0]);
set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['V2_VoltageSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(49);
% bar(order,Hac_pu)
bar(Voltage_order,[mean(abs(Sequence_voltage_positive(1:order_min,:)),2),mean(abs(Sequence_voltage_negative(1:order_min,:)),2)],2)
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',40);
ylabel('Harmonic Magnitude (per Unit)','FontSize',40); 
legend('Positive Sequence','Negative Sequence','FontSize',40);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 1e0]);
set(gcf,'Position',[300 300 1600 600]);
%set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',36);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['V1V2_VoltageSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

%current

figure(50);
% bar(order,Hac_pu)
bar(Voltage_order,abs(it1(1:order_min,:)),0.5)
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (A)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 5e0]);
set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['I1_CurrentSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

figure(51);
% bar(order,Hac_pu)
bar(Voltage_order,abs(it2(1:order_min,:)),0.5,"FaceColor","#D95319")
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (A)','FontSize',30); 
%legend('Prediction','Measurement','FontSize',30);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 5e0]);
set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['I2_CurrentSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 

Phase_current = [ita_abs,itb_abs,itc_abs];
figure(52);
% bar(order,Hac_pu)
bar(Voltage_order,abs([Phase_current(1:order_min,:)]),1)
set(gca,'yscale','log')
xlabel('Harmonic Order','FontSize',30);
ylabel('Harmonic Magnitude (A)','FontSize',30); 
legend('I_a','I_b','I_c','FontSize',20);
% title('Machine Phase Current Spectrum');
% xlim([0 Ht].*fo);
xlim([0 order_shown]);
ylim([1e-3 5e0]);
set(gcf,'Position',[200 200 1050 800]);
set(gca,'FontSize',24);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
f = gcf;
filename = ['Iabc_CurrentSpectra_' num2str(fo) 'Hz_' num2str(fc) 'Hz_' num2str(Vdc) 'V_' num2str(Test_Harmonic(1),'%4.2f') 'A.eps'];
exportgraphics(f,filename,'Resolution',600);
%print(f,filename,'-dpng','-r600'); 
%}