
%===============================================================================
% Version JMCC PLUS (DOI: TBD)
% Author: Amin Forouzandehmehr
%===============================================================================

%%
clear; close all; clc
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);
tic
%% Myofilament original initial state
Nxb = 0.9997;  % Non-permissive fraction in overlap region
XBpreR = 0.0001;   % Fraction of XBs in pre-rotated state (Stiffness generation)
XBpostR = 0.0001;  % Fraction of XBS in post-rotated state (force generation)
x_XBpreR = 0;    % Strain of XBs in the pre-rotated state
x_XBpostR = 0.007;  % Strain of XBs in the post-rotated state
TropCaL = 0.0144; % Ca bound to the low affinity troponin site
TropCaH = 0.1276; % Ca bound to the high affinity troponin sit
IntegF = 0;
SL = 1.9; % Sarcomere length

init = [Nxb XBpreR XBpostR x_XBpreR x_XBpostR TropCaL TropCaH IntegF SL];

% Y=[-0.070  0.32    0.0002  0    0    1     1     1      0      1      0   0.75  0.75  0   0.1    1    0    9.2    0     0.75    0.3     0.9     0.1    0.97    0.01    0.01    0.01    1.9     0.0070  0   0   0   0   0   0];
% Yn = Y(1:23);
% Yn(24:32) = init;
% Yn(33) = [0.133];

% Loading steady state initial conditions
load 1.5Hz_ss_control_yfinal.mat
%load 1Hz_ss_control_yfinal.mat
%load spnt_ss_control_yfinal.mat
Yn(33) = [0.1772];    
Yn(34) = [0.5]; Yn(35) = [0.5]; 

%% For sensitivity analsyses or GA optimizations otherwise 1
GAC = ones(1,25);

%% 0: Temperature = 37°C , 1: Temperature = 21°C
classicVSoptic = 0;

%% 0: spontaneous beating, 1: paced (currently 1.5 Hz; to change go to Stimulation section of hiMCES2025)
stimFlag = 1;
% Set the twitch protocol
contrflag = 1;     % 0 = isometric  1 = active contraction  2 = No contraction

%% ischemia setting
% For hiMCES_SET version (arrhythmic ischemia-reperfusion simulations) ischemiaFlag must be 1.

ischemiaFlag = 0; % 0 = no ischemia, 1 = SEV1, 2 = SEV2 
myFlags = zeros(1,7); % 1: Hyperkalemia OFF:0 / ON:1, 2: Acidosis OFF:0 / ON:1, 3: Hypoxia IKatp OFF:0 / ON:1, 4: Hypoxia on INaK OFF:0 / ON:1, 5: Hypoxia on IpCa OFF:0 / ON:1, 6: Hypoxia on Jup OFF:0 / ON:1, 7: Hypoxia on CE OFF:0 / ON:1.
myFlags = [0  0  0  0  0  0  0];
IKatpSelectionString = 'Kazbanov'; % it should take a string: Ledezma or Kazbanov

syms x
% z takes the extent of inhibition of INaK & IpCa and calculates
% corresponding normalized O2s
% fnp takes the level of oxygantion (0 to 1) and calculates the corresponding
% inhibition for INaK & IpCa

z = 1; fnp = 0;                   % Physoxia: z = 1; fnp = 0;    

if ischemiaFlag == 1
z = solve (1-rho(x) == 0.2, x);
z = double(z);
fnp = 1 - rho(z);

else if ischemiaFlag == 2
z = solve (1-rho(x) == 0.31, x);
z = double(z);
fnp = 1 - rho(z);

end
end

INaFRedMed=1; INaLRedMed=1; INaKRedMed=1; ICaLRedMed=1; IKrRedMed=1; IKsRedMed=1; 
INaCaRedMed=1; IfRedMed=1; SERCARedMed=1; SERCAKmRedMed=1; RyRRedMed=1; RyRPoRedMed=1; ItoRedMed = 1; fkon = 1; flv = 1;
fmyo = ones (1,15); tDrugApplication = 1e5; BlebFlag = 0;
%% Integration
% @hiMCES is the standard model
gma = 22; % std = 22;
[t,Yc] = ode15s(@hiMCES2025,[0 400],Yn, options, contrflag, GAC, stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, BlebFlag);

Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
Nai  = Yc(:,18);
mpi_ef = Yc(:,34);

                                               
%% Isoproterenol (10 µM) block.
tDrugApplication = 177770;
BlebFlag = 0;
% Turn on at tISO; isoLevel ~ 1 at 10 µM using a Hill activation.
%tISO    = 10;                     % s
dose_nM = 1e4;  EC50_nM = 30;  nH = 1.2;
isoLevel = 1 * (dose_nM^nH) / (EC50_nM^nH + dose_nM^nH);
isoLevel = 1;               % saturated effect

% --- Tunable gains (typical starting values) ---
s.ICaL = 0.3;           % +50% ICaL conductance    0.5
s.RyR  = 0.3;           % +50% release gain (RyR)  0.5
s.IKs  = 0.5;           % +50% IKs conductance    0.5
s.IKr  = 0.3;           % +25% IKr conductance
s.SERCAVmax = 0.3;        % +50% SERCA Vmax      0.5
s.SERCaKm   = 0.3;      % -50% SERCA Km (higher affinity) --> you can decrease kd,cai in the code
s.RyRPo = 0.3;          % -50% RyR Po half activation constant (higher affinity)
s.INaK = 0.1;          % +50% Na/K pump (via phospholemman/PKA) Def: 0.10; 0 corresponds to fpca_controlvsISO_annotated_v2.fig
s.If = 0.1;            % +25% If

s.NCX = 0.0;
s.NaF = 0.0;
s.Ito = 0.0;

% Myofilament (TnI/MyBP-C): ↓Ca sensitivity, faster cycling/lusitropy (fpca_controlvsISO_annotated_v3.fig)
s.kon   = 0.1;  s.koffL = 0.1;  s.koffH = 0.1;  s.perm50 = 0.1;
s.fapp  = 0.3;  s.gxb   = 0.0;   s.hb    = 0.0;   s.titin  = 0.1;

% --- Apply ISO scaling (create *_iso) ---
% ICaL, RyR, IKs
ICaLRedMed  = 1 * (1 + s.ICaL*isoLevel);
RyRRedMed = 1 * (1 + s.RyR *isoLevel);
IKsRedMed   = 1  * (1 + s.IKs *isoLevel);
IKrRedMed = 1  * (1 + s.IKr *isoLevel);
ItoRedMed = 1  * (1 - s.Ito *isoLevel);
IfRedMed = 1  * (1 + s.If *isoLevel);
INaCaRedMed = 1  * (1 - s.NCX *isoLevel);
INaFRedMed = 1  * (1 - s.NaF *isoLevel);

% SERCA (choose param or direct flux scaling depending on your code)
SERCARedMed = 1 * (1 + s.SERCAVmax*isoLevel);
SERCAKmRedMed = 1 * (1 - s.SERCaKm*isoLevel);
RyRPoRedMed =  1 * (1 - s.RyRPo*isoLevel);
% Km_up_iso   = 1   * (1 - s.SERCaKm  *isoLevel);
% i_up_iso    = 1    * (1 + s.SERCAVmax*isoLevel);   % fallback if you don't expose Vmax/Km

% INaK (affects both membrane current and ATP use)
INaKRedMed = (1 - s.INaK*isoLevel);

% Myofilament params
%fmyo = ones (1,9);
kon_iso    = 1    * (1 - s.kon*isoLevel);
koffL_iso  = 1  * (1 + s.koffL*isoLevel);
koffH_iso  = 1  * (1 + s.koffH*isoLevel);
perm50_iso = 1 * (1 + s.perm50*isoLevel);
fapp_iso   = 1   * (1 + s.fapp*isoLevel);
gxb_iso    = 1    * (1 + s.gxb*isoLevel);
hb_iso     = 1     * (1 + s.hb *isoLevel);
PExp_titin_iso = 1 * (1 - s.titin*isoLevel);
PCon_titin_iso = 1 * (1 - s.titin*isoLevel);
fmyo(1:9) = [kon_iso, koffL_iso, koffH_iso, perm50_iso, fapp_iso, gxb_iso, hb_iso, PExp_titin_iso, PCon_titin_iso];

if tDrugApplication > 18000
    fmyo = ones(1,15);
end
fpca2(fmyo)

%% Blebbistatin (2.5 µM) block
% BlebFlag = 1;
% tDrugApplication = 10;
% fmyo = ones (1,15);
% INaFRedMed=1; INaLRedMed=1; INaKRedMed=1; ICaLRedMed=1; IKrRedMed=1; IKsRedMed=1; INaCaRedMed=1; IfRedMed=1; SERCARedMed=1; SERCAKmRedMed =1;  RyRRedMed=1; RyRPoRedMed =1;
% 
% % To simulate 2.5 uM Blebbistatin adapted from Forouzandehmehr et al., 2022 (DOI: 10.3389/fphys.2022.1010786)
% F1 = 5.015;
% F2 = 0.1;
% TropReg = 1.0;
% ap2 = 0.03;          % ap2 = 0.03035714285714286;
% ap3 = 1.0;                         % ap3 = 0.3;
% am2 = 1.0;                         % am2 = 0.25; 
% perm500 = 1.4;                     % Def: perm50 = 1.4;
% 
% fmyo(10:15) = [F1, F2, TropReg, ap2, ap3, am2];
% fmyo(4) = perm500;
% %fmyo = ones (1,15);
% fpca2(fmyo)

%% Run
Yn = Yc (end,:); Yn(33) = [0.1437];  %Yn(1) = -0.0752; % for clamp

[t,Yc] = ode15s(@hiMCES2025,[0 400],Yn, options, contrflag, GAC, stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, BlebFlag);

Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
Nai  = Yc(:,18);
mpi_ef = Yc(:,34);

for i= 1:size(Yc,1)
[~, dati]    = hiMCES2025(t(i), Yc(i,:), contrflag, GAC, stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, BlebFlag);
    INa(i)   = dati(1);
    If(i)    = dati(2);
    ICaL(i)   = dati(3);
    Ito(i)   = dati(4);
    IKs(i)   = dati(5);
    IKr(i)   = dati(6);
    IK1(i)   = dati(7);
    INaCa(i) = dati(8);
    INaK(i)  = dati(9);
    IpCa(i)  = dati(10);
    IbNa(i)  = dati(11);
    IbCa(i)  = dati(12);
    Irel(i)  = dati(13);
    Iup(i)   = dati(14);
    Ileak(i) = dati(15); 
    Istim(i) = dati(16);
    E_K(i)   = dati(17);
    E_Na(i)  = dati(18);
    INaL(i)  = dati(19);
    Fps(i)  = dati(20);
    AT(i)  = dati(21);
    ATPase(i)  = dati(22);
    Lsarc(i) = dati(23);
    Velo(i) = dati(24);
    JCB(i) = dati(25);
    PSP(i) = dati(26);          % SERCA phosphorylation
    IKATP(i) = dati(27);
    ronak(i) = dati(28);         
    O2s(i) = dati(29);          % Bath/source O2 concentration   
    O2e(i) = dati(30);          % Extracellular O2 concentration
    vcy(i) = dati(31);          % SERCA cycle rate
    do2dt(i) = dati(32);
    fhib(i) = dati(33);         % Returns finhib
    Koo(i) = dati(34);
    ro1(i) = dati(35);
    ro2(i) = dati(36);
    fkatp(i) = dati(37);
    kon(i) = dati(38);
    IPiezo1(i) = dati(39);
    JPiezo1(i) = dati(40);
    tau_Piezo1(i) = dati(41);
    Ep_Piezo1(i) = dati(42);
    minf_Piezo1(i) = dati(43);
    OCR(i) = dati(44);
    OCR_XB(i) = dati(45);
    OCR_nonXB(i) = dati(46);
    PRF(i) = dati(47);
    
    pmol_min_cell(i) = dati(48);
    pmol_min_well(i) = dati(49);

    OCR_NaK(i) = dati(50);
    OCR_leak(i) = dati(51);
    OCR_SERCA(i) = dati(52);
    OCR_PMCA(i) = dati(53);
    
end
result       = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; E_K; E_Na; INaL];
mat_currents = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; INaL];
I_tot=sum(mat_currents);

toc
%uEnd

%% Plots

xlimit = [0 100]; fs = 14; lw = 1.2;          % xlimit = [780 820];
figure
subplot(3,2,1)
plot(t,1000*Vm), ylabel('Em (mV)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs), 

subplot(3,2,3)
plot(t, 1e6*Cai), ylabel('Cai (nM)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs), 

subplot(3,2,5)
plot(t, AT), ylabel('Tension (kPa)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs), 

subplot(3,2,2)
plot(t, O2e), ylabel('[O_2]_e (mM)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs)

subplot(3,2,6)
plot(t, OCR), ylabel('OCR (mM/s)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),xlabel('Time (s)')

subplot(3,2,4)
conv_uM_per_pct = 9.58;                 % µM per % O2 (your value)
mM_to_pct       = @(mM) (1000/conv_uM_per_pct) .* mM;   % mM -> %O2
mMs_to_pct_s    = @(mMs) (1000/conv_uM_per_pct) .* mMs; % mM/s -> %O2/s
% ---- wherever you log/plot ----
O2_peri_pct        = mM_to_pct(O2e);      % pericellular %O2
O2_bulk_pct        = mM_to_pct(O2s);        % (optional) bulk/source %O2
dO2_peri_pct_per_s = mMs_to_pct_s(OCR);% rate in %O2/s
dO2_peri_pct_per_m = 60 * dO2_peri_pct_per_s; % %O2/min (optional)
%OBL_pct = mM_to_pct(0.166);
plot(t,O2_peri_pct),ylabel('Pericellular O2 (%)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),xlabel('Time (s)')
set(gca,'box','off','tickdir','out','fontsize',fs), 
lines = findobj(gcf, 'Type', 'line');
set(lines, 'LineWidth', lw);


set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 30, 18]);  % [x y width height]
%exportgraphics(gcf, 'Spont_NewPiezo_F.png', 'Resolution', 800);

figure
subplot(3,2,1)
plot(t,1000*Vm), ylabel('Em (mV)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(3,2,2)
plot(t,1e6*Cai), ylabel('Cai (nM)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(3,2,3)
plot(t,AT), ylabel('Tension (kPa)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(3,2,4)
plot(t,Lsarc), ylabel('Sarcomere length (\mum)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(3,2,5)
plot(t,abs(OCR)), ylabel('OCR (mM/s)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(3,2,6)
plot(t,caSR), ylabel('CaSR (mM)'), xlim([390 400]), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

lines = findobj(gcf, 'Type', 'line');
set(lines, 'LineWidth', lw);

%% Biomarker calculation
stimFlag

Frc = AT; Lsrc = Lsarc; svl = Velo;
window = extractWindowForBiomarkers(Vm,I_tot, Istim ,t,Cai,Nai,Frc,Lsrc,svl, 30);

stimFlag  = stimFlag;                     % spontaneous
showPlots = true;                  % optional
[featureVectorAP,  featureVectorContr] = extractBiomarkers1(window, stimFlag, showPlots); % Retruns action potential and contractile biomarkers

MDP = featureVectorAP(1)
UV = featureVectorAP(2)
APA = featureVectorAP(3)
APPeak = featureVectorAP(4)
APD10 =  featureVectorAP(5);
APD20 =  featureVectorAP(6);
APD30 = featureVectorAP(7);
APD50 = featureVectorAP(8);
APD70 = featureVectorAP(9);
APD80 = featureVectorAP(10)
APD90 = featureVectorAP(11)
Rate_AP = featureVectorAP(12)

peakTension = featureVectorContr(1)
cellShortPerc = featureVectorContr(2)
relaxTime50 = featureVectorContr(3)


biomk = extractBiomarkers2(window)         % Retruns CaT biomarkers

Hz = Rate_AP/60;
CL = 1/Hz;

OCRcntr = OCR_XB;
OCRelec = OCR_nonXB;
% Define the time window for one cycle or it can be the last N seconds
jdx = t >= 93 & t <= 93+CL;    % for 0.75Hz: 1.3333 for 1.5Hz: 0.6667 | For spnt: idx = t >= 0.12 & t <= 1.82; For ISO (OLD): t >= 239-0.8333 & t <= 239;
t_segment = t(jdx);
OCR_segment = OCR(jdx);
PRF_segment = PRF(jdx);
OCRcntr_segment = OCRcntr(jdx);
OCRelec_segment = OCRelec(jdx);
OCR_NaK_segment = OCR_NaK(jdx);
OCR_leak_segment = OCR_leak(jdx);
OCR_SERCA_segment = OCR_SERCA(jdx);
OCR_PMCA_segment = OCR_PMCA(jdx);

IpCa_segment = IpCa(jdx);
INaK_segment = INaK(jdx);
AT_segment = AT(jdx);
pmol_min_cell_seg = pmol_min_cell(jdx);


% Integrate OCR over that segment
intOCR = cumtrapz(t_segment, OCR_segment);        % result in mM
intPRF = cumtrapz(t_segment, PRF_segment);        % result in mM
intOCRC = cumtrapz(t_segment, OCRcntr_segment); 
intOCRE = cumtrapz(t_segment, OCRelec_segment); 
intOCRNaK = cumtrapz(t_segment, OCR_NaK_segment);
intOCRleak = cumtrapz(t_segment, OCR_leak_segment);
intOCRSERCA = cumtrapz(t_segment, OCR_SERCA_segment);
intOCRPMCA = cumtrapz(t_segment, OCR_PMCA_segment);

intPMCA = cumtrapz(t_segment, IpCa_segment);
intINaK = cumtrapz(t_segment, INaK_segment);
intAT = cumtrapz(t_segment, AT_segment);

total_OCR = abs(intOCR (end))
total_contractile_OCR = abs(intOCRC (end));
total_elec_OCR = abs(intOCRE (end));
total_OCR_NaK = abs(intOCRNaK (end));
total_OCR_leak = abs(intOCRleak (end));
total_OCR_SERCA = abs(intOCRSERCA (end));
total_OCR_PMCA = abs(intOCRPMCA (end));

totalPMCA = abs(intPMCA (end));
totalINaK = abs(intINaK (end));

totalAT = abs(intAT (end))

OCR_NaK_percent = 100 * (total_OCR_NaK/total_OCR)
OCR_leak_percent = 100 * (total_OCR_leak/total_OCR)
OCR_myo_percent = 100 * (total_contractile_OCR/total_OCR)
OCR_SERCA_percent = 100 * (total_OCR_SERCA/total_OCR)
OCR_PMCA_percent = 100 * (total_OCR_PMCA/total_OCR)

