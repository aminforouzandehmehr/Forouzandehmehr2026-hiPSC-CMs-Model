%===============================================================================
% Version JMCC PLUS (DOI: TBD)
% Author: Amin Forouzandehmehr
%
% PURPOSE:
%   Integrates hiMCES2026 model with optional ischemia, Isoproterenol and Blebbistatin.
%
% MAIN FLAGS:
%   stimFlag        : 0 spontaneous, 1 paced
%   freq            : 0.6 - 1.7 Hz
%   contrflag       : 0 isometric, 1 active contraction, 2 no contraction
%   classicVSoptic  : 0 = 37C, 1 = 21C
%   ischemiaFlag    : 0 none, 1 SEV1, 2 SEV2
%   ISOFlag         : 0 OFF, 1 ON
%   BlebFlag        : 0 OFF, 1 ON
%
% INITIAL CONDITIONS:
%   - Whole-cell state loaded from MAT file (Yn)
%===============================================================================

%% Housekeeping
clear; close all; clc
tic
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);

%% ============================================================================
%  1) FLAGS 
% ============================================================================
stimFlag        = 1;          % 0 spontaneous, 1 paced
freq            = 1.0;        % Frequency (Hz) only works if stimFlag = 1
contrflag       = 1;          % 0 isometric, 1 active contraction, 2 no contraction
classicVSoptic  = 0;          % Temperature: 0: 37C, 1: 21C

ischemiaFlag    = 0;        % 0 none, 1 SEV1, 2 SEV2
myFlags         = [0 0 0 0 0 0 0];   % [HyperK, Acidosis, IKatp, INaK, IpCa, Jup, CE]
IKatpSelectionString = 'Kazbanov';   % 'Ledezma' or 'Kazbanov'

% Drug flags: Isoproterenol and Blebbistatin
ISOFlag         = 0;        % 0 OFF, 1 ON 
BlebFlag        = 0;        % 0 OFF, 1 ON 

% Drug timing (seconds)
tDrugApplication_ISO  = 0;      % used only if ISOFlag == 1
tDrugApplication_Bleb = 0;      % used only if BlebFlag == 1

% Global application time:
tDrugApplication = 1e5;          % default 

%% ============================================================================
%  2) INITIAL CONDITIONS
% ============================================================================

% Load whole-cell steady-state initial conditions

  load 1Hz_ss_control_yfinal.mat
% load spnt_ss_control_yfinal.mat
% load ISO_ss_1Hz_yfinal.mat
% load ISO_ss_spnt_yfinal.mat
% load BLEB_ss_1Hz_yfinal.mat
% load BLEB_ss_spnt_yfinal.mat

Yn = Yn(1:33);    

%% ============================================================================
%  3) PARAMETER DEFAULTS
% ============================================================================
GAC = ones(1,25);               % for sensitivity / GA; keep ones for standard

% Reduction/scaling placeholders (start at 1 = no effect)
INaFRedMed=1; INaLRedMed=1; INaKRedMed=1; ICaLRedMed=1; IKrRedMed=1; IKsRedMed=1;
INaCaRedMed=1; IfRedMed=1; SERCARedMed=1; SERCAKmRedMed=1; RyRRedMed=1; RyRPoRedMed=1;
ItoRedMed = 1;

fkon = 1;
flv  = 1;
fmyo = ones(1,15);

gma  = 22;                      % std = 22

%% ============================================================================
%  4) ISCHEMIA PRE-COMPUTE
% ============================================================================
syms x
z   = 1;
fnp = 0;                         % Physioxia: z=1, fnp=0

if ischemiaFlag == 1
    z   = double(solve(1 - rho(x) == 0.2, x));
    fnp = 1 - rho(z);
elseif ischemiaFlag == 2
    z   = double(solve(1 - rho(x) == 0.31, x));
    fnp = 1 - rho(z);
end

%% ============================================================================
%  5) DRUG SIMULATION (ISO / Bleb)
% ============================================================================
% ---- ISO (10 µM) block
if ISOFlag == 1
    tDrugApplication = tDrugApplication_ISO;

    dose_nM = 1e4;  EC50_nM = 30;  nH = 1.2;
    isoLevel = (dose_nM^nH) / (EC50_nM^nH + dose_nM^nH);
    isoLevel = 1;   % saturated effect (keep if that's what you want)

    % Tunable gains
    s.ICaL = 0.3;  s.RyR  = 0.3;  s.IKs  = 0.5;  s.IKr = 0.3;
    s.SERCAVmax = 0.3;  s.SERCaKm = 0.3;  s.RyRPo = 0.3;
    s.INaK = 0.1;  s.If   = 0.1;

    s.NCX = 0.0;  s.NaF = 0.0;  s.Ito = 0.0;

    % Myofilament (TnI/MyBP-C)
    s.kon=0.1; s.koffL=0.1; s.koffH=0.1; s.perm50=0.1;
    s.fapp=0.3; s.gxb=0.0; s.hb=0.0; s.titin=0.1;

    % Apply ISO scaling
    ICaLRedMed   = 1 * (1 + s.ICaL*isoLevel);
    RyRRedMed    = 1 * (1 + s.RyR *isoLevel);
    IKsRedMed    = 1 * (1 + s.IKs *isoLevel);
    IKrRedMed    = 1 * (1 + s.IKr *isoLevel);

    ItoRedMed    = 1 * (1 - s.Ito *isoLevel);
    IfRedMed     = 1 * (1 + s.If  *isoLevel);
    INaCaRedMed  = 1 * (1 - s.NCX *isoLevel);
    INaFRedMed   = 1 * (1 - s.NaF *isoLevel);

    SERCARedMed     = 1 * (1 + s.SERCAVmax*isoLevel);
    SERCAKmRedMed   = 1 * (1 - s.SERCaKm  *isoLevel);
    RyRPoRedMed     = 1 * (1 - s.RyRPo    *isoLevel);

    INaKRedMed      = 1 * (1 - s.INaK     *isoLevel);

    % Myofilament vector update (1:9 only)
    kon_iso          = 1 * (1 - s.kon  *isoLevel);
    koffL_iso        = 1 * (1 + s.koffL*isoLevel);
    koffH_iso        = 1 * (1 + s.koffH*isoLevel);
    perm50_iso       = 1 * (1 + s.perm50*isoLevel);
    fapp_iso         = 1 * (1 + s.fapp *isoLevel);
    gxb_iso          = 1 * (1 + s.gxb  *isoLevel);
    hb_iso           = 1 * (1 + s.hb   *isoLevel);
    PExp_titin_iso   = 1 * (1 - s.titin*isoLevel);
    PCon_titin_iso   = 1 * (1 - s.titin*isoLevel);

    fmyo(1:9) = [kon_iso, koffL_iso, koffH_iso, perm50_iso, ...
                fapp_iso, gxb_iso, hb_iso, PExp_titin_iso, PCon_titin_iso];

    %fpca(fmyo)
end

% ---- Blebbistatin (2.5 µM) block
if BlebFlag == 1
    tDrugApplication = tDrugApplication_Bleb;

    % Reset all reductions (Bleb isolated)
    INaFRedMed=1; INaLRedMed=1; INaKRedMed=1; ICaLRedMed=1; IKrRedMed=1; IKsRedMed=1;
    INaCaRedMed=1; IfRedMed=1; SERCARedMed=1; SERCAKmRedMed=1; RyRRedMed=1; RyRPoRedMed=1;
    ItoRedMed=1;

    fmyo = ones(1,15);

    % Bleb parameters (Forouzandehmehr et al., 2022; 10.3389/fphys.2022.1010786)
    F1 = 5.015;
    F2 = 0.1;
    TropReg = 1.0;
    ap2 = 0.03;
    ap3 = 1.0;
    am2 = 1.0;
    perm500 = 1.4;

    fmyo(10:15) = [F1, F2, TropReg, ap2, ap3, am2];
    fmyo(4) = perm500;

    %fpca(fmyo)
end

%% ============================================================================
%  6) INTEGRATION 
% ============================================================================

[t,Yc] = ode15s(@hiMCES2026,[0 400],Yn, options, ...
    contrflag, GAC, stimFlag, freq, classicVSoptic, tDrugApplication, ...
    INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, ...
    INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ...
    ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, ISOFlag);

Yn = Yc(end,:); Yn(33) = 0.1772;

[t,Yc] = ode15s(@hiMCES2026,[0 400],Yn, options, ...
    contrflag, GAC, stimFlag, freq, classicVSoptic, tDrugApplication, ...
    INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, ...
    INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ...
    ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, ISOFlag);

Yn = Yc(end,:); 
%% ============================================================================
%  7) POST-PROCESS 
% ============================================================================
Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
Nai  = Yc(:,18);

for i= 1:size(Yc,1)
[~, dati]    = hiMCES2026(t(i), Yc(i,:), contrflag, GAC, stimFlag, freq, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, ISOFlag);
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
    OCR(i) = dati(39);
    OCR_XB(i) = dati(40);
    OCR_nonXB(i) = dati(41);
    PRF(i) = dati(42);
    pmol_min_cell(i) = dati(43);
    pmol_min_well(i) = dati(44);
    OCR_NaK(i) = dati(45);
    OCR_leak(i) = dati(46);
    OCR_SERCA(i) = dati(47);
    OCR_PMCA(i) = dati(48);
    
end

mat_currents = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; INaL];
I_tot=sum(mat_currents);

toc

%% Plots

fs = 14; lw = 1.2; 
xlimit = [391 398];
figure
subplot(4,1,1)
plot(t,1000*Vm), ylabel('Vm (mV)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(4,1,2)
plot(t,1e6*Cai), ylabel('Cai (nM)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(4,1,3)
plot(t,AT), ylabel('Tension (kPa)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

subplot(4,1,4)
plot(t,abs(OCR)), ylabel('OCR (mM/s)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),
set(gca,'box','off','tickdir','out','fontsize',fs),

lines = findobj(gcf, 'Type', 'line');
set(lines, 'LineWidth', lw);


xlimit = [0 400];          
figure
conv_uM_per_pct = 9.58;                                 % µM per % O2 
mM_to_pct       = @(mM) (1000/conv_uM_per_pct) .* mM;   % mM -> %O2
mMs_to_pct_s    = @(mMs) (1000/conv_uM_per_pct) .* mMs; % mM/s -> %O2/s

O2_peri_pct        = mM_to_pct(O2e);                    % pericellular %O2
O2_bulk_pct        = mM_to_pct(O2s);                    % bulk/source %O2
dO2_peri_pct_per_s = mMs_to_pct_s(OCR);                 % rate in %O2/s
dO2_peri_pct_per_m = 60 * dO2_peri_pct_per_s;           % %O2/min 

plot(t,O2_peri_pct),ylabel('Pericellular O2 (%)'), xlim(xlimit), set(gca,'box','off','tickdir','out','fontsize',fs),xlabel('Time (s)')
set(gca,'box','off','tickdir','out','fontsize',fs), 
lines = findobj(gcf, 'Type', 'line');
set(lines, 'LineWidth', lw);

figure
fpca(fmyo)

%% Biomarker calculation

if stimFlag == 1
    disp(['Paced at ' num2str(freq) ' Hz']);
elseif stimFlag == 0
    disp('Spontaneous beating');
else
    disp('Error: stimFlag can only take 0 (spontaneous) or 1 (paced)');
end


Frc = AT; Lsrc = Lsarc; svl = Velo;
window = extractWindowForBiomarkers(Vm,I_tot, Istim ,t,Cai,Nai,Frc,Lsrc,svl, 30);

stimFlag  = stimFlag;                     % spontaneous
showPlots = false;                         % optional
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
Rate_AP = featureVectorAP(12);           % bpm

peakTension = featureVectorContr(1)
cellShortPerc = featureVectorContr(2)
relaxTime50 = featureVectorContr(3)


biomk = extractBiomarkers2(window)         % Retruns CaT biomarkers

Hz = Rate_AP/60;
CL = 1/Hz;

OCRcntr = OCR_XB;
OCRelec = OCR_nonXB;
% Define the time window for one cycle (jdx = t >= 393 & t <= 393+CL) or a minute
jdx = t > 340 & t <= 400;    
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
