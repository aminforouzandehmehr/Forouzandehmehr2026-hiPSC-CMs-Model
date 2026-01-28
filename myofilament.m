function [dY, Tension] = myofilament(time, Y, tDrugApplication, Cai_um, fmyo)

% Model parameters

 contrflag = 0;

% Metabolite reference concentrations
Pi_ref = 2; %mM             Def: 2
MgATP_ref = 5; %mM          Def: 5
MgADP_ref = 0.036; %mM      Def: 0.036

% MgADP dissociation constant
kdADP = 0.004; % uM

% Set metabolite concentrations and pH level
MgATP = 5; %mM              Def: 5
MgADP = 36e-3; %mM          Def: 36e-3
Pii = 2; %mM                Def: 2
pH = 7.15; %                Def: 7.15    baraye levo calib: 6.7.             

% Competitive H binding parameters
kdHCa = 1e-5;
m = 1;
pH_ref = 7.15;              %Def: 7.15
Href = 10^(-pH_ref)*1e3;

% Gas constant
R = 8.314;
TmpC = 37;

% Sarcomere geometry (um)
SL_max = 2.4;
SL_min = 1.4;
L_thick = 1.65;
L_hbare = 0.1;
L_thin = 1.2;

% Temperature dependence
%TmpC = 25;
Q_kon = 1.5;
Q_koff = 1.3;
Q_kn_p = 1.6;
Q_kp_n = 1.6;
Q_fapp = 6.25;
Q_gapp = 2.5;
Q_hf = 6.25;
Q_hb = 6.25;
Q_gxb = 6.25;

% Species parameter
Xsp = 0.2; 

% XB density
rho = 0.25;                 % rho = 0.25;

% No of ATP consumed per XB 
phi = 1;

% Ca binding to troponin and thin filament
% fmyo = [kon_iso, koffL_iso, koffH_iso, perm50_iso, fapp_iso, gxb_iso, hb_iso, PExp_titin_iso, PCon_titin_iso];
kon = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(1))*1.25*50e3;   % (mM-1 s-1) fkon for 0.3 uM Levo: 1.1, 1 uM Levo: 1.21, 2 uM Levo: 1.32, 3 uM Levo: 1.39, 5 uM Levo: 1.57, 10 uM: 1.88, 15 uM: 2.2, 20 uM: 2.51.
koffL = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(2)) * 200; % (s-1) Def: 250e-3
koffH = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(3)) * 25;  % (s-1)
perm50 = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(4)) * 0.6; 
n_perm = 0.77*15;
kn_p = 550; % (s-1) % 
kp_n = 50; % (s-1)  % 
koffmod = 0.5;        % mod to change species

% Thin filament regulation and XB cycling
fapp = 1.0*500 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(5)); % (s-1)
gapp = 1.0*70;  % (s-1)
gslmod = 1.0*6;
hf = 1*2000;  % (s-1)
hfmdc = 5;
hb = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(7)) * 1.0*400;   % (s-1)
gxb = ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(6)) * 1.0*70;   % (s-1)
sigma_p = 8;
sigma_n = 1;   

% Mean strain of strongly-bond states
x0 = 0.007; % (um)
Psi = 2;

% Normalised active and passive force
SL_rest = 1.9;  % (um)
PCon_titin = 1.0*0.002; %(normalized force) Def: 0.002
PExp_titin = 1.0*10; %(um-1)
SL_collagen = 2.25; %(uM)
PCon_collagen = 1*0.02; %(normalized force)
PExp_collagen  = 1*70; %(um-1)

% Calculation of complete muscle response
mass = 4e-1*5e-5; % for Rat (normalised force s^2 um-1)
viscosity = 0.003;  % (normalised force s um-1)
    
    % Assign the state variables
    Nxb = Y(24);  % Non-permissive fraction in overlap region
    XBpreR = Y(25);   % Fraction of XBs in pre-rotated state (Stiffness generation)
    XBpostR = Y(26);  % Fraction of XBS in post-rotated state (force generation)
    x_XBpreR = Y(27);    % Strain of XBs in the pre-rotated state
    x_XBpostR = Y(28);  % Strain of XBs in the post-rotated state
    TropCaL = Y(29); % Ca bound to the low affinity troponin site
    TropCaH = Y(30); % Ca bound to the high affinity troponin site
    IntegF = Y(31); % Integral of force
    SL = Y(32); % Sarcomere length
 
    Afterload = 0; % Afterload for work-loop contraction. Def: 0
    loop = 0; % Boolean to indicate if it is a work-loop contraction
    preload_SL = 1.9; % Initial SL for work-loop contractions
    T_loop(1) = 2000; % Time for start of relaxation phasea (end of isotonic)
    T_loop(2) = 2000; % Tme for end of relaxation phase 
    passive = 1;      % Boolean for presence of passive force
    kxb = 13.1047;    % Def: 13.1047;

    % Prevent over extension
    if SL > preload_SL
        dSL = 0;
    end

  % *****************************************************    
    % Ca binding to troponin
    
    % Call function to get Ca
    Cai = Cai_um; %uM
    
    % convert pH to mM
    H = 10^(-pH)*1e3; % H concentration in mM
    
    % Ca binding to troponin C
    konT = kon*1*Q_kon^((TmpC-37)/10);
    konT_app = konT*(kdHCa^m + Href^m)/(kdHCa^m + H^m);
    
    koffLT = koffL*koffmod*1*Q_koff^((TmpC-37)/10); 
    koffHT = koffH*koffmod*1*Q_koff^((TmpC-37)/10);
    
    d_TropCaL = konT_app*Cai*(1-TropCaL) - koffLT*TropCaL;
    d_TropCaH = konT_app*Cai*(1-TropCaH) - koffHT*TropCaH;
    
    
    % *****************************************************
    % Thin filament activation rates
    
    % Sarcomere geometry
    sovr_ze = min(L_thick/2, SL/2);
    sovr_cle = max(SL/2 - (SL-L_thin),L_hbare/2);
    L_sovr = sovr_ze - sovr_cle; % Length of single overlap region
    
    % Overlap fraction for thick filament
    sov_thick = L_sovr*2/(L_thick - L_hbare);
    % Overlap fraction for thin filament
    sov_thin = L_sovr/L_thin;
    
    TropReg = ((1-sov_thin)*TropCaL + sov_thin*TropCaH) * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(12));
    permtot = (1/(1+(perm50/TropReg)^n_perm))^0.5;
    permtot_p_n = min(100,1/permtot);
    
    % Rate constants governing the transition btw Permissive and Non
    kn_pT = kn_p*permtot*Q_kn_p^((TmpC-37)/10);
    kp_nT = kp_n*permtot_p_n*Q_kp_n^((TmpC-37)/10);

    
    % *****************************************************
    % Cross-bridge cycling rates
    
    % Pxb to XBpreR
    fappT = fapp*Xsp*Q_fapp^((TmpC-37)/10);

    % Pi-dependent transition XBpreR to Pxb
    gappslmd = 1 + (1-sov_thick)*gslmod;
    gappT = gapp*gappslmd*Xsp*Q_gapp^((TmpC-37)/10);
    gappT_true = gappT/Pi_ref;  % True first order rate constant
    
    % XBpreR to XBpostR
    hfmd = exp(-sign(x_XBpreR)*hfmdc*(x_XBpreR/x0)^2);
    hfT = hf*hfmd*Xsp*Q_hf^((TmpC-37)/10);
    
    % H-dependent transition XBpostR to XBpreR
    hbT = hb*Xsp*Q_hb^((TmpC-37)/10);
    hbT_true = hbT/Href;
    hbT_app = ((kdADP+MgADP_ref)/MgADP_ref)*(MgADP/(kdADP+MgADP))*hbT_true;  

    % MgATP-dependent transition from XBpostR to Pxb
    if (x_XBpostR < x0)
        gxbmd = exp(sigma_p*((x0-x_XBpostR)/x0)^2);
    else
        gxbmd = exp(sigma_n*((x0-x_XBpostR)/x0)^2);
    end
    
    gxbT = gxb*max(gxbmd,1)*Xsp*Q_gxb^((TmpC-37)/10);
    gxbT_true = gxbT/MgATP_ref;
    gxbT_app = ((kdADP+MgADP_ref)/(kdADP+MgADP))*gxbT_true;

    % Pxb to XBpostR - Introduced for thermodynamic efficiency
    G0 = -29600 + log(10)*R*(TmpC+273)*(-log10(1e-7));
    K_MgATP = exp(-G0./(R.*(TmpC+273)))*1e6;  % 1e6 to convert from M2 to mM2
    fxbT = (kdADP*fappT*hfT*gxbT_true)/(gappT_true*hbT_true*K_MgATP);
    fxb = (kdADP*fapp*hf*(gxb/MgATP_ref))/((gapp/Pi_ref)*(hb/Href)*K_MgATP); % Used for calculating max
    
    ap1 = fappT;
    ap2 = hfT * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(13));
    ap3 = MgATP*gxbT_app * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(14));
    am1 = Pii*gappT_true;
    am2 = H*hbT_app * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(15));
    am3 = fxbT;
    F1 = 1 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(10));
    F2 = 1 * ((time<tDrugApplication)*1+(time >= tDrugApplication)*fmyo(11));
    
    % *****************************************************
    %  Cross_bridge transitions
    
    % Update all RUs   
    Pxb = 1 - Nxb - XBpreR - XBpostR;
    d_Nxb = -kn_pT*Nxb + kp_nT*Pxb;
    
    d_XBpreR = ap1*F1*Pxb - am1*F2*XBpreR - ap2*XBpreR + am2*XBpostR;
    d_XBpostR = ap2*XBpreR + am3*Pxb - am2*XBpostR - ap3*XBpostR;
    
    ATPase = rho*phi*sov_thick*(ap3*XBpostR - am3*Pxb);
    
    % *****************************************************
    % Mean strain of strongly-bound states
           
    % Steady State occupancies of the bound states - Duty fractions
    Sum_duty = am3*am2 + ap3*ap1 + am2*ap1 + ap1*ap2 + am3*am1 + am3*ap2...
        + ap2*ap3 + am3*am1 + ap3*am1;

    dutyPreR = (am3*am2 + ap3*ap1 + am2*ap1)/Sum_duty;
    dutyPostR = (ap1*ap2 + am3*am1 + am3*ap2)/Sum_duty;
    
    % No shortening during isotonic period
    % This is to hold the SL at max contraction to get a loop
    if (time>T_loop(1)) & (time<T_loop(2)) & loop
        dSL = 0;
    end

    % Compute the active force
    F_active = sov_thick*(x_XBpreR*XBpreR + x_XBpostR*XBpostR);

    % Maximal state occupancies under optimal conditions
    max_XBpreR = (fapp*(hb+gxb)+fxb*fapp)/(gxb*hf + fapp*hf + gapp*hb + ...
        gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));

    max_XBpostR = (fapp*hf+fxb*(gapp+hb))/(gxb*hf + fapp*hf + gapp*hb + ...
        gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));

    % Factor to normalise force
    FnormD = x0*max_XBpostR;
    
    % Normalised Active force
    Fnorm = F_active/FnormD;
    
    % Normalised Passive force
    F_passive = passiveForces(SL,time,tDrugApplication,fmyo);
    F_after = Afterload;
    
    % No shortening during isotonic period
    % This is to hold the SL at max contraction to get a loop
    if (time>T_loop(1)) & (time<T_loop(2)) & loop
        dSL = 0;
    end
    
    % Difference in force
    d_Force =  F_after -  Fnorm - F_passive;
    
    % Total normalised force
    F_total = Fnorm + abs(F_passive);
    
    % Used in Fevents script 
    F_net = F_total - F_after;
    
    % SL dynamics
    
    dSL = contrflag*((IntegF + viscosity*(SL_rest - SL))/mass);
    %F_passive = passiveForces(SL,passive);
    %dIntegF = (-F_passive-Fnorm+Afterload);
    % Rate of change of the average distortions
    d_x_XBpreR = dSL/2 + (Psi/dutyPreR)*(-ap1*x_XBpreR) + am2*(x_XBpostR-x0-x_XBpreR);
    d_x_XBpostR = dSL/2 + (Psi/dutyPostR)*(ap2*(x_XBpreR + x0 - x_XBpostR));
    
    % Ca buffering by low-affinity troponin C (LTRPNCa)
    Trop_conc = 70e-3;       % (mM) troponin concentration
    
    FrSBXB    = (XBpostR+XBpreR)/(max_XBpostR + max_XBpreR);         %Done
    dFrSBXB   = (d_XBpostR +d_x_XBpreR)/(max_XBpostR + max_XBpreR);  %Done

    dsovr_ze  = -dSL/2*heav(L_thick-SL);
    dsovr_cle = -dSL/2*heav((2*L_thin-SL)-L_hbare);
    dlen_sovr = dsovr_ze-dsovr_cle;
    dSOVFThin = dlen_sovr/L_thin;

    dTropTot= Trop_conc*(-dSOVFThin*TropCaL+(1-sov_thin)*d_TropCaL + ...        
    dSOVFThin*(FrSBXB*TropCaL+(1-FrSBXB)*TropCaL) + ...
    sov_thin*(dFrSBXB*TropCaL+FrSBXB*d_TropCaH-dFrSBXB*TropCaL+...
    (1-FrSBXB)*d_TropCaL));
     
    % Assembling the derivative vector
    dY(24:32,1) = [d_Nxb, d_XBpreR, d_XBpostR, d_x_XBpreR, d_x_XBpostR, d_TropCaL, d_TropCaH, d_Force, dSL];
    
    JltrpnLH = dTropTot;
    JCaBMyo = JltrpnLH;
    frc = kxb*Fnorm;           % Only reporting active
    frc = kxb*F_total;         % Reporting active + passive
    Lsarc = SL;
    cvelo   = dSL;

    Tension = frc;
end