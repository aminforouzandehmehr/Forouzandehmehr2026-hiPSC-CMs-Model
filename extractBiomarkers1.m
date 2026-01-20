function [featureVectorAP, featureVectorContr] = extractBiomarkers1(window, stimFlag, showPlots)
% extractBiomarkersNEW  Biomarkers for hiPSC-CMs (spontaneous or paced), raw (no smoothing).
% All biomarkers are returned as the MEAN OF THE LAST 5 BEATS.
% The figure (if requested) annotates ONLY the last beat.
%
% AP method for all APDxx:
%   - beat windowing (see below),
%   - MDP-referenced percent levels,
%   - dV/dt-gated ascending crossing,
%   - linear interpolation for crossings.
%
% CaT method:
%   RISE metrics (identical to CaT plot intent):
%     - Pre-peak baseline = minimum FROM (stimulus or window start) TO peak
%       (i.e., stim→peak in paced mode; windowStart→peak in spontaneous mode)
%     - FIRST crossings after that minimum within [iMinPre : iPk]:
%         t10_from_min, t90_from_min
%     - RT10topeak = t(peak) − t10_from_min
%     - RT1090     = t90_from_min − t10_from_min
%   DECAY/WIDTH metrics:
%     - Decay thresholds use POST-diastolic baseline (Cmin_post).
%     - DT9010 = 90% (post) → 10% (post), with the 10% fall using the LAST crossing.
%     - CTDxx ("time above x%") = rise at x% (pre) → fall at x% (post). (Not plotted here.)
%     - DURATION (first CaT output) = CTD10. (Kept for compatibility; computed below.)
%
% Outputs:
%   featureVectorAP   = [MDP, dV_dt_max, APA, Peak, APD10, APD20, APD30, APD50, APD70, APD80, APD90, Rate_AP, RAPP_APD, CL]
%   featureVectorContr= [peakTension, cellShortPerc, relaxTime50]

if nargin < 3, showPlots = false; end

% ---------------- Config ----------------
SCALE_TO_mV = true;   % set false if window.V_ode is already in mV
PRECL_FRAC  = 0.40;   % fallback fraction of CL before peak for the 1st beat's window (spontaneous only)
lastNBeats  = 5;      % average over last 5 beats

% ---------- Pull raw inputs ----------
t     = window.Time(:);                     % ms
Vm    = window.V_ode(:);
if SCALE_TO_mV, Vm = 1000*Vm; end           % V -> mV
Cai   = window.Cai(:);                      % mM
Frc   = window.Frc(:);                      % mN/mm^2
Istim = [];  if isfield(window,'Istim'), Istim = window.Istim(:); end
Lsrc  = [];  if isfield(window,'Lsrc'),  Lsrc  = window.Lsrc(:);  end

dt = median(diff(t));
assert(all(isfinite(dt)) && dt>0,'Time vector must be increasing with finite positive step.');

% ---------- Beat segmentation (paced vs spontaneous) ----------
pacedMode = (stimFlag==1) && ~isempty(Istim);

if pacedMode
    % Stimulus edges (rising of |Istim| > thr)
    thr = 0.2*max(abs(Istim)); if isempty(thr) || thr==0, thr = eps; end
    stimLogic = abs(Istim) > thr;
    stimStarts = find(diff([false; stimLogic])==1);

    % Enforce minimal spacing (avoid multiple samples per stim)
    minBeatMs = 250; minPkDist = max(1, round(minBeatMs/dt));
    keep = true(numel(stimStarts),1);
    if numel(stimStarts)>1, keep(2:end) = diff(stimStarts) >= minPkDist; end
    stimStarts = stimStarts(keep);

    if numel(stimStarts) < (lastNBeats+1)
        featureVectorAP    = NaN(1,13);
        featureVectorCaT   = NaN(1,11);
        featureVectorContr = NaN(1,3);
        if showPlots
            figure('Name','Beat detection'); plot(t,Vm,'k'); title('Insufficient paced beats'); grid on
        end
        return;
    end

    % Define windows stim(k) → stim(k+1)-1
    win0 = stimStarts(1:end-1);
    win1 = stimStarts(2:end)-1;

    % For AP peak reference, pick Vm peak within each window
    apIdx = zeros(numel(win0),1);
    for k = 1:numel(win0)
        [~,loc] = max(Vm(win0(k):win1(k)));
        apIdx(k) = win0(k) + loc - 1;
    end

    % Representative CL in samples
    CLsamp = round(median(diff(stimStarts), 'omitnan'));
    if ~isfinite(CLsamp) || CLsamp<=0, CLsamp = round(1000/dt); end  % ~1 Hz fallback

else
    % Spontaneous
    apIdx = detect_beats_Vm(Vm, dt);
    nIdx = numel(apIdx);
    if nIdx < (lastNBeats+1)
        featureVectorAP    = NaN(1,13);
        featureVectorCaT   = NaN(1,11);
        featureVectorContr = NaN(1,3);
        if showPlots
            figure('Name','Beat detection'); plot(t,Vm,'k'); title('Insufficient beats'); grid on
        end
        return;
    end
    CLsamp = round(median(diff(apIdx), 'omitnan'));
    if ~isfinite(CLsamp) || CLsamp<=0, CLsamp = round(300/dt); end
end

% ---------- AP metrics per beat ----------
levelsPerc = [10 20 30 40 50 70 80 90]; % fixed order (APD90 is width at 10% of APA above MDP)
idx10 = 1;

if pacedMode
    nB = numel(win0);
else
    nB = numel(apIdx)-1;
end

MDP_b   = NaN(nB,1); APA_b   = NaN(nB,1); Peak_b  = NaN(nB,1); dVdt_b  = NaN(nB,1);
APD10_b = NaN(nB,1); APD20_b = NaN(nB,1); APD30_b = NaN(nB,1); APD40_b = NaN(nB,1);
APD50_b = NaN(nB,1); APD70_b = NaN(nB,1); APD80_b = NaN(nB,1); APD90_b = NaN(nB,1);

% Keep exact windows for plotting the last beat
i0_used = zeros(nB,1); i1_used = zeros(nB,1);
tUp_last = []; tDn_last = []; Vlevels_last = []; ap_last_dbg = [];

for k = 1:nB
    if pacedMode
        i0 = win0(k); i1 = win1(k);
    else
        % previous-peak → next-peak; first beat uses pre-CL fallback
        if k > 1, i0 = apIdx(k-1) + 1; else, i0 = max(1, apIdx(k) - round(PRECL_FRAC*CLsamp)); end
        i1 = apIdx(k+1) - 1;
    end
    i0_used(k) = i0; i1_used(k) = i1;

    ap = ap_metrics_all_levels_raw(t, Vm, i0, i1, levelsPerc); % unified method

    % Scalars
    MDP_b(k)  = ap.MDP; APA_b(k)  = ap.APA; Peak_b(k) = ap.Vpeak; dVdt_b(k) = ap.dVdtmax;

    % APDxx (MDP-ref, dV/dt-gated, linear interpolation)
    APD10_b(k) = ap.dur( ap.map(10) );
    APD20_b(k) = ap.dur( ap.map(20) );
    APD30_b(k) = ap.dur( ap.map(30) );
    APD40_b(k) = ap.dur( ap.map(40) );
    APD50_b(k) = ap.dur( ap.map(50) );
    APD70_b(k) = ap.dur( ap.map(70) );
    APD80_b(k) = ap.dur( ap.map(20) );
    APD90_b(k) = ap.dur( ap.map(10) );  % APD90 = duration at 10% level

    if k == nB
        tUp_last     = ap.tUp; tDn_last = ap.tDn; Vlevels_last = ap.levelVals; ap_last_dbg  = ap;
    end
end

% Beat-to-beat CL & rate (last-N)
if pacedMode
    CLs = diff(t(stimStarts));          % ms, stim-to-stim
else
    CLs = diff(t(apIdx));               % ms, peak-to-peak
end
lastN   = (nB-lastNBeats+1):nB;
Rate_AP = 60/(mean(CLs(lastN),'omitnan')/1000);

% ---------- Average LAST 5 BEATS ----------
featureVectorAP = [ ...
    mean(MDP_b(lastN),'omitnan'), ...
    mean(dVdt_b(lastN),'omitnan'), ...
    mean(APA_b(lastN),'omitnan'), ...
    mean(Peak_b(lastN),'omitnan'), ...
    mean(APD10_b(lastN),'omitnan'), ...
    mean(APD20_b(lastN),'omitnan'), ...
    mean(APD30_b(lastN),'omitnan'), ...
    mean(APD50_b(lastN),'omitnan'), ...
    mean(APD70_b(lastN),'omitnan'), ...
    mean(APD80_b(lastN),'omitnan'), ...
    mean(APD90_b(lastN),'omitnan'), ...
    Rate_AP, ...
    mean((APD30_b(lastN)-APD40_b(lastN))./(APD70_b(lastN)-APD80_b(lastN)),'omitnan'), ...
    mean(CLs(lastN),'omitnan') ...
];

% ---------- Contractility (beatwise, then last N) ----------
if pacedMode
    % Build synthetic indices so the function windows are stim->stim
    apIdx_for_contr = [win0; win1(end)+1];
else
    apIdx_for_contr = apIdx;
end
contrBeat = contractility_metrics_raw_beats(t, Frc, Lsrc, apIdx_for_contr);
featureVectorContr = [ ...
    mean(contrBeat.peakTension(lastN),'omitnan'), ...
    mean(contrBeat.cellShortPerc(lastN),'omitnan'), ...
    mean(contrBeat.relaxTime50(lastN),'omitnan') ...
];

% ---------- OPTIONAL PLOTS (annotate ONLY the last beat) ----------
if showPlots
    kLast = nB; 
    i0p = i0_used(kLast); i1p = i1_used(kLast);
    segT  = t(i0p:i1p); segVm = Vm(i0p:i1p); segCa = Cai(i0p:i1p);

    % Last-beat APD90 items (cached)
    V10_last   = Vlevels_last(idx10);
    tUp10_last = tUp_last(idx10);
    tDn10_last = tDn_last(idx10);
    APD90_last = tDn10_last - tUp10_last;

    % Plot 1: APD90 (Vm + dV/dt), last beat ONLY
    draw_apd90_plot(segT, segVm, ap_last_dbg.preT, ap_last_dbg.preV, ap_last_dbg.preDV, ...
        ap_last_dbg.iPkLocal, ap_last_dbg.MDP, ap_last_dbg.Vpeak, V10_last, ...
        ap_last_dbg.t_takeoff, tUp10_last, tDn10_last, APD90_last, ap_last_dbg.thrFrac);

    % Plot 3: MDP (same window), last beat ONLY
    %draw_mdp_plot_full(segT, segVm, ap_last_dbg.iPkLocal, ap_last_dbg.MDP);
end


end % ===== end main =====


% ===================== Helper functions (no nesting) =====================

% ================= CaT Helpers =================
function out = empty_out()
out = struct('mean_RT10topeak',NaN,'mean_RT1090',NaN,'mean_DT9010',NaN, ...
             'mean_MaxCa',NaN,'mean_MinCa',NaN);
end

function apIdx = detect_beats_Ca(Ca, dt)
% Simple & safe Ca peak detector
Ca = Ca(:); N = numel(Ca);
rngC = max(Ca)-min(Ca);
if ~isfinite(rngC) || rngC < 1e-9, apIdx = []; return; end
promFrac  = 0.08;                   % 8% of range
minBeatMs = 250;                    % lower bound for CaTs
prom      = max(1e-6, promFrac*rngC);
minPkDist = round(0.6*minBeatMs/dt);
minPkDist = max(1, min(minPkDist, max(1,N-1)));  % keep scalar, < signal span
try
    [~, apIdx] = findpeaks(Ca, 'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);
catch
    [~, apIdx] = findpeaks(Ca);
end
apIdx = apIdx(:);
% Need at least 2 peaks for 1 complete beat
if numel(apIdx) < 2
    apIdx = [];
end
end

function tcr = crossing_time2(tx, yx, level, dir)
% First crossing with linear interpolation; safe for 1-point input
tx = tx(:); yx = yx(:);
tcr = NaN; if isempty(tx), return; end
if strcmp(dir,'asc'), logic = yx >= level; else, logic = yx <= level; end
k = find(diff([false; logic])==1, 1, 'first');
if isempty(k), return; end
k0 = max(1,k-1); k1 = k;
if k1==k0 || numel(tx)<2
    % fallback to nearest sample if interpolation impossible
    tcr = tx(k1);
    return;
end
den = (yx(k1)-yx(k0)); if abs(den)<eps, den = sign(den)*eps + (den==0)*eps; end
frac = (level - yx(k0))/den;
tcr  = tx(k0) + frac*(tx(k1)-tx(k0));
end

function tcr = crossing_time_last(tx, yx, level)
% Last descending crossing with linear interpolation; safe for 1-point input
tx = tx(:); yx = yx(:);
tcr = NaN; if isempty(tx), return; end
logic = yx <= level;
k = find(diff([false; logic])==1, 1, 'last');
if isempty(k), return; end
k0 = max(1,k-1); k1 = k;
if k1==k0 || numel(tx)<2
    tcr = tx(k1);
    return;
end
den = (yx(k1)-yx(k0)); if abs(den)<eps, den = sign(den)*eps + (den==0)*eps; end
frac = (level - yx(k0))/den;
tcr  = tx(k0) + frac*(tx(k1)-tx(k0));
end

% ================= AP Helpers =================

function apIdx = detect_beats_Vm(Vm, dt)
rngV = max(Vm)-min(Vm);
if rngV < 1e-3, apIdx = []; return; end
prom = max(1, 0.05*rngV);
minBeatMs = 250;
minPkDist = max(1, round(0.6*minBeatMs/dt));
[~, apIdx] = findpeaks(Vm, 'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);
end

% ----- Unified AP metric (ALL levels) -----
function ap = ap_metrics_all_levels_raw(t, Vm, i0, i1, levelsPerc)
segT = t(i0:i1);
segV = Vm(i0:i1);

% Peak & baseline (MDP pre-peak)
[~, iPkLocal] = max(segV);
Vpeak = segV(iPkLocal);
if iPkLocal>1, MDP = min(segV(1:iPkLocal)); else, MDP = min(segV); end
APA  = Vpeak - MDP;

% dV/dt and takeoff (pre-peak only; raw)
dV = diff(segV)./diff(segT); 
if isempty(dV)
    dV = 0;
else
    dV(end+1) = dV(end);
end
preT  = segT(1:iPkLocal);
preV  = segV(1:iPkLocal);
preDV = dV(1:iPkLocal); preDV(preDV<0) = 0;
dVdtmax = max(preDV); if ~isfinite(dVdtmax) || dVdtmax<=0, dVdtmax = max(dV); end
thrFrac = 0.20;
takeoff_k = find(preDV >= thrFrac*dVdtmax, 1, 'first');
if isempty(takeoff_k)
    takeoff_k = find(preDV >= 0.1*dVdtmax, 1, 'first');
    if isempty(takeoff_k)
        dt_ms = median(diff(segT));
        if ~isfinite(dt_ms) || dt_ms<=0, dt_ms = 1; end
        takeoff_k = max(1, iPkLocal - round(10/dt_ms));
    end
end
t_takeoff = preT(takeoff_k);

% Levels (absolute values)
levelsPerc = levelsPerc(:).';
levelVals  = MDP + (levelsPerc/100).*APA;

% Post-peak segment
postT = segT(iPkLocal:end);
postV = segV(iPkLocal:end);

% Crossings (safe interpolation)
tUp = NaN(size(levelVals)); tDn = NaN(size(levelVals));
for m = 1:numel(levelVals)
    preT2 = preT(takeoff_k:end); preV2 = preV(takeoff_k:end);
    tUp(m) = crossing_time(preT2, preV2, levelVals(m), 'asc');
    tDn(m) = crossing_time(postT, postV, levelVals(m), 'desc');
end
dur = tDn - tUp;

% Mapping percent → index
ap.map = containers.Map('KeyType','double','ValueType','double');
for ii=1:numel(levelsPerc), ap.map(levelsPerc(ii)) = ii; end

% Pack
ap.MDP      = MDP; ap.APA = APA; ap.Vpeak = Vpeak; ap.dVdtmax = dVdtmax;
ap.levelVals= levelVals; ap.tUp = tUp; ap.tDn = tDn; ap.dur = dur;
ap.preT = preT; ap.preV = preV; ap.preDV = preDV; ap.iPkLocal = iPkLocal;
ap.t_takeoff = t_takeoff; ap.thrFrac = thrFrac;
end

% --------- Crossing helpers (first-cross) ---------
function tcr = crossing_time(tx, yx, level, dir)
% first crossing (ascending or descending) with interpolation; safe for short vectors
tx = tx(:); yx = yx(:);
tcr = NaN;
if isempty(tx), return; end
if strcmp(dir,'asc'), logic = yx >= level; else, logic = yx <= level; end
k = find(diff([false; logic])==1, 1, 'first');
if isempty(k), return; end
k0 = max(1,k-1); k1 = k;
if k1==k0 || numel(tx)<2
    tcr = tx(k1);
    return;
end
den = (yx(k1)-yx(k0)); if abs(den) < eps, den = sign(den)*eps + (den==0)*eps; end
frac = (level - yx(k0))/den;
tcr  = tx(k0) + frac*(tx(k1)-tx(k0));
end


% ----------------- Contractility metrics (beatwise) -----------------
function out = contractility_metrics_raw_beats(t, Frc, Lsrc, apIdx_like)
nB = numel(apIdx_like)-1;
peakT  = NaN(nB,1); RT50 = NaN(nB,1); shortP = NaN(nB,1);

for k=1:nB
    i0 = apIdx_like(k); i1 = apIdx_like(k+1)-1;
    segT = t(i0:i1); segF = Frc(i0:i1);

    if isempty(segT), continue; end

    [Fp, iPk] = max(segF); peakT(k)  = Fp;
    F50 = 0.5*Fp;
    if iPk < numel(segF)
        t50 = crossing_time(segT(iPk:end), segF(iPk:end), F50, 'desc');
        if ~isnan(t50), RT50(k) = t50 - segT(iPk); end
    end

    if ~isempty(Lsrc)
        segL = Lsrc(i0:i1);
        if ~isempty(segL)
            L0   = max(segL);   % diastolic approx
            Lmin = min(segL);   % systolic approx
            if L0>0, shortP(k) = (L0 - Lmin)/L0 * 100; end
        end
    end
end

out = struct('peakTension',peakT,'relaxTime50',RT50,'cellShortPerc',shortP);
end

% ===================== Plotting helpers (raw) =====================

function draw_apd90_plot(segT, segVm, preT, preV, preDV, iPkLocal, ...
        MDP, Vpeak, V10, t_takeoff, tUp10, tDn10, APD90_last, thrFrac)

cBand = [0 0.45 0.74]; cMDP = [0 .5 0]; cPk = [.5 0 .5];

figure('Name','APD90 (last beat, raw)','Color','w');
tStart = segT(1); tEnd = segT(end);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- Top: Vm with APD90 band (last beat ONLY) ---
nexttile; hold on
patch([tStart tEnd tEnd tStart], [min(segVm) min(segVm) max(segVm) max(segVm)], ...
    [0.96 0.96 0.96], 'EdgeColor','none','DisplayName','Beat window');

plot(segT, segVm, 'k-', 'LineWidth',1.3, 'DisplayName','V_m (raw)');

yline(MDP, '--', sprintf(' MDP = %.1f mV', MDP), 'Color', cMDP, 'LabelHorizontalAlignment','left');
yline(Vpeak,'--', sprintf(' Peak = %.1f mV', Vpeak), 'Color', cPk,  'LabelHorizontalAlignment','left');
yline(V10,  ':',  ' 10% level', 'Color', cBand, 'LabelHorizontalAlignment','left');

if all(isfinite([tUp10 tDn10])) && tDn10>tUp10
    patch([tUp10 tDn10 tDn10 tUp10], [V10 V10 MDP MDP], ...
        cBand, 'FaceAlpha',0.12, 'EdgeColor',cBand, 'LineWidth',1.2, 'DisplayName','APD90 band');
    plot([tUp10 tDn10], [V10 V10], '-', 'Color', cBand, 'LineWidth',2);
    text(mean([tUp10 tDn10]), V10, sprintf(' APD90 = %.2f ms', APD90_last), ...
        'VerticalAlignment','bottom','Color',cBand,'FontWeight','bold');
end

% Robust Vm value at t_takeoff (avoid interp1 errors with <2 points)
vm_takeoff = local_interp_safe(segT, segVm, t_takeoff);

plot(t_takeoff, vm_takeoff, 'x', ...
    'Color',[0 0.4 0.8],'LineWidth',1.5,'MarkerSize',9,'DisplayName','Takeoff (dV/dt gated)');

xlabel('Time (ms)'); ylabel('Vm (mV)');
title('APD90 from MDP-referenced 10% level (raw)');
legend('Location','best'); grid on; box on; xlim([tStart tEnd]);

% --- Bottom: dV/dt with threshold and takeoff ---
nexttile; hold on
plot(preT, preDV, '-', 'Color',[0.2 0.2 0.2], 'LineWidth',1.2, 'DisplayName','dV/dt (pre-peak, raw)');

% Threshold line
if ~isempty(preDV)
    dvmax = max(preDV);
else
    dvmax = 0;
end
thr   = thrFrac*dvmax;
yline(thr, ':', sprintf(' threshold = %.2f (%.0f%% of max)', thr, 100*thrFrac), ...
    'Color',[0.25 0.25 0.7], 'LabelHorizontalAlignment','left');

% Robust dV/dt at t_takeoff
dv_takeoff = local_interp_safe(preT, preDV, t_takeoff);

if isfinite(dv_takeoff)
    plot(t_takeoff, dv_takeoff, 'o', ...
        'Color',[0.25 0.25 0.7], 'MarkerFaceColor','w', 'LineWidth',1.2, 'DisplayName','Takeoff');
end

xline(segT(iPkLocal), '--', 'Peak Vm', 'Color',[.4 .4 .4]);
xlabel('Time (ms)'); ylabel('dV/dt (mV/ms)');
title('Upstroke gating by dV/dt (raw)');
legend('Location','best'); grid on; box on; xlim([tStart tEnd]);

end

function vq = local_interp_safe(x, y, xq)
% 1-point/2-point safe sampler for plotting annotations
x = x(:); y = y(:);
if isempty(x)
    vq = NaN; return
elseif numel(x)==1
    vq = y(1); return
else
    % linear interpolation with extrap allowed
    vq = interp1(x, y, xq, 'linear', 'extrap');
end
end

function draw_mdp_plot_full(segT, segVm, iPkLocal, MDP)
% Full-window MDP plot matched to the APD90 plot window (raw)
cMDP  = [0 .5 0]; tStart = segT(1); tEnd = segT(end);
[~, iMinPre] = min(segVm(1:iPkLocal));  % MDP point (pre-peak minimum)
tMDP = segT(iMinPre); Vmdp = segVm(iMinPre);

figure('Name','MDP (last beat, raw)','Color','w'); hold on
patch([tStart segT(iPkLocal) segT(iPkLocal) tStart], ...
      [min(segVm) min(segVm) max(segVm) max(segVm)], ...
      [0.96 0.96 0.96], 'EdgeColor','none', 'DisplayName','pre-peak window');

plot(segT, segVm, 'k-', 'LineWidth',1.3, 'DisplayName','V_m (raw)');
plot(tMDP, Vmdp, 'v', 'Color', cMDP, 'MarkerFaceColor', cMDP, 'MarkerSize',7, 'DisplayName','MDP');

txt = sprintf(' MDP = %.1f mV\n t = %.1f ms', MDP, tMDP);
text(tMDP, Vmdp, txt, 'VerticalAlignment','top','Color',cMDP,'FontWeight','bold');

yline(MDP,'--','Color',cMDP,'DisplayName','MDP level');
xline(segT(iPkLocal),'--','Color',[.5 .5 .5],'Label',' Peak','LabelVerticalAlignment','bottom');
xlabel('Time (ms)'); ylabel('Vm (mV)');
title('Maximum Diastolic Potential (raw) in the same window as APD90 plot');
legend('Location','best'); grid on; box on; xlim([tStart tEnd]);

end
