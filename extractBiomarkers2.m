function biomk = extractBiomarkers2(window)

% Computes Ca-transient biomarkers
%   minCa, maxCa, CaT_Duration (CTD10), CTD80,
%   RT10topeak, RT1050, RT1090, DT9010, Rate_Cai (Hz)
%
% Inputs (required fields in 'window'):
%   window.Time : time in ms
%   window.Cai  : Ca concentration (mM)
%
% Output (struct):
%   biomk.minCa
%   biomk.maxCa
%   biomk.CaT_Duration    % CTD10 = time above 10% (pre-rise → last 10% fall)
%   biomk.CTD80           % time above 80% (pre-rise → first 80% fall)
%   biomk.RT10topeak      % t(peak) - first 10% rise
%   biomk.RT1050          % first 50% rise - first 10% rise
%   biomk.RT1090          % first 90% rise - first 10% rise
%   biomk.DT9010          % first 90% fall (post) → last 10% fall (post)
%   biomk.Rate_Cai        % Hz, = 1000 / CL_ms


% ---------- Defaults ----------
biomk = struct('minCa',NaN,'maxCa',NaN,'CaT_Duration',NaN,'CTD80',NaN, ...
               'RT10topeak',NaN,'RT1050',NaN,'RT1090',NaN,'DT9010',NaN,'Rate_Cai',NaN);

% ---------- Inputs ----------
t  = window.Time(:);
Ca = window.Cai(:);
if isempty(t) || isempty(Ca) || numel(t) ~= numel(Ca), return; end

% ---------- NaN handling ----------
finiteMask = isfinite(Ca);
if ~all(finiteMask)
    if mean(~finiteMask) > 0.25, return; end           % too corrupt
    Ca = fillmissing(Ca,'linear','EndValues','nearest');
end

% ---------- Detect troughs (minima) ----------
rngC = max(Ca) - min(Ca);
if ~isfinite(rngC) || rngC < 1e-9, return; end

dt   = median(diff(t)); if ~isfinite(dt) || dt<=0, dt = 1; end
minPkDist = max(1, round(0.6*250/dt));  % >= ~150 ms in samples
promMin   = max(1e-6, 0.05*rngC);       % minima prominence (5% of range)

try
    [~, minIdx] = findpeaks(-Ca, 'MinPeakProminence', promMin, 'MinPeakDistance', minPkDist);
catch
    [~, minIdx] = findpeaks(-Ca);
end
minIdx = minIdx(:);
if numel(minIdx) < 2, return; end

% Last complete beat: trough -> trough
iMin0 = minIdx(end-1);
iMin1 = minIdx(end);
segT  = t(iMin0:iMin1);
segC  = Ca(iMin0:iMin1);

% Cycle length & rate (Hz)
CL_ms = segT(end) - segT(1);
if CL_ms > 0, biomk.Rate_Cai = 1000 / CL_ms; end

% Must contain a peak between troughs
[ Cpk, iPkLocal ] = max(segC);
if ~(isfinite(Cpk) && iPkLocal>=1 && iPkLocal<=numel(segC)), return; end

biomk.maxCa = Cpk;
biomk.minCa = min(segC);

% Baselines: use actual troughs at both ends
Cmin_pre  = segC(1);      % start trough
Cmin_post = segC(end);    % end trough

amp_pre  = Cpk - Cmin_pre;
amp_post = Cpk - Cmin_post;
if ~(amp_pre > 0 && amp_post > 0), return; end

% Segments
preT  = segT(1:iPkLocal);
preC  = segC(1:iPkLocal);
postT = segT(iPkLocal:end);
postC = segC(iPkLocal:end);
tPk   = segT(iPkLocal);

% Levels
lev10_pre  = Cmin_pre  + 0.10*amp_pre;
lev50_pre  = Cmin_pre  + 0.50*amp_pre;
lev80_pre  = Cmin_pre  + 0.80*amp_pre;
lev90_pre  = Cmin_pre  + 0.90*amp_pre;

lev10_post = Cmin_post + 0.10*amp_post;
lev80_post = Cmin_post + 0.80*amp_post;
lev90_post = Cmin_post + 0.90*amp_post;

% ---- Rise crossings ----
t10r = first_cross(preT, preC, lev10_pre, 'asc');
if isnan(t10r) && ~isempty(preC) && preC(1) >= lev10_pre, t10r = preT(1); end

t50r = first_cross(preT, preC, lev50_pre, 'asc');
if isnan(t50r) && ~isempty(preC) && preC(1) >= lev50_pre, t50r = preT(1); end

t80r = first_cross(preT, preC, lev80_pre, 'asc');
if isnan(t80r) && ~isempty(preC) && preC(1) >= lev80_pre, t80r = preT(1); end

t90r = first_cross(preT, preC, lev90_pre, 'asc');
if isnan(t90r) && ~isempty(preC) && preC(1) >= lev90_pre, t90r = preT(1); end

% If rise never reaches a level, keep NaN
if max(preC) < lev10_pre, t10r = NaN; end
if max(preC) < lev50_pre, t50r = NaN; end
if max(preC) < lev80_pre, t80r = NaN; end
if max(preC) < lev90_pre, t90r = NaN; end

% ---- Fall crossings (clamp missing to end to avoid NaN widths) ----
t10f = last_cross (postT, postC, lev10_post);
if isnan(t10f) && ~isempty(postT), t10f = postT(end); end

t80f = first_cross(postT, postC, lev80_post, 'desc');
if isnan(t80f) && ~isempty(postT), t80f = postT(end); end

t90f = first_cross(postT, postC, lev90_post, 'desc');
if isnan(t90f) && ~isempty(postT), t90f = postT(end); end

% ---- Durations / RTs / DT ----
% DURATION (CTD10)
if isfinite(t10r) && isfinite(t10f) && t10f > t10r
    biomk.CaT_Duration = t10f - t10r;
end

% CTD80
if isfinite(t80r) && isfinite(t80f) && t80f > t80r
    biomk.CTD80 = t80f - t80r;
end

% RT10topeak, RT1050, RT1090 (only if rise crossings exist)
if isfinite(t10r)
    biomk.RT10topeak = tPk - t10r;
    if isfinite(t50r), biomk.RT1050 = t50r - t10r; end
    if isfinite(t90r), biomk.RT1090 = t90r - t10r; end
end

% DT9010 = first 90% fall → last 10% fall (post-anchored)
if isfinite(t90f) && isfinite(t10f) && t10f > t90f
    biomk.DT9010 = t10f - t90f;
end

end % ===== main =====


% ===== helpers =====
function tcr = first_cross(tx, yx, level, dir)
% First crossing with linear interpolation; dir = 'asc' or 'desc'
tx = tx(:); yx = yx(:); tcr = NaN; 
if numel(tx) < 2, return; end
switch lower(dir)
    case 'asc',  logic = (yx >= level);
    case 'desc', logic = (yx <= level);
    otherwise, error('dir must be ''asc'' or ''desc''.');
end
k = find(diff([false; logic])==1, 1, 'first');
if isempty(k), return; end
k0 = max(1,k-1); k1 = k;
dy = (yx(k1)-yx(k0)); if abs(dy) < eps, dy = sign(dy)*eps + (dy==0)*eps; end
tcr = tx(k0) + (level - yx(k0))/dy * (tx(k1)-tx(k0));
end

function tcr = last_cross(tx, yx, level)
% Last descending crossing with linear interpolation
tx = tx(:); yx = yx(:); tcr = NaN; 
if numel(tx) < 2, return; end
logic = (yx <= level);
k = find(diff([false; logic])==1, 1, 'last');
if isempty(k), return; end
k0 = max(1,k-1); k1 = k;
dy = (yx(k1)-yx(k0)); if abs(dy) < eps, dy = sign(dy)*eps + (dy==0)*eps; end
tcr = tx(k0) + (level - yx(k0))/dy * (tx(k1)-tx(k0));
end
