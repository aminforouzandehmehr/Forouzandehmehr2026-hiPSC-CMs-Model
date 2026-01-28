function fpca(fmyo)
% Simulate Tension-pCa relationship
% Integrates each pCa to steady state, reads Tension, fits/annotates a Hill curve.
% clear; close all; clc
%% ---------- user settings ----------
pCa_vec      = 9.0:-0.1:4.5;     % conventional: left->right (axis will be reversed)
T_end        = 300;               % seconds per point
tDrugApplication = 0;          % no drug switch during sweep    inf
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

% Initial state (length >= 32; your function uses indices 24..32)
Y0 = zeros(32,1);
Y0(24) = 0.70;  % Nxb
Y0(25) = 0.15;  % XBpreR
Y0(26) = 0.15;  % XBpostR
Y0(27) = 0.00;  % x_XBpreR
Y0(28) = 0.00;  % x_XBpostR
Y0(29) = 0.00;  % TropCaL
Y0(30) = 0.00;  % TropCaH
Y0(31) = 0.00;  % IntegF
Y0(32) = 1.90;  % SL (µm) isometric at SL_rest

%% ---------- sweep pCa ----------
F_ss = zeros(size(pCa_vec));   % steady-state Tension from your function

for k = 1:numel(pCa_vec)
    pCa    = pCa_vec(k);
    Cai_uM = 10.^(-pCa) * 1e3;   % pCa -> [Ca] in µM

    % ODE RHS: integrate using your dY
    rhs = @(t,y) myofilament_rhs_wrapper(t,y,tDrugApplication,Cai_uM,fmyo);

    % integrate to steady state
    [~,Y] = ode15s(rhs, [0 T_end], Y0, opts);
    Yss = Y(end,:).';

    % read out tension at steady state
    [~, Tension] = myofilament(T_end, Yss, tDrugApplication, Cai_uM, fmyo);

    F_ss(k) = Tension;
    Y0 = Yss;   % warm start
end

% normalize and fit Hill (4-parameter)
%Fplot = F_ss / 9.34677;
Fplot = F_ss / max(F_ss);

hill4 = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^( p(3) .* (x - p(4)) ));
p0    = [min(Fplot), max(Fplot), 3.0, 5.6];
obj   = @(p) sum((Fplot - hill4(p,pCa_vec)).^2);
opts2 = optimset('Display','off','TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',5e4,'MaxIter',5e4);
pFit  = fminsearch(obj, p0, opts2);

Fhat = hill4(pFit, pCa_vec);
SSE  = sum((Fplot - Fhat).^2);
SST  = sum((Fplot - mean(Fplot)).^2);
R2   = 1 - SSE/SST;

nH = pFit(3); pCa50 = pFit(4);
EC50_uM = 10^(-pCa50) * 1e6;

% plot
%figure('Color','w');
%openfig('pCa_control.fig');                     % not normalized
%openfig('pCa_control_normalized.fig');          % normalized
%hold on;
plot(pCa_vec, Fplot, 'o', 'MarkerFaceColor',[0.2 0.4 0.9], 'MarkerEdgeColor','k', 'LineWidth',1.0);
pFine = linspace(min(pCa_vec), max(pCa_vec), 400);
plot(pFine, hill4(pFit, pFine), '-', 'LineWidth', 2);
set(gca,'XDir','reverse'); grid on; box on;
xlabel('pCa = -log_{10}([Ca^{2+}]_i / M)'); ylabel('Tension (kPa)');
title('Tension–pCa');

txt = {
    sprintf('n_H = %.2f', nH)
    sprintf('pCa_{50} = %.3f', pCa50)
    sprintf('EC_{50} = %.3g \\muM', EC50_uM)
    sprintf('R^2 = %.3f', R2)
    };
annotation('textbox',[0.15 0.15 0.25 0.2], 'String', txt, ...
    'FitBoxToText','on', 'BackgroundColor',[0.97 0.97 0.97], ...
    'EdgeColor',[0.6 0.6 0.6], 'FontSize',11);
%legend({'Model (steady state)', 'Hill fit'}, 'Location','southwest');

end

% --------- helper: wraps the myofilament for ODE15s ----------
function dY = myofilament_rhs_wrapper(t, Y, tDrugApplication, Cai_um, fmyo)
    [dY, ~] = myofilament(t, Y, tDrugApplication, Cai_um, fmyo);
end
