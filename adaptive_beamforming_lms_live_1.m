%% Adaptive beamforming via LMS — LIVE, staged animated demo
%  6-element half-wavelength ULA. One SOI (known Gold-code waveform, unknown
%  direction) + two tone jammers. Pure LMS (gamma = 1) learns weights that
%  pass the SOI and null the jammers.
%
%  The animation is STAGED to tell the feedback story, one beat at a time:
%     (1) weights adjust  ->  (2) beampattern responds  ->  (3) error is
%     measured, which feeds back into the next weight update.
%
%  Press Run. No toolboxes required.

clear; clc; close all;
set(groot, 'DefaultFigureVisible', 'on');
% rng(0);                 % uncomment for a reproducible run

%% ---- Parameters ----
theta_soi_deg = 30;
sample_rate = 1e6;  d = 0.5;  N = 10000;  Nr = 6;
t = (0:N-1) / sample_rate;

theta_soi = deg2rad(theta_soi_deg);
theta2    = deg2rad(60);
theta3    = deg2rad(-50);

mu    = 1e-4;             % LMS step size (the learning rate)
gamma = 1;               % 1 = pure LMS.  <1 = leaky LMS (regularizing leak)

% ---- animation pacing ----
n_beats     = 20;        % number of staged story-beats across the run
stage_pause = 0.4;       % seconds per sub-step; total ~ n_beats*3*stage_pause
smooth_win  = 200;       % moving-average window for the learning curve

refresh = round(N / n_beats);

%% ---- Steering vectors ----
elem = (0:Nr-1).';                       % element indices 0..Nr-1
s1 = exp(-2j*pi*d*elem*sin(theta_soi));
s2 = exp(-2j*pi*d*elem*sin(theta2));
s3 = exp(-2j*pi*d*elem*sin(theta3));

%% ---- SOI: length-127 Gold code, oversampled and repeated ----
gold_code = [-1 1 1 -1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 1 1 ...
             -1 1 1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 -1 ...
             -1 -1 -1 -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 -1 -1 1 1 1 -1 1 -1 1 ...
             1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 -1 -1 -1 1 1 1 1 -1 1 1 1 -1 1 -1 ...
             -1 -1 -1 1 -1 1 1 -1 -1 -1 -1 1 -1 1 1 -1 -1 -1 -1 -1 -1 1 1];
soi = repelem(gold_code, 8);
soi = repmat(soi, 1, floor(N/numel(soi)) + 1);
soi = soi(1:N);

%% ---- Received signal ----
r = s1*soi + s2*exp(2j*pi*0.02e6*t) + s3*exp(2j*pi*0.03e6*t);
r = r + 0.5*(randn(Nr,N) + 1j*randn(Nr,N));

%% ---- Beampattern helper ----
N_fft = 1024;
theta_bins = asin(linspace(-1, 1, N_fft));
pattern_dB = @(w) 10*log10(abs(fftshift(fft([conj(w(:).'), zeros(1,N_fft-Nr)]))).^2);
idx = @(th) find(abs(theta_bins - th) == min(abs(theta_bins - th)), 1);
iS = idx(theta_soi); i2 = idx(theta2); i3 = idx(theta3);

%% ---- Figure & layout : 2 x 4, polar occupies the central 2x2 block ----
%   col1: amplitude / phase      centre 2x2: polar      col4: 3D / error
cmap = turbo(Nr);                          % element index -> color (shared)
fig = figure('Name','LMS adaptive beamforming (staged live)', ...
             'Position',[60 70 1480 820]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing','compact', 'Padding','compact');

% col1 top: amplitude
axA = nexttile(tl, 1);
hAmp = scatter(axA, elem, nan(Nr,1), 80, cmap, 'filled');
grid(axA,'on'); xlim(axA,[-0.5 Nr-0.5]); ylim(axA,[0 1.05]);
xticks(axA, 0:Nr-1); xlabel(axA,'Element index'); ylabel(axA,'|w| (relative)');
title(axA,'Weight amplitude');

% col1 bottom: phase
axP = nexttile(tl, 5);
hPh = scatter(axP, elem, nan(Nr,1), 80, cmap, 'filled');
grid(axP,'on'); xlim(axP,[-0.5 Nr-0.5]); ylim(axP,[-180 180]);
xticks(axP, 0:Nr-1); yticks(axP,-180:90:180);
xlabel(axP,'Element index'); ylabel(axP,'\angle w (deg)');
title(axP,'Weight phase');

% centre: polar beampattern (2x2 block -> large & round)
pax = polaraxes(tl); pax.Layout.Tile = 2; pax.Layout.TileSpan = [2 2];
hPol = polarplot(pax, theta_bins, nan(1,N_fft), 'LineWidth', 2); hold(pax,'on');
polarplot(pax, [theta_soi theta_soi], [-40 0], 'g--', 'LineWidth', 1.5);
polarplot(pax, [theta2 theta2],       [-40 0], 'r--', 'LineWidth', 1.5);
polarplot(pax, [theta3 theta3],       [-40 0], 'r--', 'LineWidth', 1.5);
pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
pax.ThetaLim = [-90 90]; pax.RLim = [-40 2]; pax.RAxisLocation = 55;
pax.FontSize = 11;
title(pax,'Beampattern (green = SOI, red = jammers)');

% col4 top: 3D complex weights
ax3 = nexttile(tl, 4);
h3 = scatter3(ax3, elem, nan(Nr,1), nan(Nr,1), 90, cmap, 'filled');
grid(ax3,'on'); xlim(ax3,[-0.5 Nr-0.5]); ylim(ax3,[-180 180]); zlim(ax3,[0 1.05]);
xticks(ax3,0:Nr-1); yticks(ax3,-180:90:180);
xlabel(ax3,'Element'); ylabel(ax3,'Phase (deg)'); zlabel(ax3,'|w|');
title(ax3,'Complex weights'); view(ax3,[-35 22]);

% col4 bottom: learning curve
axE = nexttile(tl, 8);
hErr = plot(axE, nan, nan, 'LineWidth', 1.6);
grid(axE,'on'); xlim(axE,[1 N]); ylim(axE,[-15 2]);
xlabel(axE,'Iteration'); ylabel(axE,'Smoothed MSE (dB)');
title(axE,'Learning curve');

sgtitle(fig, 'LMS — press Run', 'FontWeight','bold', 'FontSize',13);

%% ---- LMS loop with STAGED live updates ----
w_lms = zeros(Nr,1);
error_log = zeros(N,1);
for i = 1:N
    rs   = r(:, i);
    err  = soi(i) - w_lms' * rs;            % error uses the CURRENT weights
    error_log(i) = abs(err)^2;
    w_lms = gamma*w_lms + mu*conj(err)*rs;  % ... then weights update (gamma=1: pure LMS)

    if mod(i, refresh) == 0 || i == N
        if ~ishandle(fig), break; end       % stop cleanly if window closed

        % ---- beat step 1: WEIGHTS adjust (amplitude, phase, 3D together) ----
        a  = abs(w_lms); a = a / max(a + eps);   % relative amplitude
        ph = rad2deg(angle(w_lms));
        set(hAmp, 'YData', a);
        set(hPh,  'YData', ph);
        set(h3, 'YData', ph, 'ZData', a);
        sgtitle(fig, sprintf('Iter %d/%d   —   (1) Weights adjust', i, N), ...
                'FontWeight','bold','FontSize',13);
        drawnow; pause(stage_pause);
        if ~ishandle(fig), break; end

        % ---- beat step 2: PATTERN responds ----
        dB = pattern_dB(w_lms); dB = dB - max(dB);   % normalize peak to 0 dB
        set(hPol, 'RData', dB);
        sgtitle(fig, sprintf(['Iter %d/%d   —   (2) Pattern responds   |   ' ...
                'nulls: 60°=%.0f dB, -50°=%.0f dB'], i, N, dB(i2), dB(i3)), ...
                'FontWeight','bold','FontSize',13);
        drawnow; pause(stage_pause);
        if ~ishandle(fig), break; end

        % ---- beat step 3: ERROR measured -> feeds back to weights ----
        smdB = 10*log10(movmean(error_log(1:i), smooth_win));
        set(hErr, 'XData', 1:i, 'YData', smdB);
        sgtitle(fig, sprintf('Iter %d/%d   —   (3) Error measured  ->  feeds back to weights', ...
                i, N), 'FontWeight','bold','FontSize',13);
        drawnow; pause(stage_pause);
    end
end
drawnow;

dB = pattern_dB(w_lms); dB = dB - max(dB);
fprintf('Final:  SOI gain %.1f dB (peak),  null 60° %.1f dB,  null -50° %.1f dB\n', ...
        dB(iS), dB(i2), dB(i3));
