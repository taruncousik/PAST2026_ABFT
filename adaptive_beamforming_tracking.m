function adaptive_beamforming_tracking
%% Adaptive beamforming — MOVING-JAMMER tracking sandbox (interactive)
%  6-element ULA, one SOI (Gold code) at 30 deg, a FIXED jammer at -50 deg,
%  and a MOVING jammer that sweeps toward the SOI and back. Pure LMS adapts
%  continuously. Live sliders let you explore the tracking tradeoff:
%
%     * Jammer speed  — how fast the moving jammer sweeps in angle.
%     * mu (step size)— the LMS learning rate (log scale).
%
%  Two effects to watch, separated by the two panels:
%     (A) GEOMETRIC LIMIT  — as the jammer nears the SOI their steering
%         vectors align; even the OPTIMAL beamformer cannot null it. The grey
%         "optimal" ceiling in the suppression panel drops toward 0 dB.
%     (B) TRACKING LAG     — a fast jammer + small mu makes the LMS null lag.
%         The blue marker sits BELOW the ceiling; the gap is the lag. Raise mu
%         to close it — but too high and the pattern jitters / diverges.
%
%  Buttons: Pause/Resume, Reset (re-zero weights, restart the sweep).
%  No toolboxes required.

close all;

%% ---- fixed parameters ----
d = 0.5;  Nr = 6;  fs = 1e6;
theta_soi = deg2rad(30);  theta3 = deg2rad(-50);
sig2 = 0.5;                         % per-element noise variance
elem = (0:Nr-1).';
s1 = exp(-2j*pi*d*elem*sin(theta_soi));
s3 = exp(-2j*pi*d*elem*sin(theta3));
f2 = 0.02e6;  f3 = 0.03e6;         % jammer tone frequencies

% sweep bounds for the moving jammer (deg)
aj_far = 72;  aj_near = -45;

% SOI waveform buffer (length-127 Gold code, oversampled)
gold = [-1 1 1 -1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 1 1 -1 1 ...
        1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 -1 -1 -1 -1 ...
        -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 -1 -1 1 1 1 -1 1 -1 1 1 1 1 1 -1 -1 ...
        -1 -1 1 1 1 -1 1 -1 -1 -1 1 1 1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 1 -1 1 1 ...
        -1 -1 -1 -1 1 -1 1 1 -1 -1 -1 -1 -1 -1 1 1];
soi_buf = repelem(gold, 8);  Nsoi = numel(soi_buf);

% beampattern grid
N_fft = 1024;
theta_bins = asin(linspace(-1,1,N_fft));
patt = @(w) 10*log10(abs(fftshift(fft([conj(w(:).'), zeros(1,N_fft-Nr)]))).^2);
iAng = @(th) find(abs(theta_bins-th)==min(abs(theta_bins-th)),1);
iS = iAng(theta_soi);

% optimal (Wiener) weights + suppression for a jammer at angle th2 (rad)
    function [wopt, supp] = optimal(th2)
        s2 = exp(-2j*pi*d*elem*sin(th2));
        R  = s1*s1' + s2*s2' + s3*s3' + sig2*eye(Nr);
        wopt = R \ s1;
        supp = 20*log10( abs(wopt'*s1) / (abs(wopt'*s2)+eps) );   % exact null depth (no FFT-bin snap)
    end

% precompute the optimal "ceiling" curve over the sweep range
ang_grid = aj_near:0.5:aj_far;
sup_ceiling = zeros(size(ang_grid));
for gg = 1:numel(ang_grid)
    [~, sup_ceiling(gg)] = optimal(deg2rad(ang_grid(gg)));
end
sup_ceiling = max(sup_ceiling, 0);

%% ---- figure & manual layout ----
fig = figure('Name','Adaptive beamforming — moving-jammer tracking', ...
             'Color','w','Position',[50 50 1500 830]);

axA = axes('Parent',fig,'Position',[0.05 0.57 0.20 0.35]);   % amplitude
axP = axes('Parent',fig,'Position',[0.05 0.09 0.20 0.35]);   % phase
pax = polaraxes(fig); pax.Position = [0.29 0.09 0.42 0.84];   % polar (big)
axS = axes('Parent',fig,'Position',[0.75 0.09 0.22 0.35]);   % suppression vs angle
pnl = uipanel('Parent',fig,'Units','normalized','Position',[0.735 0.55 0.245 0.40], ...
              'Title','Controls','FontWeight','bold','BackgroundColor','w');

cmap = turbo(Nr);

% amplitude / phase scatters
hAmp = scatter(axA, elem, nan(Nr,1), 80, cmap, 'filled');
grid(axA,'on'); xlim(axA,[-0.5 Nr-0.5]); ylim(axA,[0 1.05]); xticks(axA,0:Nr-1);
xlabel(axA,'Element index'); ylabel(axA,'|w| (relative)'); title(axA,'Weight amplitude');
hPh = scatter(axP, elem, nan(Nr,1), 80, cmap, 'filled');
grid(axP,'on'); xlim(axP,[-0.5 Nr-0.5]); ylim(axP,[-180 180]); xticks(axP,0:Nr-1);
yticks(axP,-180:90:180); xlabel(axP,'Element index'); ylabel(axP,'\angle w (deg)');
title(axP,'Weight phase');

% polar pattern + reference lines (jammer2 line handle updates)
hPat = polarplot(pax, theta_bins, nan(1,N_fft), 'LineWidth', 2); hold(pax,'on');
polarplot(pax, [theta_soi theta_soi], [-40 0], 'g--', 'LineWidth', 1.5);      % SOI
polarplot(pax, [theta3 theta3], [-40 0], 'r--', 'LineWidth', 1.5);            % fixed jammer
hJam = polarplot(pax, [aj_far aj_far]*pi/180, [-40 0], 'm-', 'LineWidth', 2); % moving jammer
pax.ThetaZeroLocation='top'; pax.ThetaDir='clockwise';
pax.ThetaLim=[-90 90]; pax.RLim=[-40 2]; pax.RAxisLocation=55; pax.FontSize=11;
title(pax,'Beampattern   (green=SOI, red=fixed jammer, magenta=moving jammer)');

% suppression-vs-angle panel: grey ceiling + moving markers
hCeil = plot(axS, ang_grid, sup_ceiling, 'Color',[0.5 0.5 0.5], 'LineWidth',2); hold(axS,'on');
xline(axS, -50, 'g--','LineWidth',1.5);
hTrail = plot(axS, nan, nan, '-', 'Color',[0.60 0.78 0.95], 'LineWidth',1.2);  % LMS trail (clears each sweep)
hGap  = plot(axS, [nan nan],[nan nan],'-','Color',[0.7 0.7 0.7]);   % lag gap
hOpt  = plot(axS, nan, nan, 'o','MarkerSize',8,'MarkerEdgeColor',[0.4 0.4 0.4],'LineWidth',1.5);
hLMS  = plot(axS, nan, nan, 'o','MarkerSize',9,'MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k');
grid(axS,'on'); xlim(axS,[aj_near-10 aj_far+10]); ylim(axS,[0 60]);
xlabel(axS,'Moving-jammer angle (deg)   [approaching SOI \rightarrow]');
ylabel(axS,'Jammer suppression (dB)');
title(axS,'Suppression:  grey = optimal ceiling,  blue = LMS');
legend(axS, [hCeil hOpt hLMS], {'optimal ceiling','optimal @ angle','LMS @ angle'}, ...
       'Location','northwest','FontSize',8);

%% ---- controls (sliders / buttons / readouts) ----
uicontrol(pnl,'Style','text','Units','normalized','Position',[0.05 0.86 0.9 0.10], ...
          'String','Jammer speed (deg / update)','BackgroundColor','w','FontSize',9,'HorizontalAlignment','left');
sSpd = uicontrol(pnl,'Style','slider','Units','normalized','Position',[0.05 0.78 0.9 0.08], ...
                 'Min',0,'Max',2,'Value',0.35);
uicontrol(pnl,'Style','text','Units','normalized','Position',[0.05 0.62 0.9 0.10], ...
          'String','mu  (log scale)','BackgroundColor','w','FontSize',9,'HorizontalAlignment','left');
sMu  = uicontrol(pnl,'Style','slider','Units','normalized','Position',[0.05 0.54 0.9 0.08], ...
                 'Min',log10(1e-5),'Max',log10(2e-3),'Value',log10(1e-4));
btnPause = uicontrol(pnl,'Style','togglebutton','Units','normalized','Position',[0.05 0.36 0.42 0.12], ...
                     'String','Pause','FontSize',10,'Callback',@onPause);
uicontrol(pnl,'Style','pushbutton','Units','normalized','Position',[0.53 0.36 0.42 0.12], ...
          'String','Reset','FontSize',10,'Callback',@onReset);
txt = uicontrol(pnl,'Style','text','Units','normalized','Position',[0.05 0.02 0.9 0.30], ...
                'String','','BackgroundColor','w','FontSize',9,'HorizontalAlignment','left');

%% ---- shared run state ----
w      = zeros(Nr,1);
theta2 = aj_far;         % moving-jammer angle (deg)
sdir   = -1;             % sweeping toward SOI first
nGlob  = 0;              % global sample counter (keeps tones phase-continuous)
running = true;
resetReq = false;
blk = 40;                % samples processed per redraw
trailX = []; trailY = [];  % LMS suppression trail for the current sweep

%% ---- main live loop ----
while ishandle(fig)
    if resetReq
        w(:) = 0; theta2 = aj_far; sdir = -1; trailX = []; trailY = []; resetReq = false;
    end
    if running
        mu    = 10^sMu.Value;
        speed = sSpd.Value;
        s2 = exp(-2j*pi*d*elem*sin(deg2rad(theta2)));      % current moving-jammer steering
        for k = 1:blk                                      % run a block of LMS samples
            nGlob = nGlob + 1;
            dsoi = soi_buf(mod(nGlob-1,Nsoi)+1);
            x = s1*dsoi + s2*exp(2j*pi*f2*nGlob/fs) + s3*exp(2j*pi*f3*nGlob/fs) ...
                + 0.5*(randn(Nr,1) + 1j*randn(Nr,1));
            e = dsoi - w'*x;
            w = w + mu*conj(e)*x;                          % pure LMS
        end
        % advance the jammer (bounce between bounds)
        theta2 = theta2 + sdir*speed;
        if theta2 <= aj_near, theta2 = aj_near; sdir = +1; trailX = []; trailY = []; end
        if theta2 >= aj_far,  theta2 = aj_far;  sdir = -1; trailX = []; trailY = []; end

        % ---- update plots ----
        if any(~isfinite(w)), w(:) = 0; end                % guard divergence (too-large mu)
        a = abs(w); a = a/max(a+eps); ph = rad2deg(angle(w));
        set(hAmp,'YData',a); set(hPh,'YData',ph);
        dB = patt(w); dB = dB - max(dB);
        set(hPat,'RData',dB);
        set(hJam,'ThetaData',[1 1]*deg2rad(theta2));
        supLMS = 20*log10( abs(w'*s1) / (abs(w'*s2)+eps) );   % exact null depth (no FFT-bin snap)
        [~, supOpt] = optimal(deg2rad(theta2));
        set(hLMS,'XData',theta2,'YData',max(supLMS,0));
        set(hOpt,'XData',theta2,'YData',max(supOpt,0));
        set(hGap,'XData',[theta2 theta2],'YData',[max(supLMS,0) max(supOpt,0)]);
        trailX(end+1) = theta2; trailY(end+1) = max(supLMS,0);
        set(hTrail,'XData',trailX,'YData',trailY);
        set(txt,'String',sprintf(['mu = %.1e\nspeed = %.2f deg/upd\njammer = %.1f deg\n' ...
              'coherence |s_1^H s_2|/N = %.2f\nLMS supp = %.1f dB\nopt supp = %.1f dB\n' ...
              'tracking lag = %.1f dB'], mu, speed, theta2, ...
              abs(s1'*s2)/Nr, max(supLMS,0), max(supOpt,0), max(supOpt-supLMS,0)));
    end
    drawnow;                                               % also processes slider/button UI
    pause(0.01);
end

%% ---- nested callbacks (share the workspace above) ----
    function onPause(src,~)
        running = ~src.Value;
        if src.Value, src.String = 'Resume'; else, src.String = 'Pause'; end
    end
    function onReset(~,~)
        resetReq = true;
    end
end
