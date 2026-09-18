function beamforming_lab
%% Beamforming Lab — unified interactive harness (ArrayCon26 adaptive-SP course)
%  One app, many methods behind a common interface. Pick a method from the
%  dropdown; run it stationary or against a moving jammer; Store a run and
%  Load past runs into the Comparison tab to compare methods head-to-head.
%
%  Method status:  LMS, MVDR (SMI), Conventional = working.
%                  LCMV, RLS = stubbed (return previous weights) — added next.
%
%  Shared "optimal ceiling" in the suppression panel = ideal MVDR from the TRUE
%  covariance (same bound for every method). Realized suppression rides below it.
%
%  Code-only uifigure. No toolboxes required.  Press Run.

%% ================= scenario (fixed for now) =================
S.d = 0.5;  S.Nr = 8;  S.fs = 1e6;  S.sig2 = 0.5;  S.Aj = 5;   % jammer amplitude (JNR ~17 dB)
S.theta_soi = deg2rad(30);  S.theta3 = deg2rad(-50);   % SOI, fixed jammer
S.f2 = 0.02e6;  S.f3 = 0.03e6;                          % jammer tone freqs
S.aj_far = 72;  S.aj_near = -45;                        % moving-jammer sweep bounds (deg)
S.statAngle = 55;                                       % parked angle in Stationary mode
S.elem = (0:S.Nr-1).';
S.a = @(th) exp(-2j*pi*S.d*S.elem*sin(th));             % steering vector (rad)
S.s1 = S.a(S.theta_soi);  S.s3 = S.a(S.theta3);

% SOI waveform (length-127 Gold code, oversampled x8)
gold = [-1 1 1 -1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 1 1 1 -1 1 ...
        1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 -1 -1 -1 -1 ...
        -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 -1 -1 1 1 1 -1 1 -1 1 1 1 1 1 -1 -1 ...
        -1 -1 1 1 1 -1 1 -1 -1 -1 1 1 1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 1 -1 1 1 ...
        -1 -1 -1 -1 1 -1 1 1 -1 -1 -1 -1 -1 -1 1 1];
S.soi = repelem(gold,8);  S.Nsoi = numel(S.soi);

% beampattern grid + suppression helper
S.Nfft = 1024;
S.theta_bins = asin(linspace(-1,1,S.Nfft));
S.patt = @(w) 10*log10(abs(fftshift(fft([conj(w(:).'), zeros(1,S.Nfft-S.Nr)]))).^2);

% ideal-MVDR "optimal ceiling" over the sweep range (true R)
S.ang_grid = S.aj_near:0.5:S.aj_far;
S.ceiling  = arrayfun(@(ad) ceilingAt(S, deg2rad(ad)), S.ang_grid);
S.suppYmax = max(60, ceil(max(S.ceiling)/10)*10);   % show up to the ideal peak; ticks every 20 dB

% run state
S.Lmax = 256;                       % snapshot ring-buffer length
S.Xbuf = 0.1*(randn(S.Nr,S.Lmax)+1j*randn(S.Nr,S.Lmax));
S.w = zeros(S.Nr,1);  S.P = [];     % P for RLS later
S.Racc = zeros(S.Nr); S.smiN = 0;   % batch accumulator for SMI-MVDR / LCMV
S.rlsDelta = 0.01;                   % RLS init: P = (1/delta) I
S.theta2 = S.aj_far;  S.sdir = -1;  S.nGlob = 0;
S.blk = 40;  S.running = false;
S.trailX = [];  S.trailY = [];      % LMS/realized suppression trail (one sweep)
S.rec = struct('ang',[],'supp',[]); % recorded supp-vs-angle for Store

%% ================= UI =================
S.fig = uifigure('Name','Beamforming Lab','Position',[40 40 1560 860]);
S.fig.CloseRequestFcn = @onClose;
tg = uitabgroup(S.fig,'Units','normalized','Position',[0 0 1 1]);
tabLive = uitab(tg,'Title','Live Lab');
tabCmp  = uitab(tg,'Title','Comparison');

% ---- Live Lab layout ----
gl = uigridlayout(tabLive,[2 4]);
gl.ColumnWidth = {240,'1x','1.4x','1x'};  gl.RowHeight = {'1x','1x'};

% controls panel (col 1, both rows)
cp = uipanel(gl,'Title','Controls');  cp.Layout.Row = [1 2];  cp.Layout.Column = 1;
cg = uigridlayout(cp,[21 1]);  cg.RowHeight = [repmat({'fit'},1,20), {'1x'}];  cg.ColumnWidth = {'1x'};

uilabel(cg,'Text','Method','FontWeight','bold');
S.ddMethod = uidropdown(cg,'Items',{'LMS','MVDR (SMI)','Conventional','LCMV','RLS'}, ...
                        'Value','LMS','ValueChangedFcn',@onMethod);
uilabel(cg,'Text','Mode','FontWeight','bold');
S.ddMode = uidropdown(cg,'Items',{'Moving jammer','Stationary'},'Value','Moving jammer');
S.lblSpeed = uilabel(cg,'Text','Jammer speed (deg/upd)');
S.slSpeed  = uislider(cg,'Limits',[0 2],'Value',0.35);

% --- per-method param controls (shown/hidden by method) ---
S.lblMu = uilabel(cg,'Text','mu = 1.0e-04 (log)');
S.slMu  = uislider(cg,'Limits',[log10(1e-5) log10(2e-3)],'Value',log10(1e-4),'ValueChangedFcn',@onMu);
S.lblL = uilabel(cg,'Text','SMI snapshots L');
S.ddL  = uidropdown(cg,'Items',{'64','256','1024','4096'},'Value','256');
S.lblLoad = uilabel(cg,'Text','diagonal loading = 0.01');
S.slLoad  = uislider(cg,'Limits',[0 0.2],'Value',0.01);
S.lblSoi2 = uilabel(cg,'Text','2nd look dir = 10 deg (LCMV)');
S.slSoi2  = uislider(cg,'Limits',[-80 80],'Value',10,'ValueChangedFcn',@onSoi2);
S.lblLam  = uilabel(cg,'Text','RLS forgetting \lambda = 0.999');
S.slLam   = uislider(cg,'Limits',[0.90 0.9999],'Value',0.999,'ValueChangedFcn',@onLam);
S.lblStub = uilabel(cg,'Text','','FontColor',[0.6 0.2 0.2]);

% buttons + status
S.btnRun   = uibutton(cg,'state','Text','Run','ValueChangedFcn',@onRun);
S.btnReset = uibutton(cg,'Text','Reset','ButtonPushedFcn',@onReset);
S.btnStore = uibutton(cg,'Text','Store run','ButtonPushedFcn',@onStore);
S.txt = uilabel(cg,'Text','','WordWrap','on');

% axes: amplitude (r1c2), phase (r2c2), polar (c3), suppression (c4)
S.cmap = turbo(S.Nr);
S.axA = uiaxes(gl); S.axA.Layout.Row=1; S.axA.Layout.Column=2;
S.hAmp = scatter(S.axA,S.elem,nan(S.Nr,1),60,S.cmap,'filled');
grid(S.axA,'on'); xlim(S.axA,[-0.5 S.Nr-0.5]); ylim(S.axA,[0 1.05]);
title(S.axA,'Weight amplitude'); xlabel(S.axA,'Element'); ylabel(S.axA,'|w| (rel)');

S.axP = uiaxes(gl); S.axP.Layout.Row=2; S.axP.Layout.Column=2;
S.hPh = scatter(S.axP,S.elem,nan(S.Nr,1),60,S.cmap,'filled');
grid(S.axP,'on'); xlim(S.axP,[-0.5 S.Nr-0.5]); ylim(S.axP,[-180 180]);
title(S.axP,'Weight phase'); xlabel(S.axP,'Element'); ylabel(S.axP,'\angle w (deg)');

S.pax = polaraxes(gl); S.pax.Layout.Row=[1 2]; S.pax.Layout.Column=3;
S.hPat = polarplot(S.pax,S.theta_bins,nan(1,S.Nfft),'LineWidth',2); hold(S.pax,'on');
polarplot(S.pax,[S.theta_soi S.theta_soi],[-40 0],'g--','LineWidth',1.5);
polarplot(S.pax,[S.theta3 S.theta3],[-40 0],'r--','LineWidth',1.5);
S.hJam = polarplot(S.pax,[1 1]*deg2rad(S.aj_far),[-40 0],'-','Color',[0.85 0.15 0.55],'LineWidth',2);
S.hSoi2 = polarplot(S.pax,[1 1]*deg2rad(10),[-40 0],'--','Color',[0 0.7 0.7],'LineWidth',1.5,'Visible','off');
S.pax.ThetaZeroLocation='top'; S.pax.ThetaDir='clockwise';
S.pax.ThetaLim=[-90 90]; S.pax.RLim=[-40 2]; S.pax.RAxisLocation=55;
title(S.pax,'Beampattern (green=SOI, red=fixed jam, magenta=moving jam)');

S.axS = uiaxes(gl); S.axS.Layout.Row=[1 2]; S.axS.Layout.Column=4;
plot(S.axS,S.ang_grid,S.ceiling,'Color',[0.5 0.5 0.5],'LineWidth',2); hold(S.axS,'on');
S.hTrail = plot(S.axS,nan,nan,'-','Color',[0.60 0.78 0.95],'LineWidth',1.2);
S.hGap = plot(S.axS,[nan nan],[nan nan],'-','Color',[0.7 0.7 0.7]);
S.hOpt = plot(S.axS,nan,nan,'o','MarkerSize',8,'MarkerEdgeColor',[0.4 0.4 0.4],'LineWidth',1.5);
S.hReal= plot(S.axS,nan,nan,'o','MarkerSize',9,'MarkerFaceColor',[0.85 0.15 0.55],'MarkerEdgeColor','k');
grid(S.axS,'on'); xlim(S.axS,[S.aj_near-3 S.aj_far+3]);
ylim(S.axS,[0 S.suppYmax]); yticks(S.axS,0:20:S.suppYmax);
title(S.axS,'Moving-jammer suppression (magenta) vs angle'); xlabel(S.axS,'Jammer angle (deg)'); ylabel(S.axS,'Suppression (dB)');

% ---- Comparison tab ----
gc = uigridlayout(tabCmp,[1 2]); gc.ColumnWidth = {260,'1x'};
cpanel = uipanel(gc,'Title','Loaded runs'); cpanel.Layout.Column=1;
cpg = uigridlayout(cpanel,[3 1]); cpg.RowHeight={'fit','1x','fit'};
uibutton(cpg,'Text','Load run(s)...','ButtonPushedFcn',@onLoad);
S.lst = uilistbox(cpg,'Items',{},'Multiselect','on');
uibutton(cpg,'Text','Clear','ButtonPushedFcn',@onClearCmp);
S.axCmp = uiaxes(gc); S.axCmp.Layout.Column=2;
plot(S.axCmp,S.ang_grid,S.ceiling,'Color',[0.5 0.5 0.5],'LineWidth',2); hold(S.axCmp,'on');
grid(S.axCmp,'on'); xlim(S.axCmp,[S.aj_near-3 S.aj_far+3]); ylim(S.axCmp,[0 S.suppYmax]); yticks(S.axCmp,0:20:S.suppYmax);
title(S.axCmp,'Comparison — suppression vs angle'); xlabel(S.axCmp,'Jammer angle (deg)');
ylabel(S.axCmp,'Suppression (dB)'); legend(S.axCmp,{'optimal ceiling'},'Location','northwest');
S.cmpData = {};  % stored loaded runs

% timer drives the animation
S.timer = timer('ExecutionMode','fixedSpacing','Period',0.03, ...
                'BusyMode','drop','TimerFcn',@onTick);

onMethod();          % set initial param visibility

%% ================= callbacks =================
    function onRun(src,~)
        if src.Value
            src.Text = 'Pause'; S.running = true;
            if strcmp(S.timer.Running,'off'), start(S.timer); end
        else
            src.Text = 'Run'; S.running = false;
        end
    end

    function onReset(~,~)
        S.w(:) = 0; S.P = []; S.Racc = zeros(S.Nr); S.smiN = 0; S.theta2 = S.aj_far; S.sdir = -1;
        S.trailX = []; S.trailY = []; S.rec = struct('ang',[],'supp',[]);
        S.Xbuf = 0.1*(randn(S.Nr,S.Lmax)+1j*randn(S.Nr,S.Lmax));
    end

    function onMethod(~,~)
        m = S.ddMethod.Value;
        set([S.lblMu S.slMu],'Visible',strcmp(m,'LMS'));
        set([S.lblL S.ddL S.lblLoad S.slLoad],'Visible',any(strcmp(m,{'MVDR (SMI)','LCMV'})));
        set([S.lblSoi2 S.slSoi2],'Visible',strcmp(m,'LCMV'));
        set([S.lblLam S.slLam],'Visible',strcmp(m,'RLS'));
        set(S.hSoi2,'Visible',strcmp(m,'LCMV'));
        onReset();
    end

    function onMu(src,~),   S.lblMu.Text=sprintf('mu = %.1e (log)',10^src.Value); end
    function onSoi2(src,~), S.lblSoi2.Text=sprintf('2nd look dir = %.0f deg (LCMV)',src.Value); end
    function onLam(src,~),  S.lblLam.Text=sprintf('RLS forgetting \\lambda = %.4f',src.Value); end

    function onTick(~,~)
        if ~S.running || ~isvalid(S.fig), return; end
        moving = strcmp(S.ddMode.Value,'Moving jammer');
        if moving, th2 = S.theta2; else, th2 = S.statAngle; end
        s2 = S.a(deg2rad(th2));

        % ---- generate a block of snapshots at current geometry ----
        Xblk = zeros(S.Nr,S.blk); dblk = zeros(1,S.blk);
        for k = 1:S.blk
            S.nGlob = S.nGlob + 1;
            dblk(k) = S.soi(mod(S.nGlob-1,S.Nsoi)+1);
            Xblk(:,k) = S.s1*dblk(k) + S.Aj*s2*exp(2j*pi*S.f2*S.nGlob/S.fs) ...
                      + S.Aj*S.s3*exp(2j*pi*S.f3*S.nGlob/S.fs) + 0.5*(randn(S.Nr,1)+1j*randn(S.Nr,1));
        end
        S.Xbuf = [S.Xbuf(:,S.blk+1:end), Xblk];   % push into ring buffer

        % ---- build context & call the selected module ----
        ctx.Nr=S.Nr; ctx.a_soi=S.s1; ctx.w_prev=S.w; ctx.mu=10^S.slMu.Value;
        ctx.x_block=Xblk; ctx.d_block=dblk; ctx.load=S.slLoad.Value;
        switch S.ddMethod.Value
            case 'LMS'
                S.w = mod_lms(ctx);                        % adaptive: every sample
            case 'RLS'                                     % adaptive: recursive least squares
                lam = S.slLam.Value;
                if isempty(S.P), S.P = (1/S.rlsDelta)*eye(S.Nr); end
                for kk = 1:S.blk
                    u = Xblk(:,kk); pih = S.P*u; gain = pih/(lam + u'*pih);
                    e = dblk(kk) - S.w'*u;  S.w = S.w + gain*conj(e);
                    S.P = (S.P - gain*(u'*S.P))/lam;
                end
            case 'Conventional'
                S.w = mod_conv(ctx);                       % fixed weights
            case {'MVDR (SMI)','LCMV'}                     % batch: solve once per L snapshots
                L = str2double(S.ddL.Value);
                S.Racc = S.Racc + Xblk*Xblk';  S.smiN = S.smiN + S.blk;
                if S.smiN >= L
                    ctx.R_hat = S.Racc / S.smiN;
                    if strcmp(S.ddMethod.Value,'LCMV')
                        ctx.a_soi2 = S.a(deg2rad(S.slSoi2.Value));
                        S.w = mod_lcmv(ctx);
                    else
                        S.w = mod_mvdr(ctx);
                    end
                    S.Racc = zeros(S.Nr);  S.smiN = 0;     % start next batch
                end
        end
        if any(~isfinite(S.w)), S.w(:)=0; end

        % ---- advance moving jammer ----
        if moving
            S.theta2 = S.theta2 + S.sdir*S.slSpeed.Value;
            if S.theta2<=S.aj_near, S.theta2=S.aj_near; S.sdir=+1; S.trailX=[]; S.trailY=[]; end
            if S.theta2>=S.aj_far,  S.theta2=S.aj_far;  S.sdir=-1; S.trailX=[]; S.trailY=[]; end
        end

        % ---- update plots ----
        aw = abs(S.w); aw = aw/max(aw+eps); ph = rad2deg(angle(S.w));
        set(S.hAmp,'YData',aw); set(S.hPh,'YData',ph);
        set(S.hPat,'RData', S.patt(S.w)-max(S.patt(S.w)));
        set(S.hJam,'ThetaData',[1 1]*deg2rad(th2));
        if strcmp(S.ddMethod.Value,'LCMV'), set(S.hSoi2,'ThetaData',[1 1]*deg2rad(S.slSoi2.Value)); end
        supReal = 20*log10(abs(S.w'*S.s1)/(abs(S.w'*s2)+eps));
        supOpt  = interp1(S.ang_grid, S.ceiling, th2, 'linear','extrap');  % ride the drawn curve
        set(S.hReal,'XData',th2,'YData',max(supReal,0));
        set(S.hOpt, 'XData',th2,'YData',max(supOpt,0));
        set(S.hGap, 'XData',[th2 th2],'YData',[max(supReal,0) max(supOpt,0)]);
        S.trailX(end+1)=th2; S.trailY(end+1)=max(supReal,0);
        set(S.hTrail,'XData',S.trailX,'YData',S.trailY);
        S.rec.ang(end+1)=th2; S.rec.supp(end+1)=max(supReal,0);
        S.txt.Text = sprintf('%s | jam %.1f deg | supp %.1f dB | opt %.1f dB | lag %.1f dB', ...
            S.ddMethod.Value, th2, max(supReal,0), max(supOpt,0), max(supOpt-supReal,0));
    end

    function onStore(~,~)
        if isempty(S.rec.ang), uialert(S.fig,'Run first — nothing recorded yet.','Store'); return; end
        run.method = S.ddMethod.Value; run.mode = S.ddMode.Value;
        switch run.method
            case 'LMS',        run.paramStr = sprintf('mu=%.1e',10^S.slMu.Value);
            case 'MVDR (SMI)', run.paramStr = sprintf('L=%s,load=%.2f',S.ddL.Value,S.slLoad.Value);
            case 'LCMV',       run.paramStr = sprintf('L=%s,soi2=%.0f',S.ddL.Value,S.slSoi2.Value);
            case 'RLS',        run.paramStr = sprintf('lam=%.4f',S.slLam.Value);
            otherwise,         run.paramStr = '';
        end
        run.ang = S.rec.ang; run.supp = S.rec.supp; run.stamp = datestr(now,'yyyymmdd_HHMMSS');
        if ~exist('results','dir'), mkdir('results'); end
        fn = fullfile('results',sprintf('%s_%s.mat', strrep(run.method,' ',''), run.stamp));
        save(fn,'run');
        uialert(S.fig,sprintf('Stored: %s',fn),'Store','Icon','success');
    end

    function onLoad(~,~)
        [f,p] = uigetfile('*.mat','Load stored run(s)','MultiSelect','on');
        if isequal(f,0), return; end
        if ischar(f), f={f}; end
        for i=1:numel(f)
            L = load(fullfile(p,f{i})); if ~isfield(L,'run'), continue; end
            r = L.run; lbl = sprintf('%s [%s]',r.method,r.paramStr);
            S.cmpData{end+1} = r;
            plot(S.axCmp, r.ang, r.supp, '.', 'DisplayName', lbl);
            S.lst.Items = [S.lst.Items, {lbl}];
        end
        legend(S.axCmp,'Location','northwest');
    end

    function onClearCmp(~,~)
        cla(S.axCmp); hold(S.axCmp,'on');
        plot(S.axCmp,S.ang_grid,S.ceiling,'Color',[0.5 0.5 0.5],'LineWidth',2,'DisplayName','optimal ceiling');
        xlim(S.axCmp,[S.aj_near-3 S.aj_far+3]); ylim(S.axCmp,[0 S.suppYmax]); yticks(S.axCmp,0:20:S.suppYmax); grid(S.axCmp,'on');
        legend(S.axCmp,'Location','northwest');
        S.cmpData = {}; S.lst.Items = {};
    end

    function onClose(~,~)
        try, stop(S.timer); delete(S.timer); catch, end
        delete(S.fig);
    end
end

%% ================= method modules (common interface: w = mod(ctx)) =================
function w = mod_lms(ctx)
    w = ctx.w_prev;
    for k = 1:size(ctx.x_block,2)
        x = ctx.x_block(:,k);
        e = ctx.d_block(k) - w'*x;
        w = w + ctx.mu*conj(e)*x;      % pure LMS
    end
end

function w = mod_mvdr(ctx)             % SMI-MVDR with diagonal loading
    R = ctx.R_hat + ctx.load*(trace(ctx.R_hat)/ctx.Nr)*eye(ctx.Nr);
    u = R \ ctx.a_soi;
    w = u / (ctx.a_soi' * u);          % distortionless toward SOI
end

function w = mod_conv(ctx)             % conventional (delay-and-sum toward SOI)
    w = ctx.a_soi / ctx.Nr;
end

function w = mod_lcmv(ctx)             % LCMV: unit gain at SOI + 2nd look dir, min output power
    R = ctx.R_hat + ctx.load*(trace(ctx.R_hat)/ctx.Nr)*eye(ctx.Nr);
    C = [ctx.a_soi, ctx.a_soi2];  g = [1;1];
    RiC = R \ C;
    w = RiC * ((C'*RiC) \ g);
end

%% ================= helper: ideal-MVDR suppression at one angle (true R) =========
function supp = ceilingAt(S, th2)
    s2 = S.a(th2);
    R = S.s1*S.s1' + S.Aj^2*(s2*s2' + S.s3*S.s3') + S.sig2*eye(S.Nr);
    w = R \ S.s1;
    supp = max(20*log10(abs(w'*S.s1)/(abs(w'*s2)+eps)), 0);
end
