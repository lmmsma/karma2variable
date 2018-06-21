% Try rewriting to compute observability Grammians
% This should be the same as ..._ctrbtest.m, except Bn is changed to
% match the "size" of B1, i.e., Bn(2) = B1(1). This is probably better for
% distinguishing time-step-size effects from the impact of whatever scaling
% coefficient should be "realistically" applied to un(k) if there were some
% way to perturb the n-dynamics in real life.
% Standalone file for computing controllability Grammians.
% Uses Rappel units for two-variable Karma model (exc. M=10 instead of
% M=15).

clear variables
poster = 1; % use flag to make larger plots for poster
fs = 24; % poster fontsize
lw1 = 3; % poster linewidth 1
lw2 = 1.5; % poster linewidth 2
estimpflag = 1; % set to something nonzero to use the "newer" parameter set used in observer tests


%     Vd = Verror+V;
%     nd = nerror+n;
%    desfile = '1cellOL_testA_b325_3250';
%    eval(['load ' desfile ' Vd nd']);
%dirname = 'data_1cell_b325_smooththeta_ctrb_ol'; % Directory for data file
if estimpflag
%    dirname = 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_ol';
%    dirname = 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_1150_smooththeta_ol';
    dirname = 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_920_smooththeta_ol';
else
    dirname = 'data_1cell_b225_smooththeta_ctrb_ol'; % Directory for data file
end
%dirname = 'data_1cell_b325_ctrb_ol'; % Directory for data file
eval(['load ' dirname '/configinfo'])
%        trueindex = (ii-1)*writeintsteps + k;
%     desfile = '1cellOL_testA_b325';
%     eval(['load ' desfile ' Vd nd']);
% The "desired" trajectories are an OL run of the single-cell model at the
% specified BCL. The selected BCL occurs during the transient, but maybe we
% should wait until later (when the model stabilizes). In any case, there
% should be no distinction between Vd, V or nd, n in this case.

% Define scaling factors for states, inputs, and time
if estimpflag
    Vtmax = 120 %mV; approximate span of an AP
    ntmax = 1.1; %dimensionless; approximate max of n variable
else
    Vtmax = 84.9321+13.2296 %mV; approximate span of an AP
    ntmax = 1.0327; %dimensionless; approximate max of n variable
    %u1max = stimheight/deltat/cc_ext; %muA/cm^2; max for current injections (assuming values are saturated)
    %u1max = 4; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
    %u2max = 0.0017654; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)
    u1max = 11.3360; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
    u2max = 0.0080755; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)
end
Smat = diag([1/Vtmax 1/ntmax]); % scaling matrix for error states
%Smat = eye(2);
Smatinv = inv(Smat);
%Tscale = stimperiod(1); % scaling for time derivatives. For cts time it should have units of time, but for discrete time, it should be dimensionless
Tscale = 1; % scaling for time derivatives. For cts time it should have units of time, but for discrete time, it should be dimensionless

% A = Tscale*S*A*Sinv; B1 = Tscale*S*B1*u1max; B2 = Tscale*S*B2*u2max;
if 1 %  Grammians
    %    numdiv = 20; % number of subdivisions of a single AP
%    numdiv = 5; % number of subdivisions of a single AP
    numdiv = 20; % number of subdivisions of a single AP
    ograms= cell(1,numdiv);
    ogramns=cell(1,numdiv);
    svdogs=cell(1,numdiv);
    svdogns=cell(1,numdiv);
    svdogs_plot = zeros(2,numdiv);  % these are the actual values that will be plotted
    svdogns_plot = zeros(2,numdiv);
    %    Bn = [0;B1(1)]; % B matrix for n as input?
    % scaled matrices
%    B1s = Tscale*Smat*B1*u1max; % ms / mV * (mV/(ms*cur))*cur 
%    Bns = Tscale*Smat*B2*u2max; % ms * (1/(ms)) 
    % For C-matrices of the form [c1 0] and [0 c2], the scaling should have
    % no effect
    Cn = Tscale*Smat(2,2)*[0 1]*Smatinv; % C matrix for n as output?
    C1 = Tscale*Smat(1,1)*[1 0]*Smatinv; % C matrix for V as output?

    %    init=round(3.483/deltat); % start time of a selected AP
    %    increment = 0.2; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    %    init=round(450.05/deltat); % start time of a selected AP
    %    init=round(650.05/deltat); % start time of a selected AP
    %    init=round(650.05/deltat); % start time of a selected AP
    %    init=round(0.05/deltat); % start time of a selected AP
    %    init=round(975.05/deltat); % start time of a selected AP
    %    init=round(2275.05/deltat); % start time of a selected AP
    %    init=round(2600.05/deltat); % start time of a selected AP
    %    init=round(2925.05/deltat); % start time of a selected AP
    %    init=round(9750.05/deltat); % start time of a selected AP
    %    init=round(14625.05/deltat); % start time of a selected AP
    %    init=round(6750.05/deltat); % start time of a selected AP
    %    init=round(1125.05/deltat); % start time of a selected AP
    if estimpflag
        init=1
    else
        init=round(900.05/deltat); % start time of a selected AP, 225ms
%     init=round(4200.05/deltat); % start time of a selected AP, 325 ms
   end
    %    init=round(1350.05/deltat); % start time of a selected AP
    fileindex = floor(init/writeintsteps)+1
    initindex =  init - (fileindex-1)*writeintsteps; %
    init = initindex;
    % can try to use the above 3 lines if fileindex != 1
    %    fileindex = 2;
    eval(['load ' dirname '/' num2str(fileindex)])
    %    increment = stimperiod(1)/5; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    %    increment = stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    numpd = 2;
    increment = numpd*stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    if estimpflag
        windowlength = 2875%28750%575%%575=1/50 of bcl at deltat = .008 and bcl=230ms, 2875=1/10 
    else
        windowlength = 450%4500%5%6500; % interval size (number of timesteps) over which Grammians will be computed
    end
    tic
    for nn=1:numdiv
        ogram = zeros(2,2);
        ogramn = zeros(2,2);
        %        %endstep = numstep-1;
        %       startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        %       endstep = stimperiodsteps(1)+startstep;%round(stimperiodsteps(1)/2)+startstep;%1*stimperiodsteps(1);
        startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        endstep = startstep+windowlength;%+2;%Rappel-units model has very poor scaling, so shouldn't run for very many timesteps to prevent blowup of Grammians
        svdog=zeros(2,endstep-startstep+1);
        svdogn=zeros(2,endstep-startstep+1);
        %For start time 9750.05ms, varying the number of steps in the Grammian didn’t seem to affect the qualitative result:
        % For Bn = [0;1], max(svdcgns_plot(2,:))
        % ... endstep = startstep+1
        % ... startstep+2
        % ... startstep+3
        % 1.6623 ... startstep+4
        % 1.9049 ... startstep+5
        % 2.1495 ... startstep+6
        % 3.1361 ... startstep+10
        % max(svdcgs_plot(2,:))
        % 1.4716e-029... startstep+4
        % 2.2191e-029 ... startstep+5
        % 3.0718e-029 ... startstep+6
        % 1.6552e-028 ... startstep+10
        % For a later start time (5 step Grammian, Bn = [0;1]), 9750.05ms,
        % max(svdcgs_plot(2,:))
        % 2.2191e-029
        % max(svdcgns_plot(2,:))
        % 1.9049
        % min(svdcgs_plot(1,:))
        % 0.0116
        % min(svdcgns_plot(1,:))
        % 5.9940
        for kk=startstep:endstep
            % original version was wrong
            % Define state transition matrix
            stmstarttok = eye(2);
            %The computation of transition matrices is different for obsv grammians.
            %For ctrb, we want products like Afinal...Akk, where kk varies from final
            %to 0 (the endpoint is always the same but the initial time varies). For
            %obsv, we want products like Akk...A0, where the initial time is fixed but
            %the endpoint varies.
            %            for jj = kk:endstep-1 % ctrb
            for jj = startstep:kk-1 % obsv
                %            stmktoend = reshape(allAr(:,jj),2,2)*stmktoend; %A(numstep-1)...A(kk)
                Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
                A12 = [deltat*acoeff2(1,jj)];
                Auc = [1-deltat/tauN];
                A21 = [deltat*acoeff3(1,jj)];
                Atmp = [ Ac   A12
                    A21  Auc];
                %                 A = [ Ac   A12
                %                     A21  Auc];
                % scale A-matrix
                A = Tscale*Smat*Atmp*Smatinv;
                stmstarttok = A*stmstarttok; %A(numstep-1)...A(kk)
                %                 if isnan(stmktoend)
                %                     jj
                %                     kk
                %                     A
                %                    pause
                %                 end
            end

            %            cgram = cgram + stmktoend*B1*B1'*stmktoend';
            %           cgramn = cgramn + stmktoend*Bn*Bn'*stmktoend';
            ogram = ogram + stmstarttok'*C1'*C1*stmstarttok;
            ogramn = ogramn + stmstarttok'*Cn'*Cn*stmstarttok;
            svdog(:,kk) = svd(ogram); % intermediate svd results
            svdogn(:,kk) = svd(ogramn);
            if ~mod(kk,1000)
                disp(kk)
            end
            %         if kk==round(endstep/2) % save middle value of grammians
            %             cgrammid = cgram;
            %             cgramnmid = cgramn;
            %         end
        end
        ograms{nn} = ogram;
        ogramns{nn} = ogramn;
        svdogs{nn} = svdog; % intermediate and final svds of Grammians, indexed by interval
        svdogns{nn} = svdogn;
        svdogs_plot(:,nn) = svd(ograms{nn});  % these are the actual values that will be plotted
        svdogns_plot(:,nn) = svd(ogramns{nn});  % these are the actual values that will be plotted
    end
    toc
    % should probably let this run for several periods. So far, grammian results are
    % consistent with the instantaneous results, but should (1) run for a
    % longer number of bcls (improving the efficiency of the routine would
    % help), and (2) verify the LTV model. In addition, should compute a few
    % steps by hand to confirm that the calculations are sensible.
    %save ctrbgram_poster startstep endstep cgram* svd* stimperiod numpart numstep deltax deltat increment init
end
% try to pick an appropriate scale factor for the first plot
sf=floor(log10(max(svdogs_plot(2,:))));

poster = 1; % use flag to make larger plots for poster
fs = 24; % poster fontsize
lw1 = 3; % poster linewidth 1
lw2 = 1.5; % poster linewidth 2
% *************** comment back in for standalone script
if 0
    dirname = 'data_1cell_b225_smooththeta_ctrb_ol'; % Directory for data file
    eval(['load ' dirname '/configinfo'])
    %eval(['load ' dirname '/svdslong_allscale_u1_4_range450_60pts_3pd'])
    eval(['load ' dirname '/svdslong_allscale_u1_11_u2_8em3_range450_60pts_3pd'])

    fileindex = floor(init/writeintsteps)+1
    initindex =  init - (fileindex-1)*writeintsteps; %
    init = initindex;
    eval(['load ' dirname '/' num2str(fileindex)])
end

minV=min(V(:,init:endstep));
trueoffset = (fileindex-1)*writeintsteps;

h1=figure
set(gcf,'Position',[2 100 1280 400])
%subplot(2,1,1)
hold on;
p1=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogs_plot(2,:),'-b+');
p2=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogns_plot(2,:),'-gs');
set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
%axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^0])
axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^7])
%yl=ylabel('window = 0.1\timesBCL');
yl=ylabel('minimum singular value');
xl=xlabel('start time for Grammian, in BCL indices');
% %    axis([(trueoffset+init)*deltat (trueoffset+windowlength)*deltat 10^-8 10^2])
%     axis([(trueoffset+init)*deltat (trueoffset+endstep)*deltat 10^-8 10^2])
%ti=title('Minimum singular values of controllability Grammians');
legend('Case 1: observe V directly','Case 2: observe n directly')

if poster
    set(p1,'LineWidth',lw2);
    set(p2,'LineWidth',lw2);
    %    set(ti,'FontSize',fs);
    set(gca,'FontSize',fs);
    set(yl,'FontSize',fs);
    set(xl,'FontSize',fs);
end
%tl=text(6.88,0.01,'(a)')
%tl=text(6.88,10^5,'(a)')
tl=text(3.85,10^8,'(a)')
set(tl,'FontSize',fs)
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2+7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^3,'--');
p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 10)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^12,':');

%eval(['save ' dirname '/svds_allscale_5pts_2pd svd* ogram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale'])

% *************** comment back in for long time scales
if 0
    eval(['load ' dirname '/svdslong_allscale_u1_11_u2_8em3_range4500_60pts_3pd'])

    h2=figure
    set(gcf,'Position',[2 100 1280 400])
    %subplot(2,1,2)
    hold on;
    p1=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdcgs_plot(2,:),'-b+');
    p2=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdcgns_plot(2,:),'-gs');
    set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
    xlim=get(gca,'XLim');
    ylim=get(gca,'YLim');
    %p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 - 7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^7,'--');
    %p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 3)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^30,'--');
    p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2+7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^3,'--');
    p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 10)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^12,':');
    %axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^0])
    axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^7])
    %yl=ylabel('window = 1\timesBCL');
    yl=ylabel('minimum singular value');
    %xl=xlabel('start time for Grammian');
    xl=xlabel('start time for Grammian, in BCL indices');
    if poster
        set(p1,'LineWidth',lw2);
        set(p2,'LineWidth',lw2);
        %   set(ti,'FontSize',fs);
        set(gca,'FontSize',fs);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
    end
    %tl=text(6.88,10^5,'(b)')
    tl=text(3.85,10^8,'(b)')
    set(tl,'FontSize',fs)
end

