% Revise to use more physically justifiable definition of Grammian when revdef == 1,
% instead of Tomizuka's version. 
%
% Older comments: 
% Attempt to revise to include scaling factors for states, inputs, and
% time, to see how these impact the SVs and to nondimensionalize the SVs.
% With the default Rappel dimensions, the controllability matrices and
% Grammians will have inconsistent units. 
%
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
revdef = 1; % set to nonzero number to use revised ctrb Grammian definition
empirical = 0; % set to nonzero number to use empirical Jacobian coefficients 
estimpflag = 0; % set to something nonzero to use the "newer" parameter set used in observer tests

%     Vd = Verror+V; 
%     nd = nerror+n; 
%    desfile = '1cellOL_testA_b325_3250';
%    eval(['load ' desfile ' Vd nd']);
%dirname = 'data_1cell_b210_smooththeta_ctrb_ol'; % Directory for data file
%dirname = 'data_1cell_b325_smooththeta_ctrb_ol'; % Directory for data file
dirname = 'data_1cell_b225_smooththeta_ctrb_ol'; % Directory for data file
%dirname = 'data_1cell_b325_ctrb_ol'; % Directory for data file
%dirname = 'data'; % Directory for data file
%dirname = 'data_1cell_b325_ctrb_ol_wrong_hp_wrong_acoeff2'; % Directory for data file
%dirname = 'data_1cell_b325_ctrb_ol_wrong_acoeff2'; % Directory for data file
%dirname = 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_jaccompare1e-7';

if isdir(dirname) % used a shorthand mat-file format for initial comparision of analytical/empirical
eval(['load ' dirname '/configinfo'])
else
eval(['load ' dirname ])    
B1 = [-cc_ext*deltat; 0]
end 
%        trueindex = (ii-1)*writeintsteps + k;
%     desfile = '1cellOL_testA_b325'; 
%     eval(['load ' desfile ' Vd nd']);
% The "desired" trajectories are an OL run of the single-cell model at the
% specified BCL. The selected BCL occurs during the transient, but maybe we
% should wait until later (when the model stabilizes). In any case, there
% should be no distinction between Vd, V or nd, n in this case. 

if estimpflag
    Vtmax = 120 %mV; approximate span of an AP
    ntmax = 1.1; %dimensionless; approximate max of n variable
    u1max = 11.3360; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
    u2max = 0.0080755; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)
    display('WARNING: uimax have not been determined for estimator parameter set')
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
if 1 % Controllability Grammians
%    numdiv = 60; % number of subdivisions of a single AP
%    numdiv = 20; % number of subdivisions of a single AP
    numdiv = 5; % number of subdivisions of a single AP
    cgrams= cell(1,numdiv);
    cgramns=cell(1,numdiv);
    svdcgs=cell(1,numdiv);
    svdcgns=cell(1,numdiv);
    svdcgs_plot = zeros(2,numdiv);  % these are the actual values that will be plotted 
    svdcgns_plot = zeros(2,numdiv); 
    Bn = [0;B1(1)]; % B matrix for n as input?
    % scaled matrices
%    B1s = (6.8517e+059)*Tscale*Smat*B1*u1max; % ms / mV * (mV/(ms*cur))*cur % For a window size of 450 steps, increasing B by this amount ensures that minSV of case 1 dominates that of case 2, as predicted (SV's should scale as b1^2)
    B1s = Tscale*Smat*B1*u1max; % ms / mV * (mV/(ms*cur))*cur 
%    Bns = Tscale*Smat*B2*u2max; % ms * (1/(ms)) 
    Bns = Tscale*Smat*Bn*u2max; % ms * (1/(ms)) 
    C = eye(2); % output matrix
    Cs = Smat*C*Smatinv; 
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
%%    init=round(9750.05/deltat); % start time of a selected AP, 325 ms
%    init=round(4200.05/deltat); % start time of a selected AP, 325 ms (and 210?) 
%    init=round(14625.05/deltat); % start time of a selected AP
%    init=round(6750.05/deltat); % start time of a selected AP
%    init=round(1125.05/deltat); % start time of a selected AP
%    init=round(900.05/deltat); % start time of a selected AP, 225ms
%    init=round(1350.05/deltat); % start time of a selected AP
    if estimpflag
        init=1
    else
        init=round(900.05/deltat); % start time of a selected AP, 225ms
%     init=round(4200.05/deltat); % start time of a selected AP, 325 ms
   end
 fileindex = floor(init/writeintsteps)+1
 initindex =  init - (fileindex-1)*writeintsteps; %
 init = initindex; 
% can try to use the above 3 lines if fileindex != 1
%    fileindex = 2; 
if isdir(dirname) 
    eval(['load ' dirname '/' num2str(fileindex)])
end
    %    increment = stimperiod(1)/5; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
%    increment = stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
%    increment = 2*stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
%numpd = 3; 
numpd = 1; % for comparison test
%numpd = 2; 
    increment = numpd*stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
%    windowlength = 4500%5%6500; % interval size over which Grammians will be computed, 225 ms
%    windowlength = 450 % 450*deltat = 22.5ms or 1/10 of a bcl 
%    windowlength = 1125 % 1125*deltat = 56.25ms or 1/4 of a bcl 
%    windowlength = 2250 % 1125*deltat = 112.5ms or 1/2 of a bcl 
%    windowlength = 650 % for 325
%    windowlength = 6500 % 
%    windowlength = 6500*2 % 
%    windowlength = 420 % for 210
%    windowlength = 4200 % 
    if estimpflag
        windowlength = 2875%575% 1/10 bcl = 2875 for initial tests with estimator parameter set 
    else
        windowlength = 4500%450%5%6500; % interval size (number of timesteps) over which Grammians will be computed
    end

tic
    for nn=1:numdiv
        cgram = zeros(2,2);
        cgramn = zeros(2,2);
%        %endstep = numstep-1;
        %       startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        %       endstep = stimperiodsteps(1)+startstep;%round(stimperiodsteps(1)/2)+startstep;%1*stimperiodsteps(1);
        startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        endstep = startstep+windowlength;%+2;%Rappel-units model has very poor scaling, so shouldn't run for very many timesteps to prevent blowup of Grammians
        svdcg=zeros(2,endstep-startstep+1);
        svdcgn=zeros(2,endstep-startstep+1);
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
            stmktoend = eye(2);
            if revdef
                jjrange = (kk+1):endstep; 
            else
                jjrange = kk:endstep-1; 
            end
            for jj = jjrange 
%            for jj = kk:endstep-1 %kk:(endstep) % wrong
                %            stmktoend = reshape(allAr(:,jj),2,2)*stmktoend; %A(numstep-1)...A(kk)
                if empirical 
                    Atmp = reshape(adce(:,jj),2,2); 
                else
                Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
                A12 = [deltat*acoeff2(1,jj)];
                Auc = [1-deltat/tauN];
                A21 = [deltat*acoeff3(1,jj)];
                Atmp = [ Ac   A12
                    A21  Auc];
                end
%                 A = [ Ac   A12
%                     A21  Auc];
                % scale A-matrix
                A = Tscale*Smat*Atmp*Smatinv;
                stmktoend = A*stmktoend; %A(numstep-1)...A(kk)
%                 if isnan(stmktoend)
%                     jj
%                     kk
%                     A
%                    pause
%                 end

            end
            if kk == startstep
                startstep
    Aint = stmktoend
    svd(ctrb(Aint,B1s))
        svd(ctrb(Aint,Bns))
        svd(obsv(Aint,Cs(1,:)))
         svd(obsv(Aint,Cs(2,:)))
     
end
%             cgram = cgram + stmktoend*B1*B1'*stmktoend';
%             cgramn = cgramn + stmktoend*Bn*Bn'*stmktoend';
            cgram = cgram + stmktoend*B1s*B1s'*stmktoend';
            cgramn = cgramn + stmktoend*Bns*Bns'*stmktoend';
            svdcg(:,kk) = svd(cgram); % intermediate svd results
            svdcgn(:,kk) = svd(cgramn);
            if ~mod(kk,1000)
                disp(kk)
            end
            %         if kk==round(endstep/2) % save middle value of grammians
            %             cgrammid = cgram;
            %             cgramnmid = cgramn;
            %         end
        end
        cgrams{nn} = cgram;
        cgramns{nn} = cgramn;
        svdcgs{nn} = svdcg; % intermediate and final svds of Grammians, indexed by interval
        svdcgns{nn} = svdcgn;
        svdcgs_plot(:,nn) = svd(cgrams{nn});  % these are the actual values that will be plotted     
        svdcgns_plot(:,nn) = svd(cgramns{nn});  % these are the actual values that will be plotted     
    end
        toc
    % should probably let this run for several periods. So far, grammian results are
    % consistent with the instantaneous results, but should (1) run for a
    % longer number of bcls (improving the efficiency of the routine would
    % help), and (2) verify the LTV model. In addition, should compute a few
    % steps by hand to confirm that the calculations are sensible.
    %save ctrbgram_poster startstep endstep cgram* svd* stimperiod numpart numstep deltax deltat increment init

    % try to pick an appropriate scale factor for the first plot
    sf=floor(log10(max(svdcgs_plot(2,:)))); 
    trueoffset = (fileindex-1)*writeintsteps;

h1=figure
set(gcf,'Position',[2 100 1280 400])
%subplot(2,1,1)
hold on;
p1=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdcgs_plot(2,:),'-b+');
p2=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdcgns_plot(2,:),'-gs');
set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
%axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^0])
axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^7])
%yl=ylabel('window = 0.1\timesBCL');
yl=ylabel('minimum singular value');
xl=xlabel('start time for Grammian, in BCL indices');
% %    axis([(trueoffset+init)*deltat (trueoffset+windowlength)*deltat 10^-8 10^2])
%     axis([(trueoffset+init)*deltat (trueoffset+endstep)*deltat 10^-8 10^2])
%ti=title('Minimum singular values of controllability Grammians');
legend('Case 1: control V directly','Case 2: control n directly')

if poster
    set(p1,'LineWidth',lw2);
    set(p2,'LineWidth',lw2);
%    set(ti,'FontSize',fs);
    set(gca,'FontSize',fs);
    set(yl,'FontSize',fs);
    set(xl,'FontSize',fs);
end

    figure
    subplot(2,1,1)
    hold on;
    subplot(2,1,2)
    hold on;
    % for nn=1:5
    %     startstep = init+(nn-1)*round(increment/deltat);
    % %    plot(startstep*deltat,svd(cgrams{nn}),'*')
    % %    plot(startstep*deltat,svd(cgramns{nn}),'*')
    %     plot(startstep*deltat,min(svd(cgrams{nn})),'b*')
    %     plot(startstep*deltat,min(svd(cgramns{nn})),'gs')
    % end
    for nn=1:numdiv
%        startstep = init+(nn-1)*round(increment/deltat);
% revise to use absolute time index
        trueoffset = (fileindex-1)*writeintsteps; 
        startstep = trueoffset + init+(nn-1)*round(increment/deltat);
        %    plot(startstep*deltat,svd(cgrams{nn}),'*')
        %    plot(startstep*deltat,svd(cgramns{nn}),'*')
        subplot(2,1,1)
        p1(nn)=plot(startstep*deltat,max(svd(cgrams{nn})),'b*');
        subplot(2,1,2)
        p2(nn)=plot(startstep*deltat,max(svd(cgramns{nn})),'gs');
    end
    subplot(2,1,1)
    yl=ylabel('max singular value');
%    p3=plot((init:endstep)*deltat,5e-23*(V(:,init:endstep)+85)/85,'--');
    legend('control V directly')
    ti=title('maximum singular values of controllability Grammians');
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    subplot(2,1,2)
    yl=ylabel('max singular value');
    legend('control n directly','Location','southeast')
    xl=xlabel('start time for Grammian');
%    p3=plot((init:endstep)*deltat,5*(V(:,init:endstep)+85)/85,'--');
    if poster
        set(gca,'FontSize',fs);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
    end

    
    
    figure
    hold on;
         plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat,svdcgs_plot(1,:)./svdcgs_plot(2,:),'b*');
         plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat,svdcgns_plot(1,:)./svdcgns_plot(2,:),'gs');
         set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
    legend('Case 1: control V directly','Case 2: control n directly')
    p3=plot((trueoffset + (init:endstep))*deltat,30*10^5*(V(:,init:endstep)+85)/85,'--');
    p4=plot((trueoffset + (init:endstep))*deltat,1*10^5*(n(:,init:endstep)),'--');
    yl=ylabel('singular value ratios');
    xl=xlabel('start time for Grammian');
% %    axis([(trueoffset+init)*deltat (trueoffset+windowlength)*deltat 10^1 10^8])
%     axis([(trueoffset+init)*deltat (trueoffset+endstep)*deltat 10^1 10^8])
    ti=title('Condition numbers of controllability Grammians');
    if poster
        set(gca,'FontSize',fs);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end


% save plotted SVs
%    eval(['save ' dirname '/plottedsvs *_plot'])
% save more values for longer run
%    eval(['save ' dirname '/svdslong_uscale_normed svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%    eval(['save ' dirname '/svdslong_allscale_normed svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%    eval(['save ' dirname '/svdslong_uscale_u1_4_u2normed svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%    eval(['save ' dirname '/svdslong_uscale_u1_4_u2normed_range450 svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%    eval(['save ' dirname '/svdslong svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength'])
%eval(['save ' dirname '/svdslong_uscale_u1_4_u2normed_range450_40pts svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range4500_40pts svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range1125_40pts svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range2250_40pts svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range450_60pts_3pd svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range1125_60pts_3pd svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_4_range4500_60pts_3pd svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_11_u2_8em3_range450_60pts_3pd svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_11_u2_8em3_range4500_60pts_3pd svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save ' dirname '/svdslong_allscale_u1_11_u2_8em3_range4200_5pts_2pd_b210 svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save svdscheck_allscale_u1_11_u2_8em3_range450_60pts_3pd_b225 svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save svds_allscale_u1_11_u2_8em3_range2875_5pts_1pd_b230 svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
%eval(['save svds_empirical_allscale_u1_11_u2_8em3_range2875_5pts_1pd_b230 svd* cgram* *_plot numdiv fileindex writeintsteps init increment deltat endstep windowlength Smat Tscale u1max u2max'])
end %if 1

% % plot Vd and V to show they are the same
% reprange=1:numrep;
%     h31=figure;
%     hold on;
%     plot((1:numstep)*deltat,Vd)
%     for ii=reprange % number of data writing cycles
%         eval(['load ' dirname '/' num2str(ii)])        
%         tr=([((ii-1)*writeintsteps+1):1:ii*writeintsteps])*deltat;
%         plot(tr,V(:,1:1:writeintsteps),'g--')
%         clear hprime acoeff* usf Verror nerror iext V n un
%     end

reprange=1:numrep;
    h31=figure;
    hold on;
    for ii=reprange % number of data writing cycles
        eval(['load ' dirname '/' num2str(ii)])        
        tr=([((ii-1)*writeintsteps+1):1:ii*writeintsteps])*deltat;
%        imagesc(deltax*(1:numpart), tr, n(:,1:10:writeintsteps)')
         plot(tr,acoeff1(:,1:1:writeintsteps))
         plot(tr,acoeff2(:,1:1:writeintsteps),'g--')
         plot(tr,acoeff3(:,1:1:writeintsteps),'r:')
         plot(tr,-ones(size(acoeff3(:,1:1:writeintsteps)))/tauN,'c:')
         plot(tr,V(:,1:1:writeintsteps),'k-.')        
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    legend('acoeff1','acoeff2','acoeff3','-1/tauN')
title('Warning: these coefficients are valid for the cts-time 1cell eqs, not the discrete ones')

    h31=figure;
    hold on;
    for ii=reprange % number of data writing cycles
        eval(['load ' dirname '/' num2str(ii)])        
        tr=([((ii-1)*writeintsteps+1):1:ii*writeintsteps])*deltat;
        plot(tr,acoeff3(:,1:1:writeintsteps),'r:')
        plot(tr,-ones(size(acoeff3(:,1:1:writeintsteps)))/tauN,'c:')
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    legend('acoeff3','-1/tauN')
title('Warning: these coefficients are valid for the cts-time 1cell eqs, not the discrete ones')

    % overlay
    h31=figure;
    hold on;
    ii=14; 
        eval(['load ' dirname '/' num2str(ii)])        
        tr=([((1-1)*writeintsteps+1):10:1*writeintsteps])*deltat;
        plot(tr,acoeff2(:,1:10:writeintsteps))
        clear hprime acoeff* usf Verror nerror iext V n un
    ii=15; 
        eval(['load ' dirname '/' num2str(ii)])        
        tr=([((1-1)*writeintsteps+1):10:1*writeintsteps])*deltat;
        plot(tr,acoeff2(:,1:10:writeintsteps), 'g:')
        clear hprime acoeff* usf Verror nerror iext V n un
% Try to check whether scaling of Bn matters. Replacing Bn with [0;B1(1)]
% % will change the grammian by a factor of B1(1)^2
% figure
% hold on;
% for nn=1:5
%     startstep = init+(nn-1)*round(increment/deltat);
% %    plot(startstep*deltat,svd(cgrams{nn}),'*')
% %    plot(startstep*deltat,svd(cgramns{nn}),'*')
%     plot(startstep*deltat,min(svd(cgrams{nn})),'b*')
%     plot(startstep*deltat,min(svd(B1(1)^2*cgramns{nn})),'gs')
% end

% figure
% plot((1:endstep)*deltat,svdcg)
% figure
% plot((1:endstep)*deltat,svdcg(1,:)./svdcg(2,:),'r:')
% figure
% plot((1:endstep)*deltat,svdcgn)
% figure
% plot((1:endstep)*deltat,svdcgn(1,:)./svdcgn(2,:),'r:')

% check a selected value manually
%A0=[0.9778 0; 0 0.9998]
%A1=[0.9821 0; 0 0.9998]
%A2=[0.9864 0; 0.0001 0.9998]
%cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e-004 *
%    0.395979625244837                   0
%                    0                   0
% jj=1
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A0 = [ Ac   A12
%             A21  Auc];
% jj=2
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A1 = [ Ac   A12
%             A21  Auc];
%
% jj=3
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A2 = [ Ac   A12
%             A21  Auc];
% cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e-004 *
%    0.395999903558393   0.000017710454326
%    0.000017710454326   0.000000001335500
% cgram with endstep = 3 (new)
%   1.0e-004 *
%    0.395999903558393   0.000017710454326
%    0.000017710454326   0.000000001335500
% cgram with endstep=3 (old):
%   1.0e-004 *
%    0.385323995616494   0.000037221998733
%    0.000037221998733   0.000000004138828
%
% Redo manual test for Rappel units
%A0 = 1e4*[2.251320187928600  -0.000001496755894; 0.000000000000000 0.000099980000000]
%A1 = 1e4*[2.166800616612064  -0.000001438301991; 0.000000000000000   0.000099980000000]
%A2 = 1e4*[2.084614165392807  -0.000001381583702; 0.000000000000000   0.000099980000000]
% cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e+014 *
%    5.949116066528354                   0
%                    0                   0
% cgram =
%   1.0e+014 *
%    5.949116066528354   0.000000000000000
%    0.000000000000000   0.000000000000000

% Redo for corrected Rappel-units Jacobian, Bn = [0;1]
% A0 =[    0.9901   -0.0000
%     0.0000    0.9998]
% A1 = [    0.9924   -0.0001
%     0.0000    0.9998]
% A2 = [    0.9947   -0.0005
%     0.0000    0.9998]
% cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%     0.0074         0
%          0         0
% For endstep = startstep + 2:
% cgram =
%     0.0074    0.0000
%     0.0000    0.0000