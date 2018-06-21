% Compute eigenvalues of time-integrated Karma model for different observer
% gain values. Modified to test controllabilty and observability.
% July 2011: Changed dirder_mod to dirder_normscale in call to "eigs",
% since, presumably, eigs doesn't restrict itself to perturbing one state
% variable at a time. Oddly, the results are similar either way? 
clear variables
% Set observer feedback, number of cells, and bcl
%allgains(:,:,1) = [-1; 0];
%allgains(:,:,1) = [0; 0.001];
%allgains(:,:,1) = [0; -0.0005];
%for i=2:2%2:6%2:13
%%%    allgains(:,:,i) = allgains(:,:,i-1)+[0.1;0];
%    allgains(:,:,i) = allgains(:,:,i-1)+[0;0.001];
%end


if 0
    Kpgains = -0.9:0.1:-0.7;
    %Kpngains = 0.002:0.001:0.004;
    %Kpngains = 0.00005:0.00005:0.00015;
    %Kpngains = [0.0002 0.0005 0.0008];
    Kpngains = [0.0004];
    tempgains = zeros(2*length(Kpngains),length(Kpgains));
    tempgains(1:2:(2*length(Kpngains)),:) = repmat(Kpgains,length(Kpngains),1);
    tempgains(2:2:(2*length(Kpngains)),:) = repmat(Kpngains,length(Kpngains),1)';
    tempgains=reshape(tempgains,2*length(Kpgains)*length(Kpngains),1);
    %tempgains = zeros(2*length(Kpgains),length(Kpngains));
    %tempgains(1:2:(2*length(Kpgains)),:) = repmat(Kpgains,length(Kpgains),1);
    %tempgains(2:2:(2*length(Kpgains)),:) = repmat(Kpngains,length(Kpngains),1)';
    %tempgains=reshape(tempgains,2*length(Kpgains)*length(Kpgains),1);
    for i = 1:length(Kpgains)*length(Kpngains)
        allgains(:,:,i) = tempgains((2*i - 1):(2*i),:);
    end
    allgains(:,:,i+1) = [0;0]; % dummy value (won't be used in loop below)
end

if 0
    allgains(:,:,1) = [-1.5 0; 0 -1.5; 0 0; 0 0];
    for i=2:15
        allgains(:,:,i) = allgains(:,:,i-1) + 0.1*[1 0; 0 1; 0 0; 0 0];
        %    allgains(:,:,i) = allgains(:,:,i-1)+[0;0.001];
    end
    
    %i=1
    %allgains(:,:,i) = [-0.9 0; 0 -0.9; 0 0; 0 0];
    allgains(:,:,i+1) = zeros(size(allgains(:,:,i)));
    allgains(:,:,i+2) = zeros(size(allgains(:,:,i)));
end

if 0
    %allgains(:,:,1) = [-0.3;0];
    allgains(:,:,1) = [-0.35;0];
    for i=2:8%2:3
        allgains(:,:,i) = allgains(:,:,i-1)+[0.05;0];
    end
    allgains(:,:,i+1) = [0;0]; % dummy value (won't be used in loop below)
end

if 0 % 2-cells with default parameters and bcl=160
    gainrange = -0.35:0.05:0;
    %allgains = zeros(4,2,length(gainrange)+1);
    allgains(:,:,1) = [gainrange(1) 0; 0 gainrange(1); 0 0; 0 0];
    for i=2:length(gainrange)
        allgains(:,:,i) = [gainrange(i) 0; 0 gainrange(i); 0 0; 0 0];
    end
    allgains(:,:,i+1) = zeros(size(allgains(:,:,i)));
    allgains(:,:,i+2) = zeros(size(allgains(:,:,i)));
end

%L = squeeze(allgains(:,:,1));
numstate1cell = 2;
%numpart = 1; % number of cells
%numpart = 2;
numpart = 106;
numstate = numpart*numstate1cell;
bcl = 230; %bcl, ms
%bcl = 160; %ms
%bcl = 180; %ms
%bcl = 200; %ms
L = zeros(numstate,numpart); % observer feedback gain
%load L_106cell_b230_des0p85_pwrit2 Lnsoli
%L(:,16) = Lnsoli; 
load L_106cell_b230_des-0p85_pwrit2_V16n16 Lnsoli
L(16,16) = Lnsoli(1);
L(numpart+16,16) = Lnsoli(2); 
%load macresults\kfgain_Q2_R1_C16_b200_stimstart1 Lkf
%L(:,16) = Lkf; 
%L(1,1) = -0.8
%L(2,1) = 0.008
%  L(16,16) = -100;
%  L(16,16) = -10;
%  L(16,16) = -1;
%  L(16,16) = -.1;
%  L(16,16) = -25;
%  L(16,16) = -75;
%  L(16,16) = -125;
%  L(46,46) = -.1;
%  L(46,46) = -1;
%  L(46,46) = -10;
%  L(46,46) = -25;
%  L(46,46) = -75;
%  L(46,46) = -200;
%  L(46,46) = -125;
%  L(76,76) = -.1;
%  L(76,76) = -10;
%  L(76,76) = -25;
%  L(76,76) = -75;
%  L(31,31) = -25;
%  L(61,61) = -25;
%  L(91,91) = -25;
% 106-cell fiber detector locations: 16 31 46 61 76 91
% for i=1:106
%     L(i,i) = -10;
% end
%  L((-3:3)+16+numpart,16) = .00002;
%  L((-3:3)+46+numpart,46) = .00002;
%  L((-3:3)+76+numpart,76) = .00002;
%  L((-3:3)+31+numpart,31) = .00002;
%  L((-3:3)+61+numpart,61) = .00002;
%  L((-3:3)+91+numpart,91) = .00002;
%L(i,j) is feedback from measured V at cell j being applied to state i

deltat = .008; % ms
%stimstart = deltat; % stimulus start time, ms
%stimstart = deltat;%+184; % stimulus start time, ms
%0:46:230
%allstimstart= deltat+[0    46    92   138   184   230];
%allstimstart= deltat+[46    92   138   184   230];
%allstimstart= deltat+[92   138   184   230];
allstimstart= deltat;
stimstart = allstimstart(1);

save kseparams L numpart bcl stimstart % these will be read in by karma_sim_estimtest_p2p.m
epsln = 1e-5; % relative perturbation size for use with diffjac_mod
nsolitol = 1e-12; % abs & rel tolerance for Newton-Krylov solver
useolfp = 1; % set to zero to always solve for fixed point; set to 1 to re-use OL fixed point as CL fixed point
maxabsclerr = NaN; % maximum absolute value of p2p error (use to determine whether OL f.p. is a valid CL f.p.)
is2to1 = 0; % set to 1 if 2:1 block pattern
is2to2 = 0; % set to 1 if 2:2 pattern
computejac = 0; % set to 1 to compute Jacobian explicitly
if numpart == 1
    computejac = 1;
end
%fpfolder = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_smooththeta_ol'; % folder containing trajectory that passes through fixed point
%fpfolder = 'edata_106cell_b200_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_smooththeta_ol'; % folder containing trajectory that passes through fixed point


% If computing Jacobians of OL system, can assess ctrb/obsv.
% input matrix
B = deltat*eye(numstate); % This isn't right, just a placeholder. At least in the LTV system, B isn't square and its nonzero elements aren't 1.
% output matrix
C = eye(numstate); % This is also a placeholder. There is a dimensional inconsistency with L above. If L is numstate by numpart, it is sized to fit the assumption that only V's can be measured. C should technically have the transposed dimensions of L.
Vtmax = 120; %mV; approximate span of an AP
ntmax = 1.1; %dimensionless; approximate max of n variable
Smat = diag([(1/Vtmax)*ones(1,numpart) (1/ntmax)*ones(1,numpart)]); % Scaling matrix: xbar = Smat x is the scaled state. Choose diagonal elements of S so that elements of xbar have similar amplitudes (e.g. Sii = 1/|xi_max|)
Smatinv = inv(Smat); % only need to compute once
u1max = 97; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
u2max = 0.065; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)
%umax = diag(ones(1,size(B,2))); % This is also incorrect and just a placeholder. For inputs, also need to know approximate maximum values of each element of input vector. Note that u is a deviational quantity that may or may not attain the size of the stimulus.
%umax(1,1) = u1max;
%umax(2,2) = u2max;
umax = diag([u1max*ones(1,size(B,2)/2) u2max*ones(1,size(B,2)/2)]);
% Each umax diagonal element should be size of deviational max for integrated system.
Bs = Smat*B*umax; % scaled B matrix
Cs = Smat*C*Smatinv; % scaled B matrix

rankcutoff = 1e-14; % below this level, singular values don't contribute to the rank


for i=1:length(allstimstart)%1%(size(allgains,3)-1)
    clear jacback jaccd
    % Intitalize svd matrices
    svdctrb = zeros(numstate);
    svdctrb_scaled = zeros(numstate);
    svdobsv = zeros(numstate);
    svdobsv_scaled = zeros(numstate);
    %evc = cell(numstate);
    %evo = cell(numstate);
    %evnc = cell(numstate);
    %evno = cell(numstate);
    evcf = cell(numstate);
    evof = cell(numstate);
    evncf = cell(numstate);
    evnof = cell(numstate);
    rankc = zeros(1,numstate);
    ranko = zeros(1,numstate);
    rankcf = zeros(1,numstate);
    rankof = zeros(1,numstate);
    %maxeigncf = zeros(1,numstate);
    %maxeignof = zeros(1,numstate);
    %mineigncf = zeros(1,numstate);
    %mineignof = zeros(1,numstate);
    evcf_scaled = cell(numstate);
    evof_scaled = cell(numstate);
    evncf_scaled = cell(numstate);
    evnof_scaled = cell(numstate);
    rankc_scaled = zeros(1,numstate);
    ranko_scaled = zeros(1,numstate);
    rankcf_scaled = zeros(1,numstate);
    rankof_scaled = zeros(1,numstate);
    
    if numpart == 1 % Set to 1 to turn on automatic filename generation (useful for single-cell tests)
        % Select file name for saving results
        % Replace decimal "p", otherwise it can cause problems loading files:
        for ii = 1:size(L,1)
            for jj = 1:size(L,2)
                temp = num2str(L(ii,jj)); % this can probably be vectorized but I'm not sure how (direct num2str will pad with spaces)
                kk = strfind(temp,'.'); % find decimal
                if kk
                    temp(kk) = 'p';
                end
                gainstr{ii}{jj} = temp;
            end
        end
        diffstr = '';
        % single-cell:
        fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        % This naming convention will only work for single cells. Not sure what to do to summarize multicell feedback
    else  % Provide a filename suffix
        %        fnamesuffix = input('Enter filename suffix in single quotes, or the number zero to leave blank. The default suffix is ''Kp0_Kpn0'': ')
        fnamesuffix = 'Kp0_Kpn0';
        diffstr = '_d005'; % indicate diffusion coeff
        if computejac
            fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix];
            fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix];
            fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_' fnamesuffix];
            if is2to1
                fname1 = [fname1 '_2bcl'];
                fname2 = [fname2 '_2bcl'];
                fname3 = [fname3 '_2bcl'];
            elseif is2to2
                fname1 = [fname1 '_alt'];
                fname2 = [fname2 '_alt'];
                fname3 = [fname3 '_alt'];
            end
        else
            fname4 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_' fnamesuffix]; % for use with eigs
            if is2to1
                fname4 = [fname4 '_2bcl'];
            elseif is2to2
                fname4 = [fname4 '_alt'];
                
            end
        end
    end
    
    %    % multi-cell:
    %    if numpart > 1
    %        diffstr = '_d005'; % indicate diffusion coeff
    %    else
    %        diffstr = '';
    %    end
    %     %   fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
    %     %   fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
    %     %   fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
    %     %   fname4 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert' num2str(epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}]; % for use with eigs
    %     %       fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1} '_2bcl'];
    %     %       fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1} '_2bcl'];
    %     %       fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1} '_2bcl'];
    %     fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_' fnamesuffix];
    %     fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_' fnamesuffix];
    %     fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_' fnamesuffix];
    %     fname4 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert' num2str(epsln,'%1.e') '_' fnamesuffix]; % for use with eigs
    
    if ~useolfp
        tstartfp=tic;
        % Find closed-loop fixed point (ideally it should be identical to the open-loop fixed point):
        %    [sol, it_hist, ierr, x_hist] = nsoli([-90.4237;0.7456],'karma_sim_estimtest_p2pzero',[1e-12 1e-12])
        if exist('fpfolder')
            eval(['load ' fpfolder '/1 V n']);
            ivshiftindex = bcl/deltat - stimstart/deltat + 2;
            if ivshiftindex > bcl/deltat
                ivshiftindex = ivshiftindex - bcl/deltat;
            end
            initialvalue = [V(:,ivshiftindex);n(:,ivshiftindex)];
            [sol, it_hist, ierr, x_hist] = nsoli(initialvalue,'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, searching for 1:1 steady-state
        else
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.4237,numpart,1);repmat(0.7456,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, searching for 1:1 steady-state
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.47,numpart,1);repmat(0.9058,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, stimstart = 46+deltat
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(8.252,numpart,1);repmat(1.119,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, stimstart = 92+deltat
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(31.89,numpart,1);repmat(1.027,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, stimstart = 138+deltat
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(34.35,numpart,1);repmat(0.9045,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=230, stimstart = 184+deltat
            %    [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.4885,numpart,1);repmat(0.5348,numpart,1)],'karma_sim_estimtest_p2pzero',[1e-12 1e-12]) % Approx IC for default parameters, bcl=160, searching for 2:1 steady-state
            [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.474,numpart,1);repmat(0.88,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=160, searching for 1:1 steady-state
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.488,numpart,1);repmat(0.57,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=180, searching for 2:2 steady-state (this yields 2:1 traj)
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.36,numpart,1);repmat(1.05,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=180, searching for 2:2 steady-state (this yields 2:1 traj)
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.5,numpart,1);repmat(1.05,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=180, searching for 2:2 steady-state (this yields 2:1 traj)
            %        [sol, it_hist, ierr, x_hist] = nsoli([repmat(-90.488,numpart,1);repmat(0.55,numpart,1)],'karma_sim_estimtest_p2pzero',[nsolitol nsolitol]) % Approx IC for default parameters, bcl=200, searching for 1:1 steady-state (led to an aphysical solution)
        end
        %   load edata_1cell_b160_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_cdiff_Kp0_Kpn0_2bcl.mat
        %    [sol, it_hist, ierr, x_hist] = nsoli([repmat(sol(1),numpart,1);repmat(sol(2),numpart,1)],'karma_sim_estimtest_p2pzero',[1e-12 1e-12])
        tendfp = toc(tstartfp)
        clerr = karma_sim_estimtest_p2pzero(sol); % Originally only meant to compute this when skipping fixed point ID, but it is useful in general
        maxabsclerr = max(abs(clerr))
    else
                eval(['load edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_Kp0_Kpn0 sol ierr x_hist it_hist']);  % if not including 'start' in filename
%                eval(['load edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_cdiff_Kp0_Kpn0 sol ierr x_hist it_hist']); % use eig instead of eigs
%        eval(['load edata_' num2str(numpart) 'cell_b' num2str(bcl) diffstr '_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert' num2str(epsln,'%1.e') '_start' num2str(round(stimstart/deltat)) '_Kp0_Kpn0 sol ierr x_hist it_hist']); % use eigs instead of eig
        clerr = karma_sim_estimtest_p2pzero(sol);
        maxabsclerr = max(abs(clerr))
        if maxabsclerr > nsolitol
            disp('OL fixed point is not a CL fixed point')
            break
        end
    end
    
    
    if computejac
        tstartjacfw=tic;
        diffjac_mod(sol,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln); % Solve for empirical bcl-to-bcl Jacobian. 1e-5 (searching by decades) seemed to give the closest answer to the 1-cell analytical result.
        tendjacfw=toc(tstartjacfw)
        load jac
        [v,d] = eig(jac);
        eval(['save ' fname1 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist useolfp clerr maxabsclerr stimstart'])
        %save edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_relpert1e-5_Kp-0p8 *
        
        tstartjacbk=tic;
        diffjac_mod(sol,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),-epsln); % Solve for empirical bcl-to-bcl Jacobian using perturbation of opposite sign
        tendjacbk=toc(tstartjacbk)
        load jac
        [v,d] = eig(jac);
        eval(['save ' fname2 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist useolfp clerr maxabsclerr stimstart'])
        %save edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_relpert-1e-5_Kp-0p8 *
        jacback=jac;
        
        eval(['load ' fname1 '.mat jac'])
        jaccd = (jac+jacback)/2;
        [v,d] = eig(jaccd); % eigenvalues of "central difference" Jacobian
        norm(diag(d))
        eval(['save ' fname3 ' L bcl d epsilon epsln ierr it_hist jac numpart sol v x_hist jacback jaccd useolfp clerr maxabsclerr stimstart'])
        
        %save edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_relpert1e-5_cdiff_Kp-0p8 *
        if ~L
            % compute singular values of controllability and observability
            % matrices for specified B and C for both scaled and unscaled
            % systems
            rankcutoffdefault = max(size(jaccd))*eps(norm(jaccd));% This is Matlab's default rank cutoff. Can't easily control cutoffs in orth /null, so to interpret these results, would need to know default
            for kb = 1:size(B,2)
                P = ctrb(jaccd,B(:,kb));
                Ps = ctrb(Smat*jaccd*Smatinv,Bs(:,kb));
                %                 M = [orth(P) null(P)];
                %                 Abc = inv(M)*jaccd*M;
                %                 if ~sum(sum(isnan(Abc)))
                %                 evc{kb} = eig(Abc(1:size(orth(P),2),1:size(orth(P),2)),'nobalance');
                %                 evnc{kb} = eig(Abc((size(orth(P),2)+1):end,(size(orth(P),2)+1):end),'nobalance');
                %                 end
                svdctrb(:,kb) = svd(P);
                svdctrb_scaled(:,kb) = svd(Ps);
                rankc(kb) = sum(svdctrb(:,kb)>=rankcutoff);
                rankc_scaled(kb) = sum(svdctrb_scaled(:,kb)>=rankcutoff);
                [Abar,Bbar,Cbar,T,k] = ctrbf(jaccd,B(:,kb),C(1,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankcf(kb) = sum(k>0);
                rankcf(kb) = sum(k);
                evcf{kb} = eig(Abar((numstate-rankcf(kb)+1):end,(numstate-rankcf(kb)+1):end),'nobalance'); % eigenvalues of Ac
                evncf{kb} = eig(Abar(1:(numstate-rankcf(kb)),1:(numstate-rankcf(kb))),'nobalance'); % eigenvalues of Anc
                [Abars,Bbars,Cbars,Ts,ks] = ctrbf(Smat*jaccd*Smatinv,Bs(:,kb),Cs(1,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankcf_scaled(kb) = sum(ks>0);
                rankcf_scaled(kb) = sum(ks);
                evcf_scaled{kb} = eig(Abars((numstate-rankcf_scaled(kb)+1):end,(numstate-rankcf_scaled(kb)+1):end),'nobalance'); % eigenvalues of Ac
                evncf_scaled{kb} = eig(Abars(1:(numstate-rankcf_scaled(kb)),1:(numstate-rankcf_scaled(kb))),'nobalance'); % eigenvalues of Anc
            end
            for kc = 1:size(C,1)
                Q = obsv(jaccd,C(kc,:));
                Qs = obsv(Smat*jaccd*Smatinv,Cs(kc,:));
                %                 Q = ctrb(jaccd',C(kc,:)'); %I'm using the transposed version since I wasn't sure how to get orth/null into proper row form (transposing them doesn't work)
                %                 O = [orth(Q) null(Q)];
                %                 Abo = inv(O)*jaccd'*O;
                %                 if ~sum(sum(isnan(Abo)))
                %                 evo{kc} = eig(Abo(1:size(orth(Q),2),1:size(orth(Q),2)),'nobalance');
                %                 evno{kc} = eig(Abo((size(orth(Q),2)+1):end,(size(orth(Q),2)+1):end),'nobalance');
                %                 end
                svdobsv(:,kc) = svd(Q);
                svdobsv_scaled(:,kc) = svd(Qs);
                ranko(kc) = sum(svdobsv(:,kc)>=rankcutoff);
                ranko_scaled(kc) = sum(svdobsv_scaled(:,kc)>=rankcutoff);
                [Abar,Bbar,Cbar,T,k] = obsvf(jaccd,B(:,1),C(kc,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                % The rank of the observability matrix is supposed to be sum(k)? The result
                % doesn't seem to agree with that of the standard rank computation
                %                rankof(kc) = sum(k>0);
                rankof(kc) = sum(k);
                evof{kc} = eig(Abar((numstate-rankof(kc)+1):end,(numstate-rankof(kc)+1):end),'nobalance'); % eigenvalues of Ao
                evnof{kc} = eig(Abar(1:(numstate-rankof(kc)),1:(numstate-rankof(kc))),'nobalance'); % eigenvalues of Ano
                [Abars,Bbars,Cbars,Ts,ks] = obsvf(Smat*jaccd*Smatinv,Bs(:,1),Cs(kc,:),rankcutoff); % compute staircase form assuming B = I(:,1) and C = I(1,:)
                %                rankof_scaled(kc) = sum(ks>0);
                rankof_scaled(kc) = sum(ks);
                evof_scaled{kc} = eig(Abars((numstate-rankof_scaled(kc)+1):end,(numstate-rankof_scaled(kc)+1):end),'nobalance'); % eigenvalues of Ao
                evnof_scaled{kc} = eig(Abars(1:(numstate-rankof_scaled(kc)),1:(numstate-rankof_scaled(kc))),'nobalance'); % eigenvalues of Ano
            end
            eval(['save ' fname3 ' rank* svd* ev* B* C* umax Smat tend* -append'])
        end
    else % use eigs
        %    tic
        %    deigs=eigs(@(xpert) dirder_mod(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate,1);
        %    toc % can't use simple tic/toc since there are others embedded in code
        teigsstart = tic;
%        [veigs,deigs,flageigs]=eigs(@(xpert) dirder_mod(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate,1); % compute eigenvalues
        [veigs,deigs,flageigs]=eigs(@(xpert) dirder_normscale(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate,1); % compute eigenvalues
        toc(teigsstart)
        
        if ~useolfp
            eval(['save ' fname4 ' L bcl epsln ierr it_hist numpart sol x_hist deigs veigs flageigs useolfp maxabsclerr stimstart svd*'])
        else
            eval(['save ' fname4 ' L bcl epsln numpart sol deigs veigs flageigs useolfp maxabsclerr stimstart svd*'])
        end
        
    end
    
    %L = squeeze(allgains(:,:,i+1));
    %save kseparams L -append % these will be read in by karma_sim_estimtest_p2p.m
    if i < length(allstimstart)
        stimstart=allstimstart(i+1);
        save kseparams stimstart -append % these will be read in by karma_sim_estimtest_p2p.m
    end
    
end
if 0%~L
    svdctrb
    svdctrb_scaled
    svdobsv
    svdobsv_scaled
    clerr
    
    % for kk=1:numstate
    % maxeigncf(kk) = max(abs(evcf{kk}));
    % maxeignof(kk) = max(abs(evof{kk}));
    % mineigncf(kk) = min(abs(evcf{kk}));
    % mineignof(kk) = min(abs(evof{kk}));
    % end
    
    figure
    plot(1:numstate, flipud(sort(abs(eig(jaccd,'nobalance')))),'b-*')
    ylabel('eigenvalue norm')
    xlabel('eigenvalue index')
    
    figure 
    plot(real(eig(jaccd,'nobalance')),imag(eig(jaccd,'nobalance')),'*')
    axis equal
axis([-1.1 0.1 -.2 .2])
grid
xlabel('real(\lambda)','fontsize',20)
ylabel('imag(\lambda)','fontsize',20)
set(gca,'fontsize',20)

    figure
    plot(1:numstate,rankc)
    xlabel('state vector element')
    ylabel('controllability matrix rank')
    
    figure
    plot(1:numstate,ranko)
    xlabel('state vector element')
    ylabel('observability matrix rank')
    
    figure
    plot(1:numstate,rankcf)
    xlabel('state vector element')
    ylabel('controllability-f matrix rank')
    
    figure
    plot(1:numstate,rankof)
    xlabel('state vector element')
    ylabel('observability-f matrix rank')
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evcf{kk})
        p1 = plot(kk,abs(evcf{kk}),'b*');
        end
        if ~isempty(evncf{kk})
        p2 = plot(kk,abs(evncf{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (ctrbf)')
    xlabel('state vector element')
    legend([p1(1) p2(1)],'controllable','uncontrollable')
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evof{kk})
        p3 = plot(kk,abs(evof{kk}),'b*');
        end
        if ~isempty(evnof{kk})
        p4 = plot(kk,abs(evnof{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (obsvf)')
    xlabel('state vector element')
    legend([p3(1) p4(1)],'observable','unobservable')
    
    %%%%%%%%%%%%%%%%
    
    figure
    plot(1:numstate,rankc_scaled)
    xlabel('state vector element')
    ylabel('scaled controllability matrix rank')
    
    figure
    plot(1:numstate,ranko_scaled)
    xlabel('state vector element')
    ylabel('scaled observability matrix rank')
    
    figure
    plot(1:numstate,rankcf_scaled)
    xlabel('state vector element')
    ylabel('scaled controllability-f matrix rank')
    
    figure
    plot(1:numstate,rankof_scaled)
    xlabel('state vector element')
    ylabel('scaled observability-f matrix rank')
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evcf_scaled{kk})
        p1 = plot(kk,abs(evcf_scaled{kk}),'b*');
        end        
        if ~isempty(evncf_scaled{kk})
        p2 = plot(kk,abs(evncf_scaled{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (scaled ctrbf)')
    xlabel('state vector element')
    legend([p1(1) p2(1)],'controllable','uncontrollable')
    
    figure
    hold on;
    for kk=1:numstate
        if ~isempty(evof_scaled{kk})
        p3 = plot(kk,abs(evof_scaled{kk}),'b*');
        end
        if ~isempty(evnof_scaled{kk})
        p4 = plot(kk,abs(evnof_scaled{kk}),'r*');
        end
    end
    ylabel('eigenvalue norm (scaled obsvf)')
    xlabel('state vector element')
    legend([p3(1) p4(1)],'observable','unobservable')
end

% figure
% hold on;
% plot(1:numstate,maxeigncf)
% %plot(1:numstate,mineigncf)
% xlabel('state vector element')
% ylabel('controllability-f max eig norm')
% %axis([0 numstate 0 20])
%
% figure
% hold on;
% plot(1:numstate,maxeignof)
% %plot(1:numstate,mineignof)
% xlabel('state vector element')
% ylabel('observability-f matrix max eig norm')
% %axis([0 numstate 0 100])
%
% figure
% hold on;
% for kk=1:numstate
%     p1 = plot(kk,abs(evc{kk}),'b*');
%     p2 = plot(kk,abs(evnc{kk}),'r*');
% end
% ylabel('eigenvalue norm')
% xlabel('state vector element')
% legend([p1(1) p2(1)],'controllable','uncontrollable')
%
% figure
% hold on;
% for kk=1:numstate
%     p3 = plot(kk,abs(evo{kk}),'b*');
%     p4 = plot(kk,abs(evno{kk}),'r*');
% end
% ylabel('eigenvalue norm')
% xlabel('state vector element')
% legend([p3(1) p4(1)],'observable','unobservable')
%

% figure
% plot3(real(eig(jaccd,'nobalance')),imag(eig(jaccd,'nobalance')),zeros(1,numstate),'*')
% hold on
% for kk = 1:numstate
% plot3(real(eig(jaccd,'nobalance')),imag(eig(jaccd,'nobalance')),kk*ones(1,numstate),'*')
% plot3(real(evcf{kk}),imag(evcf{kk}),kk*ones(size(evcf{kk})),'r*')
% end