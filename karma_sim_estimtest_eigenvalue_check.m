% Repurpose part of karma_sim_estimtest_batch_eig to evaluate CL lead eigenvalue
% for a given gain setting. 

clear variables
numstate1cell = 2;
numpart = 106;
numstate = numpart*numstate1cell;
%bcl = 230; %bcl, ms
%bcl = 160; %ms
%bcl = 180; %ms
bcl = 200; %ms
L = zeros(numstate,numpart); % observer feedback gain
%load L_106cell_b230_des0p85_pwrit2 Lnsoli
%L(:,16) = Lnsoli; 

%load L_106cell_b230_des-0p85_pwrit2_V16n16 Lnsoli
%L(16,16) = Lnsoli(1);
%L(numpart+16,16) = Lnsoli(2); 

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
%  L(16,16) = -25;
%  L(46,16) = -25;
%  L(76,16) = -25;
%L(15,16) = -12.5;
%L(16,16) = -12.5;
%L(17,16) = -12.5;
%L(45,46) = -12.5;
%L(46,46) = -12.5; 
%L(47,46) = -12.5;
%L((numpart+1):(2*numpart),16) = 0.0001; 
%L((numpart+1):(2*numpart),46) = 0.0001; 
%L(1:numpart,16) = -3; 
%L(1:numpart,46) = -3; 
% 106-cell fiber detector locations: 16 31 46 61 76 91
% L(i,j) is feedback from measured V at cell j being applied to state i

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

%load edata_106cell_b230_d005_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_cdiff_Kp0_Kpn0 sol
load edata_106cell_b200_d005_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert1e-005_start1_Kp0_Kpn0 sol

teigsstart = tic;
[veigs,deigs,flageigs]=eigs(@(xpert) dirder_normscale(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate,1); % compute eigenvalues
toc(teigsstart)
deigs

