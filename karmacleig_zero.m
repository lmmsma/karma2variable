function out = karmacleig_zero(Lpartial)
%To stop execution if Ctrl-C doesn't work, copied this from
%http://scrolls.mafgani.net/2007/03/matlab-interrupting-function-execution/
%To use it, create a file named stopexec in this function's home directory.
%To resume execution, type return at the prompt. 
statcode=exist('stopexec');
if statcode==2,
delete('stopexec');
%keyboard % doesn't work (MCC problem?) 
dbstop in karma_cleigzero 
end 
%load edata_1cell_b230_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_start1_cdiff_Kp0_Kpn0 sol 
%load edata_106cell_b230_d005_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_cdiff_Kp0_Kpn0 sol %jaccd
load edata_106cell_b200_d005_dt008_tauV0p7_tauN170_Vstar4_eigs_relpert1e-005_start1_Kp0_Kpn0 sol
%eigopts.tol = 0.02/norm(jaccd); 
%eigopts.tol = 0.01/norm(jaccd); 
%eigopts.tol = 0.02/1.092296752166513; % 1 cell
numstate1cell = 2;
%numpart = 1; 
numpart = 106;
numstate = numpart*numstate1cell;
%bcl = 230; 
bcl = 200; 
deltat=0.008;
stimstart = deltat; 
epsln = 1e-5;
%L = [Lpartial; 0]; % 1 cell
L = zeros(numstate,numpart); % observer feedback gain
%L(16,16) = Lpartial;
%L(:,16) = Lpartial;
%L(16,16) = Lpartial(1);
%L(46,46) = Lpartial(1);
%L(numpart+16,16) = Lpartial(2); 
%L(46,16) = Lpartial(1);
%L(76,16) = Lpartial(1);
%L(15,16) = Lpartial(1);
%L(16,16) = Lpartial(1);
%L(17,16) = Lpartial(1);
%L(45,46) = Lpartial(1);
%L(46,46) = Lpartial(1);
%L(47,46) = Lpartial(1);
L((numpart+1):(2*numpart),16) = Lpartial(2); 
save kseparams L numpart bcl stimstart

%desmaxeig = .8; % 1 cell
%desmaxeig = .85; 
desmaxeig = -.85; 
%out = desmaxeig - abs(eigs(@(xpert) dirder_mod(sol,xpert,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln),numstate,1,'lm',eigopts));
% Afun, numstate = size of argument of Afun, 1 = number of eigenvalues to
% return, 'lm' = largest magnitude

% Since eigs is slow (maybe options structure is ignored in Afun mode), try power iteration instead
bprev = ones(numstate,1); % initialize vector

for i=1:2
    jactimesb=dirder_normscale(sol,bprev,'karma_sim_estimtest_p2p',feval('karma_sim_estimtest_p2p',sol),epsln); 
    % dirder_mod assumes that the input vector has one entry that is
    % exactly 1, with all the rest zero. This produces an error here, since b is initialized to ones. Restored something close to the
    % original logic in dirder_normscale. 
    bnext = jactimesb/norm(jactimesb); 
    mu = bprev'*jactimesb/(bprev'*bprev);
    bprev = bnext; 
end
%%out = desmaxeig - abs(mu); 
%out = desmaxeig - mu; % if using a zero-finder, such as nsoli
%out = mu; % if using a minimizer
out = mu^2; % if using fmincon