function out = karma_altamp(Lpartial)
%To stop execution if Ctrl-C doesn't work, copied this from
%http://scrolls.mafgani.net/2007/03/matlab-interrupting-function-execution/
%To use it, create a file named stopexec in this function's home directory.
%To resume execution, type return at the prompt. 
statcode=exist('stopexec');
if statcode==2,
delete('stopexec');
%keyboard % doesn't work (MCC problem?) 
dbstop in karma_altamp 
end 
%load edata_1cell_b230_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_start1_cdiff_Kp0_Kpn0 sol 
%load edata_106cell_b230_d005_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-005_cdiff_Kp0_Kpn0 sol %jaccd
numstate1cell = 2;
%numpart = 1; 
numpart = 106;
numstate = numpart*numstate1cell;
bcl = 230; 
deltat=0.008;
stimstart = deltat; 
epsln = 1e-5;
%L = [Lpartial; 0]; % 1 cell
L = zeros(numstate,numpart); % observer feedback gain
%L(16,16) = Lpartial;
%L(:,16) = Lpartial;
L(16,16) = Lpartial(1);
%L(numpart+16,16) = Lpartial(2); 

save kseparams L numpart bcl stimstart

%desaa = 10; %ms
%aa = karma_sim_estimtest_p2p([-90*ones(numpart,1);0.8*ones(numpart,1)]); 
out = karma_sim_estimtest_p2p([-90*ones(numpart,1);0.8*ones(numpart,1)]); 
%out = desaa - abs(aa); 