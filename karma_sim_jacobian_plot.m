% clear variables
% 
% % Create Vd and nd for use with karma_sim_jacobian_compare_1cell 
% filename = 'data_1cell_b800_18000_Vpert_0p03125'
% eval(['load ' filename '/configinfo']) % load corresponding data
% 
% x=floor(perttime/writeint)%9 % Indicate the interval before the one where the perturbation was applied
% % The end time for this interval is x*writeint (e.g., 9*1600 = 14400)
% eval(['load ' filename '/' num2str(x)]) % load corresponding data
% Vd =[V V];
% nd = [n n]; 
% save steadystate_1cell_b800 Vd nd x 

clear variables
h1=figure
h2=figure
h3=figure
h4=figure

filename = 'data_1cell_b800_pt14400p8_Vpert_0p03125'
%filename = 'data_1cell_b800_pt14560p8_Vpert_0p03125'
%filename = 'data'
%filename = 'data_2cell_b800_pt14400p8_Vpert1_0p03125'
eval(['load ' filename '/configinfo']) % load corresponding data

%perttime = 14400.8 
bcl=stimperiod(1)
x=floor(perttime/writeint)%9 % Indicate the interval before the one where the perturbation was applied
% The end time for this interval is x*writeint (e.g., 9*1600 = 14400)
eval(['load ' filename '/' num2str(x)]) % load corresponding data
st=perttime-x*writeint; % start time for perturbation, relative to data segment
relpti = round(st/deltat); % perturbation time index, relative to data segment
Aint=eye(2); % cumulative Jacobian
for k=[(1:(bcl/deltat))+relpti] % this is the Jacobian over indices 17 to 16016, for example
    Ac = [deltat*acoeff1(1,k)+1];
    A12 = [deltat*acoeff2(1,k)];
    Auc = [1-deltat/tauN];
    A21 = [deltat*acoeff3(1,k)];
    A = [ Ac   A12
        A21  Auc];
    Aint = A*Aint;
end
figure(h1)
hold on;
%plot(x*writeint+bcl+(relpti)*deltat,Aint(1,1),'p')
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(1,1),'p') 
% Not sure what time index to use: the cumulative Jacobian is supposed to
% be computed from 12800.85 up through 13600.80, and compared to empirical
% Jacobians with endpoints at 15200.85. 
figure(h2)
hold on;
%plot(x*writeint+bcl+(relpti)*deltat,Aint(2,1),'rs')
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(2,1),'rs')
figure(h3)
hold on;
%plot(x*writeint+bcl+(relpti)*deltat,Aint(1,2),'go')
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(1,2),'go')
figure(h4)
hold on;
%plot(x*writeint+bcl+(relpti)*deltat,Aint(2,2),'ch')
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(2,2),'ch')
%set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'

%'data_2cell_b800_pt14400p8_npert1_1em4'
%eval(['load data/' num2str(x)])
Vold=V;
nold=n;
eval(['load ' filename '/' num2str(x+1)]) % interval where perturbation was applied
%V(16015:16020)-Vold(16015:16020)
%ans(3)/Vpertval
%n(16015:16020)-nold(16015:16020)
%ans(3)/Vpertval
Verror_tp1 = V-Vold; 
nerror_tp1 = n-nold; 
Vpertind = find(Vpertval); 
npertind = find(npertval); 
if ~isempty(Vpertind)
for ii = 1:length(Vpertval) 
Aemp(ii,Vpertind)=Verror_tp1(ii,bcl/deltat + relpti+1)/Vpertval(Vpertind)
% This should be computing V at t= (14400 + 800 + 0.85) minus V at t = (14400 - 800 + 0.85), since writeint = 1600
% Changing the offset (+1) to +0 or +2 seems to slightly worsen agreement with Aint
Aemp(length(Vpertval)+ii,Vpertind)=nerror_tp1(ii,bcl/deltat + relpti+1)/Vpertval(Vpertind) 
end
end
figure(h1)
plot(x*writeint+bcl+(relpti+1)*deltat,Aemp(1,1),'*')
xlabel('1 BCL after perturbation time, ms')
ylabel('Jacobian element')
legend('Asim11', 'Aemp11')
figure(h2)
plot(x*writeint+bcl+(relpti+1)*deltat,Aemp(2,1),'r*')
xlabel('1 BCL after perturbation time, ms')
ylabel('Jacobian element')
legend('Asim21', 'Aemp21')

%filename2 = 'data_1cell_b800_pt14400p8_npert_1em4'
filename2 = 'data_1cell_b800_pt14560p8_npert_1em4'
eval(['load ' filename2 '/configinfo']) % load corresponding data
eval(['load ' filename2 '/' num2str(x+1)])

Verror_tp1 = V-Vold; 
nerror_tp1 = n-nold; 

Vpertind = find(Vpertval); 
npertind = find(npertval); 

if ~isempty(npertind)
for ii = 1:length(npertval) 
Aemp(ii,length(npertval)+npertind)=Verror_tp1(ii,bcl/deltat + relpti+1)/npertval(npertind)
Aemp(length(npertval)+ii,length(npertval)+npertind)=nerror_tp1(ii,bcl/deltat + relpti+1)/npertval(npertind)
end
end

figure(h3)
plot(x*writeint+bcl+(relpti+1)*deltat,Aemp(1,2),'g*')
xlabel('1 BCL after perturbation time, ms')
ylabel('Jacobian element')
legend('Asim12', 'Aemp12')
figure(h4)
plot(x*writeint+bcl+(relpti+1)*deltat,Aemp(2,2),'c*')
xlabel('1 BCL after perturbation time, ms')
ylabel('Jacobian element')
legend('Asim2', 'Aemp22')



%save Aemp_1cell_b800_sti0p85_Vpert_0p03125_npert_1em4 Aemp Aint
%save Aemp_1cell_b800_sti160p85_Vpert_0p03125_npert_1em4 Aemp Aint