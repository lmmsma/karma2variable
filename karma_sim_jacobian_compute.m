% clear variables
% 
% % Create Vd and nd for use with karma_sim_jacobian_compare_1cell 
% %filename = 'data_1cell_b800_18000_Vpert_0p03125'
% filename = 'data_2cell_b800_pt14400p8_Vpert1_0p03125'
% eval(['load ' filename '/configinfo']) % load corresponding data
% 
% x=floor(perttime/writeint)%9 % Indicate the interval before the one where the perturbation was applied
% % The end time for this interval is x*writeint (e.g., 9*1600 = 14400)
% eval(['load ' filename '/' num2str(x)]) % load corresponding data
% Vd =[V V];
% nd = [n n]; 
% %save steadystate_1cell_b800 Vd nd x 
% save steadystate_2cell_b800 Vd nd x 
% 
clear variables

filenames = {'data_2cell_b800_pt14400p8_Vpert1_0p03125','data_2cell_b800_pt14400p8_Vpert2_0p03125'}
%filename = 'data_1cell_b800_pt14560p8_Vpert_0p03125'
%filename = 'data'
eval(['load ' filenames{1} '/configinfo']) % load corresponding data

%perttime = 14400.8 
bcl=stimperiod(1)
x=floor(perttime/writeint)%9 % Indicate the interval before the one where the perturbation was applied
% The end time for this interval is x*writeint (e.g., 9*1600 = 14400)
st=perttime-x*writeint; % start time for perturbation, relative to data segment
relpti = round(st/deltat); % perturbation time index, relative to data segment

for jj=1:length(filenames)
    filename = filenames{jj}
eval(['load ' filename '/configinfo']) %
eval(['load ' filename '/' num2str(x)]) % load corresponding data

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
end

%filename2 = 'data_1cell_b800_pt14400p8_npert_1em4'
filenames2 = {'data_2cell_b800_pt14400p8_npert1_1em4','data_2cell_b800_pt14400p8_npert2_1em4'}

for jj=1:length(filenames2)
    filename2 = filenames2{jj}
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
end

%save Aemp_2cell_b800_sti0p85_Vpert_0p03125_npert_1em4 Aemp 
%save Aemp_1cell_b800_sti0p85_Vpert_0p03125_npert_1em4 Aemp Aint
%save Aemp_1cell_b800_sti160p85_Vpert_0p03125_npert_1em4 Aemp Aint