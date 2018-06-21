clear variables 

filename = 'data_1cell_b800_pt14400p8_Vpert_0p03125'
eval(['load ' filename '/configinfo']) % load corresponding data

%perttime = 14400.8 
bcl=stimperiod(1)

h1=figure
subplot(4,1,1)
hold on;
ylabel('1,1')
title(['Elements of Aint (1-BCL mapping matrix) for BCL = ' num2str(bcl) ' ms'])
subplot(4,1,2)
hold on;
ylabel('1,2')
subplot(4,1,3)
hold on;
ylabel('2,1')
subplot(4,1,4)
hold on;
ylabel('2,2')
xlabel('time, ms')
h2=figure
subplot(2,1,1)
hold on;
title(['Eigenvalues of Aint (1-BCL mapping matrix) for BCL = ' num2str(bcl) ' ms'])
subplot(2,1,2)
hold on;
xlabel('time, ms')

x=floor(perttime/writeint)%9 % Indicate the interval before the one where the perturbation was applied
% The end time for this interval is x*writeint (e.g., 9*1600 = 14400)
eval(['load ' filename '/' num2str(x)]) % load corresponding data
st=perttime-x*writeint; % start time for perturbation, relative to data segment

allptis = round(st/deltat):3200:(bcl/deltat) % divide BCL into 5 segments
for jj=1:5
relpti = allptis(jj)
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
e=eig(Aint);
figure(h1)
subplot(4,1,1)
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(1,1),'x')
subplot(4,1,2)
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(1,2),'o')
subplot(4,1,3)
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(2,1),'+')
subplot(4,1,4)
plot(x*writeint+bcl+(relpti+1)*deltat,Aint(2,2),'>')
figure(h2) 
subplot(2,1,1)
plot(x*writeint+bcl+(relpti+1)*deltat,e(1),'x')
subplot(2,1,2)
plot(x*writeint+bcl+(relpti+1)*deltat,e(2),'o')
end


alleig = zeros(2,bcl/deltat); 
relpti=allptis(1)
K = [0.01 0];
ACLint=eye(2); % cumulative CL Jacobian
for k=[(1:(bcl/deltat))+relpti] %k=1:5 %% this is the Jacobian over indices 17 to 16016, for example
    Ac = [deltat*acoeff1(1,k)+1];
    A12 = [deltat*acoeff2(1,k)];
    Auc = [1-deltat/tauN];
    A21 = [deltat*acoeff3(1,k)];
    A = [ Ac   A12
        A21  Auc];
%    A=f;
    ACLint = (A-B1*K)*ACLint;
    alleig(:,k-relpti) = eig(ACLint); 
end

% Adjusting the gain from 0.01 to 0 to -0.01 appears to reveal that both
% cumulative eigenvalues, as a function of time, are shifted upward for negative gains and downward for
% positive gains. This leads to the alternans eigenvalue being reduced for
% negative gain, since it was negative to begin with. I'm not sure whether
% a stability criterion could be formulated based on the trend rather than
% the final values. Both Aint and ACLint have a zero eigenvalue (not sure
% why this is preserved in ACLint). 