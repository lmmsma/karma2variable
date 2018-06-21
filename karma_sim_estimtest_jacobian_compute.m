clear variables
%dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_ol' % This is a 1-BCL simulation generated from the fixed-point IC
%dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_1150_smooththeta_ol' 
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_920_smooththeta_ol' 

eval(['load ' dirname '/configinfo'])
capprox = 1; % Gain used for approximation of Heaviside function 1(V-Vn) as 0.5*(1+tanh(capprox*(V-Vn)))
bcli = stimperiod(1)/deltat; 
% hprime = zeros(numpart,numstep);
% acoeff1 = zeros(numpart,numstep);
% acoeff2 = zeros(numpart,numstep);
% acoeff3 = zeros(numpart,numstep);
hprime = zeros(numpart,writeintsteps);
acoeff1 = zeros(numpart,writeintsteps);
acoeff2 = zeros(numpart,writeintsteps);
acoeff3 = zeros(numpart,writeintsteps);

%aic = zeros((2*numpart)^2,writeintsteps); % coefficients of Aint(k) 
% 
Aint=eye(2*numpart); % cumulative Jacobian
ACLint=eye(2*numpart);
% Closed-loop values
B1 = [-deltat; 0];
K = [-0.8 0]; 
%eval(['load ' dirname '/1'])
for ii=1:numrep
        eval(['load ' dirname '/' num2str(ii)])
    for k = 1:(bcli)%1:writeintsteps%1:(stimperiod(1)/deltat)
        hprime(:,k) = ch1 + 2*ch2*V(:,k) + 3*ch3*V(:,k).^2; % Rappel version
        htemp = ch0 + ch1*V(:,k) + ch2*V(:,k).^2 + ch3*V(:,k).^3;

        acoeff1(:,k) = (-1+(Vstar - (n(:,k)/nB).^M).*hprime(:,k))/epsil2;
        acoeff2(:,k) = -(M*(1/nB)*(n(:,k)/nB).^(M-1)).*htemp/epsil2;
        acoeff3(:,k) = 0.5*(1/b)*(capprox*sech(capprox*(V(:,k)-Vn)).^2)/tauN;
        %            A(1,1) = Fom + 1 + acoeff1(1,k)*deltat - 2*Fom;
        %            A(1,2) = Fom;
        %            A(1,3) = acoeff2(1,k)*deltat;
        %            A(2,1) = Fom;
        %            A(2,2) = Fom + 1 + acoeff1(2,k)*deltat - 2*Fom;
        %             A(2,4) = acoeff2(2,k)*deltat;
        %             A(3,1) = acoeff3(1,k)*deltat;
        %             A(3,3) = 1-deltat/tauN;
        %             A(4,2) = acoeff3(2,k)*deltat;
        %             A(4,4) = 1-deltat/tauN;
        %             %         B1(1,1) = -cc_ext*deltat;
        %             %         B2(2,1) = -cc_ext*deltat;
        %             Ared=A(1:2,1:2);
        %             %         B1red=B1(1:2);
        %             %         B2red=B2(1:2);
        %             AIred(1:2,3:4) = Ared;
        %             AIred(3:4,3:4) = Ared;
        %             AI(1:4,5:8) = A;
        %             AI(5:8,5:8) = A;
        if numpart == 1 & ii==1 & k <= bcli
            Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,k)+1];
            A12 = [deltat*acoeff2(1,k)];
            Auc = [1-deltat/tauN];
            A21 = [deltat*acoeff3(1,k)];

            A = [Ac   A12
                A21  Auc];
            
            Aint = A*Aint;
 %           aic(:,k) = reshape(Aint,(2*numpart)^2,1); 
            if k >= 2 % feedback was turned on at k=2
                ACLint = (A-B1*K)*ACLint;
            end
        end
    end
%    eval(['save ' dirname '/' num2str(ii) ' acoeff* hprime -append'])
end

eig(Aint)
[va,da]=eig(Aint)
load edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_ol 
d - da
[vca,dca]=eig(ACLint)
load edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_Kp0p8
d - dca

% h=figure
% hold
% plot(bcli,Aint(1,1),'s')
% plot(bcli,reshape(Aint,4,1),'*')
% plot(1:bcli, aic(:,1:bcli))

% Compute analytical LTV jacobian coefficients 
adc = zeros((2*numpart)^2,bcli); % coefficients of Ak 
adc(1,:) = (-deltat*epsil1)/deltax^2+deltat*acoeff1(1,1:bcli)+1; 
adc(2,:) = deltat*acoeff2(1,1:bcli);
adc(3,:) = deltat*acoeff3(1,1:bcli);
adc(4,:) = (1-deltat/tauN)*ones(size(adc(4,:)));


%hh=figure
%hold
%plot(1:bcli,adc)

% Load instantaneous V-pert trajectories
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_jacVpert1e-7'
eval(['load ' dirname '/configinfo'])
eval(['load ' dirname '/' num2str(1)])
fi = strfind(dirname,'pert'); 
eval(['epsilon = ' dirname(fi+4:end)])
adce(1,:) = -[Verror(:,2:end) Verrornext]/epsilon; % The direction of subtraction has to be consistent: the perturbed elements have to both be the first element in the num/den, or both the second element in num/den
adce(2,:) = -[nerror(:,2:end) nerrornext]/epsilon;

% Load instantaneous n-pert trajectories
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_jacnpert1e-7'
eval(['load ' dirname '/configinfo'])
eval(['load ' dirname '/' num2str(1)])
fi = strfind(dirname,'pert'); 
eval(['epsilon = ' dirname(fi+4:end)])
adce(3,:) = -[Verror(:,2:end) Verrornext]/epsilon; % The direction of subtraction has to be consistent: the perturbed elements have to both be the first element in the num/den, or both the second element in num/den
adce(4,:) = -[nerror(:,2:end) nerrornext]/epsilon;
% Use different indexing for adce, for more convenient use with "reshape"
%There is a jump in the adc vs. adce error at the end of the run. Is there some problem with *errornext? 
figure
subplot(411)
plot(1:bcli,adc(1,:))
hold
plot(1:bcli,adce(1,:),'r:')
ylabel('a11')
legend('analytical', 'empirical')
subplot(412)
plot(1:bcli,adc(2,:))
hold
plot(1:bcli,adce(3,:),'r:')
ylabel('a12')
subplot(413)
plot(1:bcli,adc(3,:))
hold
plot(1:bcli,adce(2,:),'r:')
ylabel('a21')
subplot(414)
plot(1:bcli,adc(4,:))
hold
plot(1:bcli,adce(4,:),'r:')
ylabel('a22')
xlabel('time index')
% max(abs(adc-adce),[],2)
% 
%   1.0e-004 *
% 
%    0.08056742278773
%    0.14198369377993
%    0.00230926389122
%    0.00231310440912

% Compare empirical LTI Jacobian to product of instantaneous empirical LTV
% Jacobians
Aprodemp=eye(2*numpart); % cumulative Jacobian
for k=1:bcli
    Atemp = reshape(adce(:,k),2,2); 
    Aprodemp = Atemp*Aprodemp; 
end
load edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_eig_ol 
% jac =
% 
%   -0.00000019159391  -0.02698687574757
%   -0.00000714403838  -1.02943565932825
% Aprodemp =
%   -0.00000018732978  -0.02699565027856
%   -0.00000714320945  -1.02939094946443
  
if 0
newdirname = [dirname(1:end-9) 'compare' dirname(fi+4:end)]; 
eval(['save ' newdirname ' epsilon adce acoeff*'])
clear variables
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_jacnpert1e-7'
fi = strfind(dirname,'pert'); 
newdirname = [dirname(1:end-9) 'compare' dirname(fi+4:end)]; 
eval(['load ' dirname '/configinfo'])
eval(['save ' newdirname ' * -append'])
end
% hold
% plot(stimperiod(1)/deltat,Aint(1,1),'s')
% plot(stimperiod(1)/deltat,reshape(Aint,4,1),'*')
% plot(1:stimperiod(1)/deltat, aic(1:2,1:stimperiod(1)/deltat))


if 0
% These are older attempts to compute products of Jacobians empirically.
% Methods for obsv grammian terms worked, but those for ctrb grammians did not

% rename variables
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_920_smooththeta_Vpert1e-7' 
eval(['load ' dirname '/configinfo'])
Veq=V(:,1:(bcli));
neq=n(:,1:(bcli)); 
eval(['load ' dirname '/' num2str(1)])
Vvp=V(:,1:(bcli));
nvp=n(:,1:(bcli)); 
Vnextvp = Vnext; 
nnextvp = nnext; 
aice = zeros((2*numpart)^2,(bcli)); % coefficients of Aemp(k) 
% extract perturbation size from file name 
fi = strfind(dirname,'pert'); 
eval(['epsilon = ' dirname(fi+4:end)])
aice(1,:) = (Vvp-Veq)/epsilon; % The direction of subtraction has to be consistent: the perturbed elements have to both be the first element in the num/den, or both the second element in num/den
aice(2,:) = (nvp-neq)/epsilon;

% rename variables
dirname= 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_230_smooththeta_npert1e-7' 
eval(['load ' dirname '/configinfo'])
eval(['load ' dirname '/' num2str(1)])
aice(3,:) = (V-Veq)/epsilon; % The direction of subtraction has to be consistent: the perturbed elements have to both be the first element in the num/den, or both the second element in num/den
aice(4,:) = (n-neq)/epsilon;

figure(h)
plot(1:bcli, aice(:,1:bcli),':')

% Plot errors
figure 
plot(1:bcli, aice(:,1:bcli)-aic(:,1:bcli))

% Compute ctrb grammian quantities
    Aback=eye(2*numpart);
    abc = zeros((2*numpart)^2,writeintsteps); % coefficients of An-1An-2...Ak
    for ii=1:numrep
        %        eval(['load ' dirname '/' num2str(ii)])
        for k = (bcli):-1:1%1:writeintsteps%1:(stimperiod(1)/deltat)
            %         hprime(:,k) = ch1 + 2*ch2*V(:,k) + 3*ch3*V(:,k).^2; % Rappel version
            %         htemp = ch0 + ch1*V(:,k) + ch2*V(:,k).^2 + ch3*V(:,k).^3;
            %
            %         acoeff1(:,k) = (-1+(Vstar - (n(:,k)/nB).^M).*hprime(:,k))/epsil2;
            %         acoeff2(:,k) = -(M*(1/nB)*(n(:,k)/nB).^(M-1)).*htemp/epsil2;
            %         acoeff3(:,k) = 0.5*(1/b)*(capprox*sech(capprox*(V(:,k)-Vn)).^2)/tauN;
            if numpart == 1 & ii==1 & k <= bcli
                Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,k)+1];
                A12 = [deltat*acoeff2(1,k)];
                Auc = [1-deltat/tauN];
                A21 = [deltat*acoeff3(1,k)];

                A = [Ac   A12
                    A21  Auc];

                Aback = Aback*A;
                abc(:,k) = reshape(Aback,(2*numpart)^2,1);
            end
        end
        %    eval(['save ' dirname '/' num2str(ii) ' acoeff* hprime -append'])
    end
    hh=figure
    hold
    plot(bcli,Aback(1,1),'s')
    plot(bcli,reshape(Aback,4,1),'*')
    plot(1:bcli, fliplr(abc(:,1:bcli)))

    % None of the empirical methods for computing ctrb-grammian terms
    % listed below seemed to work 
    acce = zeros((2*numpart)^2,writeintsteps); % empirical coefficients of An-1An-2...Ak
    for k = 1:(bcli)
        acce(1:2,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[Vvp(:,bcli)-Veq(:,bcli);V(:,bcli)-Veq(:,bcli)]; %inv(...)*
        acce(3:4,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[nvp(:,bcli)-neq(:,bcli);n(:,bcli)-neq(:,bcli)]; %inv(...)*
        %    acce(1:2,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[Vnextvp-Veq(:,bcli);Vnext-Veq(:,bcli)]; %inv(...)*
        %    acce(3:4,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[nnextvp-neq(:,bcli);nnext-neq(:,bcli)]; %inv(...)*
    end
    % the 1,2 and 2,1 entries in acce are transposed relative to abc. Other
    % than that, the algorithm is only stable for a short span of timesteps counting
    % backwards from bcli.
    plot(1:stimperiod(1)/deltat, fliplr(acce(:,1:stimperiod(1)/deltat)),':')
end


% aie = zeros((2*numpart)^2,writeintsteps); % empirical coefficients of Ak
% for k = 1:(bcli)-1
%     aie(1:2,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[Vvp(:,k+1)-Veq(:,k+1);V(:,k+1)-Veq(:,k+1)]; %inv(...)*
%     aie(3:4,k) = [Vvp(:,k)-Veq(:,k) nvp(:,k)-neq(:,k); V(:,k)-Veq(:,k) n(:,k)-neq(:,k)]\[nvp(:,k+1)-neq(:,k+1);n(:,k+1)-neq(:,k+1)]; %inv(...)*
% end
%eval(['save ' dirname '/1 acoeff* hprime -append'])