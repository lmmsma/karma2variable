% Plot Karma model eigenvalues and ctrb/obsv plots resulting from
% karma_sim_estimtest_batch_eig.m using large-size font and symbols.
% Not sure if I already wrote a plotting function, but if I did, I can't find it.

clear variables
fs = 12;
ms = 10;

dir = 'H:\dell_optiplex_backup\laura\research\projects\alternans\matlab\macresults\';
%filename = 'edata_106cell_b230_d005_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-05_start1_cdiff_Kp0_Kpn0';
filename = 'edata_106cell_b200_d005_dt008_tauV0p7_tauN170_Vstar4_eig_relpert1e-05_start1_cdiff_Kp0_Kpn0';
load([dir filename])

numstate = size(jaccd,1);

figure
plot(1:numstate, flipud(sort(abs(eig(jaccd,'nobalance')))),'b-*','markersize',ms)
ylabel('eigenvalue norm')
xlabel('eigenvalue index')

figure
plot(real(eig(jaccd,'nobalance')),imag(eig(jaccd,'nobalance')),'*','markersize',ms)
axis equal
%axis([-1.1 0.1 -.2 .2])
axis([-1.3 0.5 -.2 .2])
grid
xlabel('real(\lambda)','fontsize',fs)
ylabel('imag(\lambda)','fontsize',fs)
set(gca,'fontsize',fs)
title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)

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
        p1 = plot(kk,abs(evcf{kk}),'b*','markersize',ms);
    end
    if ~isempty(evncf{kk})
        p2 = plot(kk,abs(evncf{kk}),'r*','markersize',ms);
    end
end
ylabel('eigenvalue norm (ctrbf)','fontsize',fs)
xlabel('state vector element','fontsize',fs)
legend([p1(1) p2(1)],'controllable','uncontrollable')
set(gca,'fontsize',fs)
title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)

figure
hold on;
for kk=1:numstate
    if ~isempty(evof{kk})
        p3 = plot(kk,abs(evof{kk}),'b*','markersize',ms);
    end
    if ~isempty(evnof{kk})
        p4 = plot(kk,abs(evnof{kk}),'r*','markersize',ms);
    end
end
ylabel('eigenvalue norm (obsvf)','fontsize',fs)
xlabel('state vector element','fontsize',fs)
legend([p3(1) p4(1)],'observable','unobservable')
set(gca,'fontsize',fs)
title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)

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
        p1 = plot(kk,abs(evcf_scaled{kk}),'b*','markersize',ms);
    end
    if ~isempty(evncf_scaled{kk})
        p2 = plot(kk,abs(evncf_scaled{kk}),'r*','markersize',ms);
    end
end
ylabel('eigenvalue norm (scaled ctrbf)','fontsize',fs)
xlabel('state vector element','fontsize',fs)
legend([p1(1) p2(1)],'controllable','uncontrollable')
set(gca,'fontsize',fs)
title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)

figure
hold on;
for kk=1:numstate
    if ~isempty(evof_scaled{kk})
        p3 = plot(kk,abs(evof_scaled{kk}),'b*','markersize',ms);
    end
    if ~isempty(evnof_scaled{kk})
        p4 = plot(kk,abs(evnof_scaled{kk}),'r*','markersize',ms);
    end
end`
ylabel('eigenvalue norm (scaled obsvf)','fontsize',fs)
xlabel('state vector element','fontsize',fs)
legend([p3(1) p4(1)],'observable','unobservable')
set(gca,'fontsize',fs)
title(['BCL = ' num2str(bcl) ' ms, start at ' num2str(stimstart) ' ms'],'fontsize',fs)
