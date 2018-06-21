%Plots MSVs of obsv Grammians (with a/b labels) for CinC paper. 

clear variables

poster = 1; % use flag to make larger plots for poster
fs = 32; % poster fontsize (use at least 32 pt for a full-screen figure that will occupy 1 column; roughly equal to 9 or 10pt font in paper) 
lw1 = 3; % poster linewidth 1
lw2 = 3; % poster linewidth 2
dirname = 'edata_1cell_b230_dt008_d005_tauV0p7_tauN170_Vstar4_920_smooththeta_ol'; % Directory for data file
eval(['load ' dirname '/configinfo'])
numpd = 3;
%eval(['load ' dirname '/svdslong_allscale_u1_4_range450_60pts_3pd'])
eval(['load ' dirname '/svds_allscale_window2875_20pts_2pd'])

%init=round(900.05/deltat); % start time of a selected AP
fileindex = floor(init/writeintsteps)+1
initindex =  init - (fileindex-1)*writeintsteps; %
init = initindex;
eval(['load ' dirname '/' num2str(fileindex)])

minV=min(V(:,init:endstep));
trueoffset = (fileindex-1)*writeintsteps;

%%%%%%%%%%%
% This version plots 2 figures on one plot
h1=figure
%set(gcf,'Position',[2 100 1280 400])
subplot(2,1,1)
hold on;
p1=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogs_plot(2,:),'-b+');
p2=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogns_plot(2,:),'-gs');
set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
%axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^7])
axis([(trueoffset+init)*deltat/stimperiod(1) 1.2 10^-15 10^12])
%yl=ylabel('window = 0.1\timesBCL');
yl=ylabel('min. sing. value');
%xl=xlabel('start time for Grammian, in BCL indices');
% %    axis([(trueoffset+init)*deltat (trueoffset+windowlength)*deltat 10^-8 10^2])
%     axis([(trueoffset+init)*deltat (trueoffset+endstep)*deltat 10^-8 10^2])
%ti=title('Minimum singular values of controllability Grammians');
legend('Case 1: observe V directly','Case 2: observe n directly')

if poster
    set(p1,'LineWidth',lw2);
    set(p2,'LineWidth',lw2);
%    set(ti,'FontSize',fs);
    set(gca,'FontSize',fs);
    set(yl,'FontSize',fs);
%    set(xl,'FontSize',fs);
end
set(gca,'ytick',[10^-10 1 10^10])
%tl=text(3.85,10^8,'(a)')
tl=text(1.22,10^12,'(a)')
set(tl,'FontSize',fs)

%eval(['load ' dirname '/svdslong_allscale_u1_4_range4500_60pts_3pd'])
eval(['load ' dirname '/svds_allscale_window28750_20pts_1pt2pd'])

%h2=figure
%set(gcf,'Position',[2 100 1280 400])
subplot(2,1,2)
hold on;
p1=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogs_plot(2,:),'-b+');
p2=plot([trueoffset + init+((1:numdiv)-1)*round(increment/deltat)]*deltat/stimperiod(1),svdogns_plot(2,:),'-gs');
set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2+7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^3,'--');
p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 10)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^12,':');
%axis([(trueoffset+init)*deltat/stimperiod(1) (trueoffset+init+numpd*stimperiod(1)/deltat)*deltat/stimperiod(1) 10^-20 10^7])
axis([(trueoffset+init)*deltat/stimperiod(1) 1.2 10^-16 10^12])
%yl=ylabel('window = 1\timesBCL');
yl=ylabel('min. sing. value');
%xl=xlabel('start time for Grammian');
xl=xlabel('start time for Grammian, in BCL indices');
if poster
    set(p1,'LineWidth',lw2);
    set(p2,'LineWidth',lw2);
 %   set(ti,'FontSize',fs);
    set(gca,'FontSize',fs);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
    set(p3,'LineWidth',lw2); 
    set(p4,'LineWidth',lw2); 
end
set(gca,'ytick',[10^-10 1 10^10])
%tl=text(3.85,10^8,'(b)')
tl=text(1.22,10^12,'(b)')
set(tl,'FontSize',fs)

subplot(2,1,1)
%figure(h1)
%p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 - 7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^7,'--');
%p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 3)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^30,'--');
p3=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2+7)*(V(:,init:(init+numpd*stimperiod(1)/deltat))-minV+1).^3,'--');
p4=plot((trueoffset + (init:(init+numpd*stimperiod(1)/deltat)))*deltat/stimperiod(1),10^(log10(ylim(1))+(log10(ylim(2))-log10(ylim(1)))/2 + 10)*(n(:,init:(init+numpd*stimperiod(1)/deltat))).^12,':');
set(p3,'LineWidth',lw2); 
set(p4,'LineWidth',lw2); 
