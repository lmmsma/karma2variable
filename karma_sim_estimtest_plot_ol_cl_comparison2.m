% Adapted from karma_sim_statefbtest_plot_ol_sfi_comparison.m, converted 
% to generate observer test plots. 
% Older comments: 
% Make a subfigure contour plot for OL vs. SFI
% Currently using A. Brun's "myaa" script to try to improve appearance of
% figures and avoid aliasing streaks that appear in many of the 
% pdfs. At present, only figure08a.eps and figure08b.eps have been
% processed with this script. There are still some streaks in the resulting
% figures, but perhaps fewer than before. I could try increasing the
% resolution in myaa, at the expense of increasing the figure sizes. 

clear variables
poster = 1; % use flag to make larger plots for poster
fs = 18; % poster fontsize
%fs = 24; % poster fontsize
%fs = 32; % poster fontsize
%lw1 = 3; % poster linewidth 1
%lw2 = 1.5; % poster linewidth 2
lw1 = 2; % poster linewidth 1
lw2 = 1.5; % poster linewidth 2
%dirname = 'data_100cell_b225_28000_ol';
%dirname = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_4200_smooththeta_ol';
dirname = 'edata_106cell_b200_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_2000_smooththeta_ol';
cmap=interp1([1 64], [0 0 0; 0.8 1 1], 1:64);
%clim2 = [-92 36] %***widen range based on V-vs-t plots
clim2 = [-100 50]; 
crange = [16 31 46 61 76 91]; %detlocindices
%ratiox = 1.02; % these are used to produce consistent placement for subfigure labels. 
%ratioy = 1.1;
ratiox = 1.02; % these are used to produce consistent placement for subfigure labels. 
ratioy = 1.05;

eval(['load ' dirname '/configinfo'])
%stimlocindex = 1


h=figure;
%set(gcf,'Position',[1 35 1280 694]) % maximize
set(gcf,'Position',[5 250 1212 475]) % This better approximates a width-filling 2-col figure

subplot(3,2,1)
hold on;
for ii=1:numrep % number of data writing cycles
    eval(['load ' dirname '/' num2str(ii)])
%%    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange 1],1:writeintsteps));
%    p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,iext(crange(3),1:50:writeintsteps));
%    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,iext(1,1:50:writeintsteps),'r--');
    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un(crange,1:writeintsteps));
    clear hprime acoeff* usf Verror nerror iext V n un
end
%%yl=ylabel('current (\muA/cm^2)');
%ti=title(['Stimulus and feedback terms vs. time']);
%yl=ylabel({'feedback', 'and pacing',  '(\muA/cm^2)'});
%yl=ylabel({'perturbation', 'to n',  '(no dim.)'});
yl=ylabel({'perturbation', 'to n'});
xmin = 0;
xmax = finaltime-100;
ymin = -1e-2;
ymax = 1e-2; 
axis([xmin xmax ymin ymax])
%set(gca,'yticklabelmode','manual')
%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(a)')
%axis([0 finaltime -1200 1200])
if poster
    set(gca,'FontSize',fs);
    set(p1,'LineWidth',lw1);
%    set(p2,'LineWidth',lw1);
    set(yl,'FontSize',fs);
%    set(tl,'FontSize',fs);
end
%%set(gca,'Position',[0.166406 0.709265 0.298253 0.215735])
%set(gca,'Position',[0.166406 0.739265 0.298253 0.215735])
set(gca,'Position',[0.166406-0.03 0.739265 0.33 0.215735])
colororder=get(gca,'ColorOrder');
% Define an offset colororder
colororderoffset = colororder([3:7 1:2],:);

subplot(3,2,3)
hold on;
clear p1 p2;
for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
    set(gca,'ColorOrder',colororder);
    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange(4),1:50:writeintsteps));
    set(gca,'ColorOrder',colororderoffset);
    p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Verror(crange(4),1:50:writeintsteps)+V(crange(4),1:50:writeintsteps),'--');
    clear hprime acoeff* usf Verror nerror iext V n un
end
xl=xlabel('time (ms)');
yl=ylabel({'membrane', 'potential', '(mV)'});
%set(xl,'Position',[2044.875 -134])
%ti=title(['V and V^d vs. time']);
%legend([p1(end); p2(end)],'V','V^d')
xmin = 0;
xmax = finaltime-100;
ymin = -100;
ymax = 40;
axis([xmin xmax ymin ymax])
%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(c)')
%axis([0 finaltime -100 40])
if poster
    set(gca,'FontSize',fs);
    set(p1,'LineWidth',lw1);
    set(p2,'LineWidth',lw2);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
 %   set(tl,'FontSize',fs);
end
% %set(gca,'Position',[0.166406 0.409632 0.298253 0.215735])
% set(gca,'Position',[0.166406 0.429632 0.298253 0.215735])
set(gca,'Position',[0.166406-0.03 0.429632 0.33 0.215735])

subplot(3,2,5)
hold on;
reprange = (numrep-5):(numrep-1);%10:20%1:10%11:20%41:50%
for ii=reprange % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
    tr=([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat;
    imagesc(deltax*(1:numpart), tr, V(:,1:10:writeintsteps)')
    clear hprime acoeff* usf Verror nerror iext V n un
end
axis xy
xmin = deltax;
xmax = fiberlength;
ymin = ((reprange(1)-1)*writeintsteps+1)*deltat;
ymax = tr(end);
axis([xmin xmax ymin ymax])
%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(e)')
%axis([deltax fiberlength ((reprange(1)-1)*writeintsteps+1)*deltat tr(end)])
colormap(cmap)
yl=ylabel('time (ms)');
xl=xlabel('distance (cm)');
set(xl,'Position',[1.067 670])
%ti=title('Membrane potential, V(x,t) (mV)');
if poster
    set(gca,'FontSize',fs);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
%    set(tl,'FontSize',fs);
end
caxis(clim2)
%colorbar
set(gca,'FontSize',fs);
% %set(gca,'Position',[0.166406 0.11 0.298253 0.215735])
% set(gca,'Position',[0.166406 0.105 0.298253 0.215735])
set(gca,'Position',[0.166406-0.03 0.105 0.33 0.215735])

% Closed loop
%%dirname = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_4200_smooththeta_Kp135_25';
%%dirname = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_4160_smooththeta_Kp1_25';
%dirname = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_4160_smooththeta_Kp1_75';
dirname = 'edata_106cell_b200_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_2000_smooththeta_KpVton_V135_toalln_5e-6_5e-6_1e-5_Vesat40';
eval(['load ' dirname '/configinfo'])

subplot(3,2,2)
hold on;
clear p1
for ii=1:numrep % number of data writing cycles
    eval(['load ' dirname '/' num2str(ii)])
%%    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange 1],1:writeintsteps));
%    p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,iext(crange(1),1:50:writeintsteps));
%    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,iext(1,1:50:writeintsteps),'r--');
   p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un(crange,1:writeintsteps));
    clear hprime acoeff* usf Verror nerror iext V n un
end
%yl=ylabel('current (\muA/cm^2)');
%ti=title(['Stimulus and feedback terms vs. time']);
%yl=ylabel({'observer', 'feedback', '(\muA/cm^2)'});
xmin = 0;
xmax = finaltime-100;
ymin = -1e-2;
ymax = 1e-2;
axis([xmin xmax ymin ymax])
%set(gca,'yticklabelmode','manual')

%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(b)')
%axis([0 finaltime -1200 1200])
if poster
    set(gca,'FontSize',fs);
    set(p1,'LineWidth',lw1);
%    set(p2,'LineWidth',lw1);
%    set(yl,'FontSize',fs);
%    set(tl,'FontSize',fs);
end
% %set(gca,'Position',[0.570341 0.709265 0.298253 0.215735])
% set(gca,'Position',[0.570341 0.739265 0.298253 0.215735])
set(gca,'Position',[0.570341-0.02 0.739265 0.33 0.215735])
legend('obsv. feedback','pacing pulse')

% colororder=get(gca,'ColorOrder');
% % Define an offset colororder
% colororderoffset = colororder([3:7 1:2],:);

subplot(3,2,4)
hold on;
clear p1 p2;
for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
    set(gca,'ColorOrder',colororder);
    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange(4),1:50:writeintsteps));
    set(gca,'ColorOrder',colororderoffset);
    p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Verror(crange(4),1:50:writeintsteps)+V(crange(4),1:50:writeintsteps),'--');
    clear hprime acoeff* usf Verror nerror iext V n un
end
xl=xlabel('time (ms)');
%set(xl,'Position',[2044.875 -134])
%yl=ylabel({'membrane', 'potential', '(mV)'});
%ti=title(['V and V^d vs. time']);
xmin = 0;
xmax = finaltime-100;
ymin = -100;
ymax = 40;
axis([xmin xmax ymin ymax])
%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(d)')
%axis([0 finaltime -100 40])
if poster
    set(gca,'FontSize',fs);
    set(p1,'LineWidth',lw1);
    set(p2,'LineWidth',lw2);
    set(xl,'FontSize',fs);
%    set(yl,'FontSize',fs);
%    set(tl,'FontSize',fs);
end
% %set(gca,'Position',[0.570341 0.409632 0.298253 0.215735])
% set(gca,'Position',[0.570341 0.429632 0.298253 0.215735])
set(gca,'Position',[0.570341-0.02 0.429632 0.33 0.215735])
legend([p1(end); p2(end)],'V','V meas.')


subplot(3,2,6)
hold on;
reprange = (numrep-5):(numrep-1);%10:20%1:10%11:20%41:50%
for ii=reprange % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
    tr=([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat;
    imagesc(deltax*(1:numpart), tr, V(:,1:10:writeintsteps)')
    clear hprime acoeff* usf Verror nerror iext V n un
end
axis xy
xmin = deltax;
xmax = fiberlength;
ymin = ((reprange(1)-1)*writeintsteps+1)*deltat;
ymax = tr(end);
axis([xmin xmax ymin ymax])
%tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(f)')
%axis([deltax fiberlength ((reprange(1)-1)*writeintsteps+1)*deltat tr(end)])
colormap(cmap)
%yl=ylabel('time (ms)');
xl=xlabel('distance (cm)');
set(xl,'Position',[1.067 670])
%ti=title('Membrane potential, V(x,t) (mV)');
if poster
    set(gca,'FontSize',fs);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
%    set(tl,'FontSize',fs);
end
caxis(clim2)
cb=colorbar;
set(cb,'Ytick',[-80 -40 0 40])
set(gca,'FontSize',fs);
%set(cb,'Position',[.90 .105 0.04 0.215735]) % I think this has to go before resizing gca otherwise it will automatically resize the axes
set(cb,'Position',[.91 .105 0.03 0.215735])
% %set(gca,'Position',[0.570341 0.11 0.298253 0.215735])
% set(gca,'Position',[0.570341 0.105 0.298253 0.215735])
set(gca,'Position',[0.570341-0.02 0.105 0.33 0.215735])

%subplot(1,2,1)
%colorbar off
%***
% Plot CL APDs at a distal location
    eval(['load ' dirname '\' num2str(numrep)])
    hfig = figure
    colororder=get(gca,'ColorOrder');
    close(hfig);
    nrco = size(colororder,1); % number of rows (distinct colors)

    psyms = repmat(strvcat('o','*','s','+','x','d','^','v','>','<','p','h'),2,1);
    colorcounter = 1;
    markerindex = 1;
    hwhat=figure
    for i=crange(4)
        hold on;
        if ~isempty(apds{i})
            if length(apendindices{i}) ~= length(apds{i})
                p(i) = plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
                % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                % stimulus is below peak AP height, (c) BCL is less than initial APD
            else
                p(i) = plot(deltat*apendindices{i}, deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
                %              axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
                %                axis([0 finaltime 50 200])%axis([0 10 .12 .14])%
                axis([0 finaltime 0 275])%axis([0 10 .12 .14])%
            end
            if i == detvec(1)
                %                ti=title('AP duration vs. time');
                ti=title('AP duration vs. repolarization time');
            end
        end
        colorcounter = colorcounter + 1;
        if colorcounter > nrco
            colorcounter = 1;
            markerindex = markerindex + 1;
        end
        %        xl=xlabel('repolarization time');
        xl=xlabel('time (ms)');
        yl=ylabel('APD (ms)');
    end
    %    legend(num2str(detvec),'Location','EastOutside')
    p(end+1)=plot(deltat*apendindices{i}, dapd*ones(size(apendindices{i})),'r--');
    %    legend(strvcat(num2str(detvec), 'desired apd'))
    if poster
        set(gca,'FontSize',fs);
        set(p,'LineWidth',lw1);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
        if numpart == 1
            legend off;
        end
%        set(gca,'Position',[0.13 0.11 0.775 0.312178])
    end

    if 0
dirname = 'edata_106cell_b230_dx02_dt008_d005_tauV0p7_tauN170_Vstar4_4200_smooththeta_ol';
eval(['load ' dirname '/configinfo'])
    figure
subplot(3,2,3)
hold on;
clear p1 p2;
for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
 %   set(gca,'ColorOrder',colororder);
%    p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange(4),1:50:writeintsteps));
%    set(gca,'ColorOrder',colororderoffset);
    p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Verror(crange(4),1:50:writeintsteps)+V(crange(4),1:50:writeintsteps),'--');
    clear hprime acoeff* usf Verror nerror iext V n un
end
xl=xlabel('time (ms)');
yl=ylabel({'membrane', 'potential', '(mV)'});
set(xl,'Position',[2044.875 -134])
xmin = 0;
xmax = finaltime-100;
ymin = -100;
ymax = 40;
axis([xmin xmax ymin ymax])
if poster
    set(gca,'FontSize',fs);
    set(p2,'LineWidth',lw2);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
end
set(gca,'Position',[0.166406-0.03 0.429632 0.33 0.215735])

set(gcf,'Position',[5 250 1212 475]) % This better approximates a width-filling 2-col figure
subplot(3,2,5)
hold on;
reprange = (numrep-5):(numrep-1);%10:20%1:10%11:20%41:50%
for ii=reprange % number of data writing cycles
    eval(['load ' dirname '\' num2str(ii)])
    tr=([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat;
    imagesc(deltax*(1:numpart), tr, Verror(:,1:10:writeintsteps)'+V(:,1:10:writeintsteps)')
    clear hprime acoeff* usf Verror nerror iext V n un
end
axis xy
xmin = deltax;
xmax = fiberlength;
ymin = ((reprange(1)-1)*writeintsteps+1)*deltat;
ymax = tr(end);
axis([xmin xmax ymin ymax])
colormap(cmap)
cb=colorbar;
set(cb,'Ytick',[-60 -20 20])
yl=ylabel('time (ms)');
xl=xlabel('distance (cm)');
if poster
    set(gca,'FontSize',fs);
    set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
    set(cb,'FontSize',fs);
end
set(xl,'Position',[1.067 2690])
set(gca,'FontSize',fs);
set(gca,'Position',[0.166406-0.03 0.105 0.33 0.215735])
    end