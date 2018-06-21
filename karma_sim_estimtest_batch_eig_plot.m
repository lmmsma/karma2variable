% Plot results of gain tests
clear variables
fs = 24; % font size (use at least 32 pt for a full-screen figure that will occupy 1 column; roughly equal to 9 or 10pt font in paper)
% Set observer feedback, number of cells, and bcl
%allgains(:,:,1) = [0; -0.0005];
ratiox = 1.0; % relative placement for subfigure labels
ratioy = 1.02;

h1=figure
numsub = 1; %3 % number of subfigures
trajflag = 2; % use this flag to plot different gain ranges depending on the type of trajectory: 1=1to1, 2=2to1

for nn=1:numsub
    clear allgains
    
    if nn==1
        if trajflag == 1
            allgains(:,:,1) = [-1; 0];
            for i=2:11%12%2:13%2:2%2:6%
                allgains(:,:,i) = allgains(:,:,i-1)+[0.1;0];
            end
            ctr=i;
            allgains(:,:,ctr+1) = [-1.7; 0];
            for i=(ctr+2):(ctr+7)%2:13%2:2%2:6%
                allgains(:,:,i) = allgains(:,:,i-1)+[0.1;0];
            end
        elseif trajflag==2
            allgains(:,:,1) = [-0.35;0];
            for i=2:8%2:3
                allgains(:,:,i) = allgains(:,:,i-1)+[0.05;0];
            end
        end
        psflag = 0;
        
    elseif nn==2
        allgains(:,:,1) = [0; 0.001];
        for i=2:10%2:13%2:2%2:6%
            allgains(:,:,i) = allgains(:,:,i-1)+[0;0.001];
        end
        %        allgains(:,:,i+1) = [0; -0.0001];
        %        i=i+1;
        %allgains(:,:,i+1) = [-1.3878e-016; 0.000];
        %i=i+1;
        psflag = 0;
    elseif nn==3
        Kpgains = -0.9:0.1:-0.7;
        %Kpngains = 0.002:0.001:0.004;
        %Kpngains = 0.00005:0.00005:0.00015;
        %Kpngains = [0.0002 0.0005 0.0008];
        %        Kpngains = [0 0.0002 0.0005];
        Kpngains = [0 0.0002 0.0004];
        %Kpngains = [.007 .008 .009];
        tempgains = zeros(2*length(Kpgains),length(Kpngains));
        tempgains(1:2:(2*length(Kpgains)),:) = repmat(Kpgains,length(Kpgains),1);
        tempgains(2:2:(2*length(Kpgains)),:) = repmat(Kpngains,length(Kpngains),1)';
        tempgains=reshape(tempgains,2*length(Kpgains)*length(Kpgains),1);
        for i = 1:length(Kpgains)*length(Kpgains)
            allgains(:,:,i) = tempgains((2*i - 1):(2*i),:);
        end
        psflag = 1
    end
    
    if psflag
        psyms = {'*','+','o','*','+','o','*','+','o'};
    end
    
    allgains(:,:,i+1) = [0;0]; % dummy value (won't be used in loop below)
    
    %    allgains
    %    pause
    
    L = squeeze(allgains(:,:,1));
    numpart = 1;
    %    bcl = 230; %ms
    bcl = 160; %ms
    epsln = 1e-5; % relative perturbation size for use with diffjac_mod
    
    %h2=figure
    %hold on;
    for i=1:(size(allgains,3)-1)
        % Select file name
        % Replace decimal "p", otherwise it can cause problems loading files:
        for ii = 1:size(L,1)
            for jj = 1:size(L,2)
                temp = num2str(L(ii,jj)); % this can probably be vectorized but I'm not sure how (direct num2str will pad with spaces)
                kk = strfind(temp,'.'); % find decimal
                if kk
                    temp(kk) = 'p';
                end
                gainstr{ii}{jj} = temp;
            end
        end
        
        %    fname1 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        %    fname2 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(-epsln,'%1.e') '_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        %    fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        
        % need some modification for mac format:
        tempstr = num2str(epsln,'%1.e');
        if trajflag == 2
            epstr = tempstr;
        elseif length(tempstr) > 5
            epstr  = tempstr([1:3 5:6]);
        else
            epstr = tempstr
        end
        fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' epstr '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        if trajflag == 2
            fname3 = [fname3 '_2bcl'];
        end
        if ~exist([fname3 '.mat'])
            % use the Toshiba format
            fname3 = ['edata_' num2str(numpart) 'cell_b' num2str(bcl) '_dt008_tauV0p7_tauN170_Vstar4_eig_relpert' num2str(epsln,'%1.e') '_cdiff_Kp' gainstr{1}{1} '_Kpn' gainstr{2}{1}];
        end
        % This naming convention will only work for single cells. Not sure what to
        % do summarize multicell feedback
        
        eval(['load ' fname3])
        
        figure(h1)
        if nn == 1
            subplot(numsub,1,1)
            hold on;
            if psflag
                plot(L(1),norm(diag(d)),psyms{i})
            else
                plot(L(1),norm(diag(d)),'*')
            end
            
            % optional: use the two lines below to produce L21-variation plot
            %set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
            %legend('L_{21} = 0.0000','L_{21} = 0.0002','L_{21} = 0.0005')
            
            %figure(h2)
        elseif nn == 2
            subplot(numsub,1,2)
            hold on;
            if psflag
                plot(L(2),norm(diag(d)),psyms{i})
            else
                plot(L(2),norm(diag(d)),'*')
            end
            %    semilogy(L(2),norm(diag(d)),'*')
            
        elseif nn == 3 & numsub == 3
            subplot(3,1,3)
            hold on;
            if psflag
                plot(L(1),norm(diag(d)),psyms{i})
            else
                plot(L(1),norm(diag(d)),'*')
            end
        end
        
        L = squeeze(allgains(:,:,i+1));
        
    end
end

figure(h1)
set(gcf,'Position',[1 35 1280 694]) % maximize
subplot(numsub,1,1)
xl=xlabel('L_{11}');
%            yl=ylabel('2-norm of eigenvalue vector');
yl=ylabel('2-norm of e.v.');
%            ti = title('Closed-loop eigenvalue norm vs. gain');
set(gca,'FontSize',fs)
set(xl,'FontSize',fs)
% set(xl,'Position',[])
set(yl,'FontSize',fs)
%            set(ti,'FontSize',fs)
%            set(gca,'Position',[0.15 0.183333 0.314659 0.65]) % reduce height to make room for title
if trajflag == 1
    xmin = -1.9;
    xmax = 0.2;
    ymin = 0;
    ymax = 1.4;
    set(yl,'Position',[xmin-0.1*(xmax-xmin) (ymax-ymin)/2])
    set(xl,'Position',[-.851 -.2])
elseif trajflag == 2
    set(gcf,'Position',[192 483 491 495]) % smaller for ppt
    xmin = -.4;
    xmax = 0.05;
    ymin = 0;
    ymax = .8;
    set(gca,'ytick',[0.2 0.4 0.6 0.8])
end
axis([xmin xmax ymin ymax])
tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(a)')
set(tl,'FontSize',fs)
if numsub == 3
    set(gca,'Position',[0.182143 0.762275+0.005 0.722857 0.212365])
elseif trajflag == 1
    %            set(gca,'Position',[0.13 0.68188 0.775 0.24312]) % for small-size figure
    set(gca,'Position',[0.13 0.609498 0.775 0.315502])
else
    set(gca,'Position',[0.18 0.17 0.72 0.77])
end

subplot(numsub,1,2)
xl=xlabel('L_{21}');
yl=ylabel('2-norm of e.v.');
%            yl=ylabel('2-norm of eigenvalue vector');
%            ti=title('Closed-loop eigenvalue norm vs. gain');
set(gca,'FontSize',fs)
set(xl,'FontSize',fs)
set(yl,'FontSize',fs)
%           set(ti,'FontSize',fs)
%            set(gca,'Position',[0.589416 0.183333 0.315584 0.65]) % reduce height to make room for title
%            axis([-.0015 .011 10^-12 1])
%            tl=text(,'(b)')
xmin = 0;
xmax = .011;
ymin = 0;
ymax = 6e-10;
set(yl,'Position',[xmin-0.1*(xmax-xmin) (ymax-ymin)/2])
axis([xmin xmax ymin ymax])
%            set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
%            temp=log10(ymin)+ratioy*(log10(ymax)-log10(ymin)); % only use for log scale
%            tl=text(xmin + ratiox*(xmax-xmin),10^temp,'(b)') % only use for log scale
tl=text(xmin + ratiox*(xmax-xmin),ymin + ratioy*(ymax-ymin),'(b)')
set(tl,'FontSize',fs)
%            set(gca,'Position',[0.414638 0.183333 0.213406 0.741667])
if numsub == 3
    set(gca,'Position',[0.182143 0.442643 0.722857 0.212365])
else
    %            set(gca,'Position',[0.13 0.183333 0.775 0.24312])  % for small-size figure
    set(gca,'Position',[0.13 0.16 0.775 0.315502])  % for large-size figure
end

if numsub == 3
    subplot(3,1,3)
    xl=xlabel('L_{11}');
    yl=ylabel('2-norm of e.v.');
    %            yl=ylabel('2-norm of eigenvalue vector');
    %            ti=title('Closed-loop eigenvalue norm vs. gain');
    set(gca,'FontSize',fs)
    set(xl,'FontSize',fs)
    set(yl,'FontSize',fs)
    xmin = -.95;
    xmax = -.5;
    ymin = 10^-6;
    ymax = 1;
    temp=log10(ymin)+ratioy*(log10(ymax)-log10(ymin));
    set(yl,'Position',[xmin-0.1*(xmax-xmin) 10^(log10(ymin)+(log10(ymax)-log10(ymin))/2)])
    axis([xmin xmax ymin ymax])
    set(gca,'yscale','log') % semilogy does not appear to work with 'hold on'
    tl=text(xmin + ratiox*(xmax-xmin),10^temp,'(c)')
    set(tl,'FontSize',fs)
    
    %            set(gcf,'Position',[200 280 1000 480])
    set(gcf,'Position',[360 10 700 720])
    legend('L_{21} = 0','L_{21} = 0.0002','L_{21} = 0.0004')
    %            set(gca,'Position', [0.695435 0.183333 0.213406 0.741667])
    set(gca,'Position',[0.182143 0.11 0.722857 0.212365])
    %            set(gcf,'Position',[360 280 560 420])
end
