%% Read in and plot Obs and Ctrlb results for the paper
% Observability and Controllability of Nonlinear Networks: The Role of Symmetry
% Andrew J. Whalen, Sean N. Brennan, Timothy D. Sauer, and Steven J. Schiff
% Physical Review X, 2014
%
%We offer these codes so that others may replicate and build upon our findings. 
%Please cite the PRX paper from which these original results were published 
%if you make use of these codes, or modify them for your own data and research.
% 
% 11/18/2014
% Andrew Whalen
% Center for Neural Engineering
% Penn State University
% email: andrew.whalen@yale.edu (Updated email as of 1/1/2023)

% Ver 8.0, 03/04/2014
function PlotResults

clear all
close all

%% CONTROL PARAMETERS:
SaveMeans=0;
CalcMeans=1;
LoadData=1;
savePlots=0;
% Nm: motifs 1-8 , Nr: dynamics 1-Chaos 2-LC 3-Con
Nm=1;
Nr=1;


%% Load Obs/Ctrlb index values:
if ispc
path='MathematicaData\';
else
path='MathematicaData/';
end
fname1='FNCtrlbDELTAHet';
fname2='FNObsDELTAHet';
fname3='FNCtrlbDELTAId';
fname4='FNObsDELTAId';

% Initialize variables:
% # nodes, # coupling,  # data points
Nn=3;      Nc=20;       Np=12000;

if LoadData==1
    CIndexHET=ones(Nc,Nn,Np);
    OIndexHET=ones(Nc,Nn,Np);
    CIndexID=ones(Nc,Nn,Np);
    OIndexID=ones(Nc,Nn,Np);
    
    ii=[1 15 2 3 4 5 6 7];
    i=Nm; % Motif
    for j=1:Nc % Coupling Strengths
        for k=1:Nn % Node #
            
            fnameC=[path,fname1,'N',num2str(k),'M',num2str(ii(i)),'ChaosK',num2str(j),'IC1'];
            fnameO=[path,fname2,'N',num2str(k),'M',num2str(ii(i)),'ChaosK',num2str(j),'IC1'];
            fnameC2=[path,fname1,'N',num2str(k),'M',num2str(ii(i)),'LCK',num2str(j),'IC1'];
            fnameO2=[path,fname2,'N',num2str(k),'M',num2str(ii(i)),'LCK',num2str(j),'IC1'];
            fnameC3=[path,fname1,'N',num2str(k),'M',num2str(ii(i)),'ConsK',num2str(j),'IC1'];
            fnameO3=[path,fname2,'N',num2str(k),'M',num2str(ii(i)),'ConsK',num2str(j),'IC1'];
            
            if Nr==1
                CIndexHET(j,k,:)=load(fnameC,'-ascii');
                OIndexHET(j,k,:)=load(fnameO,'-ascii');
            elseif Nr==2
                CIndexHET(j,k,:)=load(fnameC2,'-ascii');
                OIndexHET(j,k,:)=load(fnameO2,'-ascii');
            elseif Nr==3
                CIndexHET(j,k,:)=load(fnameC3,'-ascii');
                OIndexHET(j,k,:)=load(fnameO3,'-ascii');
            end
            
            fnameC4=[path,fname3,'N',num2str(k),'M',num2str(ii(i)),'ChaosK',num2str(j),'IC1'];
            fnameO4=[path,fname4,'N',num2str(k),'M',num2str(ii(i)),'ChaosK',num2str(j),'IC1'];
            fnameC5=[path,fname3,'N',num2str(k),'M',num2str(ii(i)),'LCK',num2str(j),'IC1'];
            fnameO5=[path,fname4,'N',num2str(k),'M',num2str(ii(i)),'LCK',num2str(j),'IC1'];
            fnameC6=[path,fname3,'N',num2str(k),'M',num2str(ii(i)),'ConsK',num2str(j),'IC1'];
            fnameO6=[path,fname4,'N',num2str(k),'M',num2str(ii(i)),'ConsK',num2str(j),'IC1'];
            
            if Nr==1
                CIndexID(j,k,:)=load(fnameC4,'-ascii');
                OIndexID(j,k,:)=load(fnameO4,'-ascii');
            elseif Nr==2
                CIndexID(j,k,:)=load(fnameC5,'-ascii');
                OIndexID(j,k,:)=load(fnameO5,'-ascii');
            elseif Nr==3
                CIndexID(j,k,:)=load(fnameC6,'-ascii');
                OIndexID(j,k,:)=load(fnameO6,'-ascii');
            end
            
        end
    end
end


%% Calculate the mean and std:
if CalcMeans==1
    % Compute the mean and std dev for log-normal distribution:
    i=Nm; % Motifs
    for j=1:Nc % Coupling Strengths
        for k=1:Nn % Node #
            
            % Censor the zeros:
            ind1=squeeze(CIndexHET(j,k,:));
            ind2=squeeze(OIndexHET(j,k,:));
                        
            ind3=squeeze(CIndexID(j,k,:));
            ind4=squeeze(OIndexID(j,k,:));   
                        
            if ~isempty(ind1(ind1~=0))
                temp1 = ind1(ind1~=0);
            else
                temp1 = ones(size(ind1(ind1==0)));
            end
            if ~isempty(ind2(ind2~=0))
                temp2 = ind2(ind2~=0);
            else
                temp2 = ones(size(ind2(ind2==0)));
            end
            if ~isempty(ind3(ind3~=0))
                temp3 = ind3(ind3~=0);
            else
                temp3 = ones(size(ind3(ind3==0)));
            end
            if ~isempty(ind4(ind4~=0))
                temp4 = ind4(ind4~=0);
            else
                temp4 = ones(size(ind4(ind4==0)));
            end
            
     
            % Compute delta: proportion of zeros to non-zeros
            delta1=length(ind1(ind1==0))/length(ind1);
            delta2=length(ind2(ind2==0))/length(ind2);
            delta3=length(ind3(ind3==0))/length(ind3);
            delta4=length(ind4(ind4==0))/length(ind4);           
            
            % logn mean, upper and lower std calculation:
            m1 = mean(log(temp1)); s1=std(log(temp1)); u1=m1+s1; l1=m1-s1;
            m2 = mean(log(temp2)); s2=std(log(temp2)); u2=m2+s2; l2=m2-s2;
            m3 = mean(log(temp3)); s3=std(log(temp3)); u3=m3+s3; l3=m3-s3;
            m4 = mean(log(temp4)); s4=std(log(temp4)); u4=m4+s4; l4=m4-s4;
            
            %inverse transformation, compute mean, std upper and lower errorbars:
            
            muCH(j,k)=(1-delta1).*exp(m1+0.5*s1^2);
            muOH(j,k)=(1-delta2).*exp(m2+0.5*s2^2);
            muCI(j,k)=(1-delta3).*exp(m3+0.5*s3^2);
            muOI(j,k)=(1-delta4).*exp(m4+0.5*s4^2);

            uCH(j,k)=(1-delta1).*exp(u1+0.5*s1^2);
            uOH(j,k)=(1-delta2).*exp(u2+0.5*s2^2);
            uCI(j,k)=(1-delta3).*exp(u3+0.5*s3^2);
            uOI(j,k)=(1-delta4).*exp(u4+0.5*s4^2);
            
            lCH(j,k)=(1-delta1).*exp(l1+0.5*s1^2);
            lOH(j,k)=(1-delta2).*exp(l2+0.5*s2^2);
            lCI(j,k)=(1-delta3).*exp(l3+0.5*s3^2);
            lOI(j,k)=(1-delta4).*exp(l4+0.5*s4^2);            
        end
    end
end

if SaveMeans==1
    save('muOH','muOH');
    save('muCH','muCH');
    save('muOI','muOI');
    save('muCI','muCI');
   
    save('uOH','uOH');
    save('uCH','uCH');
    save('uOI','uOI');
    save('uCI','uCI');
    
    save('lOH','lOH');
    save('lCH','lCH');
    save('lOI','lOI');
    save('lCI','lCI');
end


%% Plots:
% Chaos Id and Het for Obs and Ctrlb:

% Coupling Strength:
k=linspace(0.004,0.9999,20);

% Correct out of range values on a log-scale plot:
muOH(muOH>=1)=1;
muCH(muCH>=1)=1;
uOH(uOH>=1)=1;
uCH(uCH>=1)=1;
lOH(lOH>=1)=1;
lCH(lCH>=1)=1;

muOH(muOH==0)=1e-14;
muCH(muCH==0)=1e-14;
uOH(uOH==0)=1e-14;
uCH(uCH==0)=1e-14;
lOH(lOH==0)=1e-14;
lCH(lCH==0)=1e-14;


%% PLOTS:

% Set up some constants for changing font sizes, etc
AxisFontsize=12;
TextFontsize=10;
FtNm='CMU Sans Serif';
ArrowWidth=10;
%Marker Sizes
ssz=6;  %sqaure
xsz=8; %x
dsz=15; %dot
bssz=3; %bounds square
bxsz=5; %bounds x
bdsz=8; %bounds dot

%parameters for figure and panel size
plotheight=4;
plotwidth=4;
subplotsx=2;
subplotsy=2;   
leftedge=0.45;
rightedge=0.04;   
topedge=.27;
bottomedge=.35;
spacex=0.1;
spacey=0.1;

% get the subplot coordinates
pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% Observability
fig5=figure(5);
hsp=subplot(2,2,1)
[ax,h1,h2]=plotyy(k,muOH(:,1),k,zeros(size(k)),'semilogy','plot'); hold(ax(1),'on'); hold(ax(2),'on');
if isempty(find(muOH(:,1),1)); hh1=ax(2); else hh1=ax(1); end;
plot(hh1,k,muOH(:,1),'g','Marker','^','LineWidth',1.7,'MarkerFaceColor','g','MarkerSize',ssz); hold on; plot(k,lOH(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');plot(k,uOH(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');
if isempty(find(muOH(:,2),1)); hh2=ax(2); else hh2=ax(1); end;
plot(hh2,k,muOH(:,2),'b','Marker','x','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',xsz); hold on; plot(k,lOH(:,2),'b:','Marker','x','MarkerSize',bxsz);plot(k,uOH(:,2),'b:','Marker','x','MarkerSize',bxsz);
if isempty(find(muOH(:,3),1)); hh3=ax(2); else hh3=ax(1); end;
plot(hh3,k,muOH(:,3),'r','Marker','.','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',dsz); hold on; plot(k,lOH(:,3),'r:','Marker','.','MarkerSize',bdsz);plot(k,uOH(:,3),'r:','Marker','.','MarkerSize',bdsz);
% Labels and Axes:
set(ax(1),'XColor','k','YColor','k'); set(ax(2),'XColor','k','YColor','k');
set(ax(1),'YLim',[1e-14,1e0]); set(ax(1),'YTick',[1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2]); YL=get(ax(2),'ylim'); set(ax(2),'ylim',[0 YL(2)]); set(ax(2),'YTick',[0,YL(2)]);
set(ax(1),'XLim',[-0.05,1.05]); set(ax(2),'XLim',[-0.05,1.05]); set(ax(1),'XTick',[0,0.5,1]);set(ax(2),'XTick',[0,0.5,1]);set(ax(2),'YTick',[]);
set(ax(1),'XTickLabel',''); set(ax(2),'XTickLabel',''); 
set(ax(1),'YTickLabel',{'0';'$10^{-12}$';'$10^{-10}$';'$10^{-8}$';'$10^{-6}$';'$10^{-4}$';'$10^{-2}$'},'FontSize',AxisFontsize); set(ax(2),'YTickLabelMode','auto','FontSize',AxisFontsize);
set(ax(1),'fontsize',AxisFontsize,'xgrid','off','ygrid','off')
set(ax(1),'Position',pos{3}); set(ax(2),'Position',pos{3});
% Line styles and Colors:
set(h1,'LineStyle','-','Color','g');
set(h2,'LineStyle','-','Color','k','LineWidth',0.1);
plot(ax(1),[-0.05 1.05],[1 1],'color','k','linewidth',1);
title(ax(1),'Observability');%,'FontName',FtNm); 
plotTickLatex2D;
text(0, 1e-1, 'a)','fontsize',AxisFontsize,'Fontname',FtNm);
hold(ax(1),'off'); hold(ax(2),'off');

%fig1=figure(6);
hsp=subplot(2,2,3)
[ax,h1,h2]=plotyy(k,muOI(:,1),k,zeros(size(k)),'semilogy','plot'); hold(ax(1),'on'); hold(ax(2),'on');
if isempty(find(muOI(:,1),1)); hh4=ax(2); else hh4=ax(1); end;
plot(hh4,k,muOI(:,1),'g','Marker','^','LineWidth',1.7,'MarkerFaceColor','g','MarkerSize',ssz); hold on; plot(k,lOI(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');plot(k,uOI(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');
if isempty(find(muOI(:,2),1)); hh5=ax(2); else hh5=ax(1); end;
plot(hh5,k,muOI(:,2),'b','Marker','x','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',xsz); hold on; plot(k,lOI(:,2),'b:','Marker','x','MarkerSize',bxsz);plot(k,uOI(:,2),'b:','Marker','x','MarkerSize',bxsz);
if isempty(find(muOI(:,3),1)); hh6=ax(2); else hh6=ax(1); end;
plot(hh6,k,muOI(:,3),'r','Marker','.','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',dsz); hold on; plot(k,lOI(:,3),'r:','Marker','.','MarkerSize',bdsz);plot(k,uOI(:,3),'r:','Marker','.','MarkerSize',bdsz);
% Labels and Axes:
set(ax(1),'XColor','k','YColor','k'); set(ax(2),'XColor','k','YColor','k');
set(ax(1),'YLim',[1e-14,1e0]); set(ax(1),'YTick',[1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2]); YL=get(ax(2),'ylim'); set(ax(2),'ylim',[0 YL(2)]); set(ax(2),'YTick',[0,YL(2)]);
set(ax(1),'XLim',[-0.05,1.05]); set(ax(2),'XLim',[-0.05,1.05]); set(ax(1),'XTick',[0,0.5,1]);set(ax(2),'XTick',[0,0.5,1]);set(ax(2),'YTick',[]);
set(ax(1),'XTickLabel','0|0.5|1'); set(ax(2),'XTickLabel','');
set(ax(1),'YTickLabel',{'0';'$10^{-12}$';'$10^{-10}$';'$10^{-8}$';'$10^{-6}$';'$10^{-4}$';'$10^{-2}$'},'FontSize',AxisFontsize); set(ax(2),'YTickLabelMode','auto','FontSize',AxisFontsize);
set(ax(1),'fontsize',AxisFontsize,'xgrid','off','ygrid','off')
set(ax(1),'Position',pos{1}); set(ax(2),'Position',pos{1});
% Line styles and Colors:
set(h1,'LineStyle','-','Color','g');
set(h2,'LineStyle','-','Color','k','LineWidth',1.0);
plot(ax(1),[-0.05 1.05],[1 1],'color','k','linewidth',1);
plotTickLatex2D('Axis',ax(1));
text(0, 1e-1, 'c)','fontsize',AxisFontsize,'Fontname',FtNm);
hold(ax(1),'off'); hold(ax(2),'off');

% Controllability
%fig1=figure(7);
hsp=subplot(2,2,2)
[ax,h1,h2]=plotyy(k,muCH(:,1),k,zeros(size(k)),'semilogy','plot'); hold(ax(1),'on'); hold(ax(2),'on');
if isempty(find(muCH(:,1),1)); hh7=ax(2); else hh7=ax(1); end;
plot(hh7,k,muCH(:,1),'g','Marker','^','LineWidth',1.7,'MarkerFaceColor','g','MarkerSize',ssz); hold on; plot(k,lCH(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');plot(k,uCH(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');
if isempty(find(muCH(:,2),1)); hh8=ax(2); else hh8=ax(1); end;
plot(hh8,k,muCH(:,2),'b','Marker','x','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',xsz); hold on; plot(k,lCH(:,2),'b:','Marker','x','MarkerSize',bxsz);plot(k,uCH(:,2),'b:','Marker','x','MarkerSize',bxsz);
if isempty(find(muCH(:,3),1)); hh9=ax(2); else hh9=ax(1); end;
plot(hh9,k,muCH(:,3),'r','Marker','.','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',dsz); hold on; plot(k,lCH(:,3),'r:','Marker','.','MarkerSize',bdsz);plot(k,uCH(:,3),'r:','Marker','.','MarkerSize',bdsz);
% Labels and Axes:
set(ax(1),'XColor','k','YColor','k'); set(ax(2),'XColor','k','YColor','k');
set(ax(1),'YLim',[1e-14,1e0]); set(ax(1),'YTick',[1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2]); YL=get(ax(2),'ylim'); set(ax(2),'ylim',[0 YL(2)]); set(ax(2),'YTick',[0,YL(2)]);
set(ax(1),'XLim',[-0.05,1.05]); set(ax(2),'XLim',[-0.05,1.05]); set(ax(1),'XTick',[0,0.5,1]);set(ax(2),'XTick',[0,0.5,1]);set(ax(2),'YTick',[]);
set(ax(1),'XTickLabel',''); set(ax(2),'XTickLabel','');
set(ax(1),'YTickLabel','','FontSize',AxisFontsize); set(ax(2),'YTickLabelMode','auto','FontSize',AxisFontsize);
set(ax(1),'fontsize',AxisFontsize,'xgrid','off','ygrid','off')
set(ax(1),'Position',pos{4}); set(ax(2),'Position',pos{4});
% Line styles and Colors:
set(h1,'LineStyle','-','Color','g');
set(h2,'LineStyle','-','Color','k','LineWidth',0.1);
plot(ax(1),[-0.05 1.05],[1 1],'color','k','linewidth',1);
title(ax(1),'Controllability');%,'FontName',FtNm); 
plotTickLatex2D;
text(0, 1e-1, 'b)','fontsize',AxisFontsize,'Fontname',FtNm);
hold(ax(1),'off'); hold(ax(2),'off');

%fig1=figure(8);
hsp=subplot(2,2,4)
[ax,h1,h2]=plotyy(k,muCI(:,1),k,zeros(size(k)),'semilogy','plot'); hold(ax(1),'on'); hold(ax(2),'on');
if isempty(find(muCI(:,1),1)); hh10=ax(2); else hh10=ax(1); end;
plot(hh10,k,muCI(:,1),'g','Marker','^','LineWidth',1.7,'MarkerFaceColor','g','MarkerSize',ssz); hold on; plot(k,lCI(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');plot(k,uCI(:,1),'g:','Marker','^','MarkerSize',bssz,'MarkerFaceColor','g');
if isempty(find(muCI(:,2),1)); hh11=ax(2); else hh11=ax(1); end;
plot(hh11,k,muCI(:,2),'b','Marker','x','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',xsz); hold on; plot(k,lCI(:,2),'b:','Marker','x','MarkerSize',bxsz);plot(k,uCI(:,2),'b:','Marker','x','MarkerSize',bxsz);
if isempty(find(muCI(:,3),1)); hh12=ax(2); else hh12=ax(1); end;
plot(hh12,k,muCI(:,3),'r','Marker','.','LineWidth',1.7,'MarkerFaceColor','none','MarkerSize',dsz); hold on; plot(k,lCI(:,3),'r:','Marker','.','MarkerSize',bdsz);plot(k,uCI(:,3),'r:','Marker','.','MarkerSize',bdsz);
% Labels and Axes:
set(ax(1),'XColor','k','YColor','k'); set(ax(2),'XColor','k','YColor','k');
set(ax(1),'YLim',[1e-14,1e0]); set(ax(1),'YTick',[1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2]); YL=get(ax(2),'ylim'); set(ax(2),'ylim',[0 YL(2)]); set(ax(2),'YTick',[0,YL(2)]);
set(ax(1),'XLim',[-0.05,1.05]); set(ax(2),'XLim',[-0.05,1.05]); set(ax(1),'XTick',[0,0.5,1]);set(ax(2),'XTick',[0,0.5,1]);set(ax(2),'YTick',[]);
set(ax(1),'XTickLabel','0|0.5|1'); set(ax(2),'XTickLabel','');
set(ax(1),'YTickLabel','','FontSize',AxisFontsize); set(ax(2),'YTickLabelMode','auto','FontSize',AxisFontsize);
set(ax(1),'fontsize',AxisFontsize,'xgrid','off','ygrid','off','fontweight','bold')
set(ax(1),'Position',pos{2}); set(ax(2),'Position',pos{2});
% Line styles and Colors:
set(h1,'LineStyle','-','Color','g');
set(h2,'LineStyle','-','Color','k','LineWidth',1.0);
plot(ax(1),[-0.05 1.05],[1 1],'color','k','linewidth',1);
plotTickLatex2D('Axis',ax(1));
text(0.93, 1e-1, 'd)','fontsize',AxisFontsize,'Fontname',FtNm);
text(-0.8, 1e-16, 'Network Connection Strength','fontsize',AxisFontsize);%,'Fontname',FtNm);
hold(ax(1),'off'); hold(ax(2),'off');



if savePlots == 1
    % Figure Save Size:
%     set(fig1,'PaperSize',[20 11]);
%     set(fig2,'PaperSize',[20 11]);
%     set(fig3,'PaperSize',[20 11]);
%     set(fig4,'PaperSize',[20 11]);
    set(fig5,'PaperSize',[plotwidth plotheight]);
%     set(fig6,'PaperSize',[20 11]);
%     set(fig7,'PaperSize',[20 11]);
%     set(fig8,'PaperSize',[20 11]);
%     set(fig9,'PaperSize',[20 11]);
%     set(fig10,'PaperSize',[20 11]);
%     set(fig11,'PaperSize',[20 11]);
%     set(fig12,'PaperSize',[20 11]);
    
%     set(fig1,'PaperPosition',[0 0 20 11]);
%     set(fig2,'PaperPosition',[0 0 20 11]);
%     set(fig3,'PaperPosition',[0 0 20 11]);
%     set(fig4,'PaperPosition',[0 0 20 11]);
    set(fig5,'PaperPosition',[0 0 plotwidth plotheight]);
%     set(fig6,'PaperPosition',[0 0 20 11]);
%     set(fig7,'PaperPosition',[0 0 20 11]);
%     set(fig8,'PaperPosition',[0 0 20 11]);
%     set(fig9,'PaperPosition',[0 0 20 11]);
%     set(fig10,'PaperPosition',[0 0 20 11]);
%     set(fig11,'PaperPosition',[0 0 20 11]);
%     set(fig12,'PaperPosition',[0 0 20 11]);
    
%     figure(1)
%     print -painters -dpdf -r600 ObsSymChaos.pdf
%     figure(2)
%     print -painters -dpdf -r600 ObsNonSymChaos.pdf
%     figure(3)
%     print -painters -dpdf -r600 ObsSymLC.pdf
%     figure(4)
%     print -painters -dpdf -r600 ObsNonSymLC.pdf
    figure(5)
    print -painters -dpdf -r600 M7Chaos.pdf
    %print '-PPDF Printer' M1Chaos.pdf
    %export_fig M1Chaos.pdf
    
%     figure(6)
%     print -painters -dpdf -r600 ObsNonSymCons.pdf
%     figure(7)
%     print -painters -dpdf -r600 CtrlbSymChaos.pdf
%     figure(8)
%     print -painters -dpdf -r600 CtrlbNonSymChaos.pdf
%     figure(9)
%     print -painters -dpdf -r600 CtrlbSymLC.pdf
%     figure(10)
%     print -painters -dpdf -r600 CtrlbNonSymLC.pdf
%     figure(11)
%     print -painters -dpdf -r600 CtrlbSymCons.pdf
%     figure(12)
%     print -painters -dpdf -r600 CtrlbNonSymCons.pdf
end
end

function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)

    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
 
    for i=1:nbx
       for j=1:nby
 
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           yfirst=bottommargin+(j-1.0)*(subysize+spacey);
 
           positions{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
 
       end
    end
end

function plotTickLatex2D(varargin)
% 
% Copyright (C) 2011 Alex Bikfalvi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%
% Optional arguments
optargin = size(varargin,2);

if mod(optargin,2) ~= 0
    error('The number of optional arguments must be even');
end

xLabelDy = 0;
yLabelDx = 0;
hAxis = gca;

i = 1;
while i <= optargin
    switch lower(varargin{i})
        case 'xlabeldy'
            xLabelDy = varargin{i+1};
        case 'ylabeldx'
            yLabelDx = varargin{i+1};
        case 'axis'
            hAxis = varargin{i+1};
        case 'fontsize'
            fontSize = varargin{i+1};
    end
    i = i + 2;
end

% Get properties for the specified axis
xLimit = get(hAxis,'XLim');
xTick = get(hAxis,'XTick');
xTickLabel = get(hAxis,'XTickLabel');
xTickLabelMode = get(hAxis,'XTickLabelMode');
xScale = get(hAxis,'XScale');
xAxisLocation = get(hAxis,'XAxisLocation');

yLimit = get(hAxis,'YLim');
yTick = get(hAxis,'YTick');
yTickLabel = get(hAxis,'YTickLabel');
yTickLabelMode = get(hAxis,'YTickLabelMode');
yScale = get(hAxis,'YScale');
yAxisLocation = get(hAxis,'YAxisLocation');

if ~exist('fontSize','var')
    fontSize = get(hAxis, 'FontSize');
end

% Get properties for the current axis
xLimitCurr = get(gca,'XLim');
xScaleCurr = get(gca,'XScale');

yLimitCurr = get(gca,'YLim');
yScaleCurr = get(gca,'YScale');

% Clear the current labels
set(hAxis,'XTickLabel',[]);
set(hAxis,'YTickLabel',[]);

% Get position of the specified axis
posAxis = get(hAxis,'Position');

% Get position of the current axis
posCurr = get(gca, 'Position');

% Convert x point between figure and axis data coordinates on linear scale
xFig2DatLinAxis = @(x)(xLimit(1) + diff(xLimit)*(x-posAxis(1))/posAxis(3));
xDat2FigLinAxis = @(x)(posAxis(1) + (x - xLimit(1))*posAxis(3)/diff(xLimit));

xFig2DatLinCurr = @(x)(xLimitCurr(1) + diff(xLimitCurr)*(x-posCurr(1))/posCurr(3));
xDat2FigLinCurr = @(x)(posCurr(1) + (x - xLimitCurr(1))*posCurr(3)/diff(xLimitCurr));

% Convert y point between figure and axis data coordinates on linear scale
yFig2DatLinAxis = @(y)(yLimit(1) + diff(yLimit)*(y-posAxis(2))/posAxis(4));
yDat2FigLinAxis = @(y)(posAxis(2) + (y - yLimit(1))*posAxis(4)/diff(yLimit));

yFig2DatLinCurr = @(y)(yLimitCurr(1) + diff(yLimitCurr)*(y-posCurr(2))/posCurr(4));
yDat2FigLinCurr = @(y)(posCurr(2) + (y - yLimitCurr(1))*posCurr(4)/diff(yLimitCurr));

% Convert x point between figure and axis data coordinates on logarithmic scale
xFig2DatLogAxis = @(x)(exp(log(xLimit(1)) + log(xLimit(2)/xLimit(1))*(x-posAxis(1))/posAxis(3)));
xDat2FigLogAxis = @(x)(posAxis(1) + posAxis(3)*log(x/xLimit(1))/log(xLimit(2)/xLimit(1)));

xFig2DatLogCurr = @(x)(exp(log(xLimitCurr(1)) + log(xLimitCurr(2)/xLimitCurr(1))*(x-posCurr(1))/posCurr(3)));
xDat2FigLogCurr = @(x)(posCurr(1) + posCurr(3)*log(x/xLimitCurr(1))/log(xLimitCurr(2)/xLimitCurr(1)));

% Convert y point between figure and axis data coordinates on logarithmic scale
yFig2DatLogAxis = @(y)(exp(log(yLimit(1)) + log(yLimit(2)/yLimit(1))*(y-posAxis(2))/posAxis(4)));
yDat2FigLogAxis = @(y)(posAxis(2) + posAxis(4)*log(y/yLimit(1))/log(yLimit(2)/yLimit(1)));

yFig2DatLogCurr = @(y)(exp(log(yLimitCurr(1)) + log(yLimitCurr(2)/yLimitCurr(1))*(y-posCurr(2))/posCurr(4)));
yDat2FigLogCurr = @(y)(posCurr(2) + posCurr(4)*log(y/yLimitCurr(1))/log(yLimitCurr(2)/yLimitCurr(1)));

% Convert x point between figure and axis [0,1] coordinates
xFig2Ax = @(x)((x - posAxis(1))/posAxis(3));
xAx2Fig = @(x)(posAxis(1) + x*posAxis(3));

% Convert y point between figure and axis [0,1] coordinates
yFig2Ax = @(y)((y - posAxis(2))/posAxis(4));
yAx2Fig = @(y)(posAxis(2) + x*posAxis(4));

switch xScale
    case 'linear'
        xFig2DatAxis = xFig2DatLinAxis;
        xDat2FigAxis = xDat2FigLinAxis;
    case 'log'
        xFig2DatAxis = xFig2DatLogAxis;
        xDat2FigAxis = xDat2FigLogAxis;
end

switch yScale
    case 'linear'
        yFig2DatAxis = yFig2DatLinAxis;
        yDat2FigAxis = yDat2FigLinAxis;
    case 'log'
        yFig2DatAxis = yFig2DatLogAxis;
        yDat2FigAxis = yDat2FigLogAxis;
end

switch xScaleCurr
    case 'linear'
        xFig2DatCurr = xFig2DatLinCurr;
        xDat2FigCurr = xDat2FigLinCurr;
    case 'log'
        xFig2DatCurr = xFig2DatLogCurr;
        xDat2FigCurr = xDat2FigLogCurr;
end

switch yScaleCurr
    case 'linear'
        yFig2DatCurr = yFig2DatLinCurr;
        yDat2FigCurr = yDat2FigLinCurr;
    case 'log'
        yFig2DatCurr = yFig2DatLogCurr;
        yDat2FigCurr = yDat2FigLogCurr;
end

if ~isempty(xTickLabel)
    % Set the X Axis
    
    xTickIndex = find((xTick >= xLimit(1)) & (xTick <= xLimit(2)));
    xTickVisible = xTick((xTick >= xLimit(1)) & (xTick <= xLimit(2)));
    xLabel = cell(length(xTickVisible),1);
    
    switch xTickLabelMode
        case 'auto'
            assert(length(xTickVisible) <= size(xTickLabel,1));
            switch xScale 
                case 'linear'
                    % Determine where there should be an exponent
                    xExp = max(abs(xLimit));
                    if xExp > 0
                        xExpLog = ceil(log10(xExp));
                        if(xExpLog > 0)
                            xExpLog = xExpLog - 1;
                        end
                        if (xExpLog > -3) && (xExpLog <= 3)
                            xExpLog = 0;
                        end
                    else
                        xExpLog = 0;
                    end
                    xExp = 10^xExpLog;
                    
                    for i=1:length(xTickVisible)
                        value = xTickVisible(i)/xExp;
                        if abs(value) <= eps
                            value = 0;
                        end
                        xLabel{i} = ['$' num2str(value) '$'];
                    end
                    
                    if (abs(xExpLog) > eps) && (abs(xExpLog) < 1/eps)
                        switch xAxisLocation
                            case 'bottom'
                                hText = text(...
                                    xFig2DatCurr(xDat2FigAxis(xLimit(2))),...
                                    yFig2DatCurr(yDat2FigAxis(yLimit(1))-0.06),...
                                    ['$\times 10^{' num2str(xExpLog) '}$'],...
                                    'HorizontalAlignment','Right',...
                                    'Interpreter','latex',...
                                    'FontSize', fontSize);
                            case 'top'
                                hText = text(...
                                    xFig2DatCurr(xDat2FigAxis(xLimit(2))),...
                                    yFig2DatCurr(yDat2FigAxis(yLimit(2))+0.06),...
                                    ['$\times 10^{' num2str(xExpLog) '}$'],...
                                    'HorizontalAlignment','Right',...
                                    'Interpreter','latex',...
                                    'FontSize', fontSize);
                        end
                        set(hText,'HitTest','off');
                    end                    
                case 'log'
                    for i=1:length(xTickVisible)
                        sgn = sign(xTickVisible(i));
                        xExp = log10(abs(xTickVisible(i)));
                        if abs(xExp) <= eps
                            xExp = 0;
                        end
                        xLabel{i} = ['$' num2str(10*sgn) '^{' num2str(xExp) '}$'];
                    end
            end
        case 'manual'
            if iscell(xTickLabel)
                for i=1:length(xTickVisible)
                    xLabel{i} = xTickLabel{1+mod(xTickIndex(i)-1,length(xTickLabel))};
                end
            else
                for i=1:length(xTickVisible)
                    xLabel{i} = xTickLabel(1+mod(xTickIndex(i)-1,size(xTickLabel,1)),:);
                end
            end
    end

    switch xAxisLocation
        case 'bottom'
            for i = 1:length(xTickVisible)
                hText = text(...
                    xFig2DatCurr(xDat2FigAxis(xTickVisible(i))),...
                    yFig2DatCurr(yDat2FigAxis(yLimit(1))-0.025),...
                    strtrim(xLabel{i}),...
                    'HorizontalAlignment','Center',...
                    'Interpreter','latex',...
                    'FontSize', fontSize);
                set(hText,'HitTest','off');
            end
        case 'top'
            for i = 1:length(xTickVisible)
                hText = text(...
                    xFig2DatCurr(xDat2FigAxis(xTickVisible(i))),...
                    yFig2DatCurr(yDat2FigAxis(yLimit(2))+0.025),...
                    strtrim(xLabel{i}),...
                    'HorizontalAlignment','Center',...
                    'Interpreter','latex',...
                    'FontSize', fontSize);
                set(hText,'HitTest','off');
            end
    end

    xLabel = get(hAxis,'XLabel');
    xLabelPos = get(xLabel,'Position');
    set(xLabel,'Position',[xLabelPos(1) yFig2DatCurr(yDat2FigAxis(yLimit(1))-0.07+xLabelDy) xLabelPos(3)]);

    xlim(hAxis, xLimit);
end

if ~isempty(yTickLabel)
    % Set the Y axis
    
    yTickIndex = find((yTick >= yLimit(1)) & (yTick <= yLimit(2)));
    yTickVisible = yTick((yTick >= yLimit(1)) & (yTick <= yLimit(2)));
    yLabel = cell(length(yTickVisible),1);
    
    switch yTickLabelMode
        case 'auto'
            assert(length(yTickVisible) <= size(yTickLabel,1));
            switch yScale
                case 'linear'
                    % Determine where there should be an exponent
                    yExp = max(abs(yLimit));
                    if yExp > 0
                        yExpLog = ceil(log10(yExp))-1;
                        if(yExpLog > 0)
                            yExpLog = yExpLog - 1;
                        end
                        if (yExpLog > -3) && (yExpLog <= 3)
                            yExpLog = 0;
                        end
                    else
                        yExpLog = 0;
                    end
                    yExp = 10^yExpLog;
                    
                    for i=1:length(yTickVisible)
                        value = yTickVisible(i)/yExp;
                        if abs(value) <= eps
                            value = 0;
                        end
                        yLabel{i} = ['$' num2str(value) '$'];
                    end
                    
                    if (abs(yExpLog) > eps) && (abs(yExpLog) < 1/eps)
                        switch yAxisLocation
                            case 'left'
                                hText = text(...
                                    xFig2DatCurr(xDat2FigAxis(xLimit(1))),...
                                    yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03),...
                                    ['$\times 10^{' num2str(yExpLog) '}$'],...
                                    'Interpreter','latex',...
                                    'FontSize', fontSize);
                            case 'right'
                                hText = text(...
                                    xFig2DatCurr(xDat2FigAxis(xLimit(2))),...
                                    yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03),...
                                    ['$\times 10^{' num2str(yExpLog) '}$'],...
                                    'HorizontalAlignment','Right',...
                                    'Interpreter','latex',...
                                    'FontSize', fontSize);
                        end
                        set(hText,'HitTest','off');
                    end                    
                case 'log'
                    for i=1:length(yTickVisible)
                        sgn = sign(yTickVisible(i));
                        yExp = log10(yTickVisible(i));
                        if abs(yExp) <= eps
                            yExp = 0;
                        end
                        yLabel{i} = ['$' num2str(10*sgn) '^{' num2str(yExp) '}$'];
                    end
            end
        case 'manual'
            if iscell(yTickLabel)
                for i=1:length(yTickVisible)
                    yLabel{i} = yTickLabel{1+mod(yTickIndex(i)-1,length(yTickLabel))};
                end
            else
                for i=1:length(yTickVisible)
                    yLabel{i} = yTickLabel(1+mod(yTickIndex(i)-1,size(yTickLabel,1)),:);
                end
            end
    end

    xMin = bitmax;
    
    switch yAxisLocation
        case 'left'
            for i = 1:length(yTickVisible)
                hText = text(...
                    xFig2DatCurr(xDat2FigAxis(xLimit(1))-0.01),...
                    yFig2DatCurr(yDat2FigAxis(yTickVisible(i))),...
                    strtrim(yLabel{i}),...
                    'HorizontalAlignment','Right',...
                    'Interpreter','latex',...
                    'FontSize', fontSize);
                set(hText,'HitTest','off');
                set(hText,'Units','normalized');
                xExt = get(hText,'Extent');
                set(hText,'Units','data');
                xMin = min(xMin, xExt(1));
            end
        case 'right'
            for i = 1:length(yTickVisible)
                hText = text(...
                    xFig2DatCurr(xDat2FigAxis(xLimit(2))+0.01),...
                    yFig2DatCurr(yDat2FigAxis(yTickVisible(i))),...
                    strtrim(yLabel{i}),...
                    'HorizontalAlignment','Left',...
                    'Interpreter','latex',...
                    'FontSize', fontSize);
                set(hText,'HitTest','off');
                set(hText,'Units','normalized');
                xExt = get(hText,'Extent');
                set(hText,'Units','data');
                xMin = min(xMin, xExt(1));
            end
    end

    yLabel = get(hAxis,'YLabel');
    yLabelPos = get(yLabel,'Position');
    set(yLabel,'Position',[xFig2DatCurr(xAx2Fig(xMin) - 0.002 + yLabelDx) yLabelPos(2) yLabelPos(3)]);

    ylim(hAxis, yLimit);
end
    

end

