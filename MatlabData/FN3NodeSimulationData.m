%% OBSERVABILITY AND CONTROLLABILITY DATA SIMULATION, FITZHUGH-NAGUMO 3-NODE NETWORK DYNAMICS
% 
% From the Paper:
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
%
%THIS FILE WILL CREATE A DATA SIMULATION OF THE MOTIFS IN THE PAPER,
% VARIOUS COUPLED 3-NODE NETWORKS WITH FN NODAL DYNAMICS FOR THE FOLLOWING OPTIONS:
%   A) LIMIT CYCLE / CHAOTIC DYNAMIC REIGMES
%   B) EQUIVALANT COUPLING STRENGTHS / HETEROGENOUS COUPLING STRENGTHS
% OUTPUT IS A SERIES OF CSV FILES WITH THE TIME COURSE OF STATE VARIABLES
% NAMED ACCORDING TO THE SPECIFIED OPTIONS ABOVE
% ****ANDREW WHALEN 04/10/2014****
% Major input function revision to admit true sigmoid coupling
% email: awhalen@psu.edu
%%
function FN3NodeSimulationData
clear all;tic
close all
global f g dt nn
orient tall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTROL PARAMETERS:
sigma=0; % Initial Condition std dev.

FNinput=1; % 2 for Chaos, 5 for LC (1 for constant input LC)
IDcoupling=0; % 1 for ID, 0 for HET coupling
saveStates=1; % 1 to save data
doPlots=1; % 1 for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN SIMULATION CODE:
% Remove Transients:
stt=2400*5*19; edd=2400*5*20;

% Initial Conditions:
v0=sigma.*randn(3,1);
w0=sigma.*randn(3,1);

%%% COUPLING STRENGTH LOOP
% LOOPS THROUGH THE COUPLING STRENGTHS
l=1;
for k=1:20

    %Dimension of Grid:
    f=3*8; g=1;
    % Dimensions: (du for input vector, dq for param., dx augmented state, dy observation)
    du=1; dq=3; dx=dq+2*f*g+du; dy=f*g*2; %5 parameters plus 2 variables each cell

    randn('state',3);
    N=2400*5*20; % number of data samples
    dt=0.04;nn=1; % Time step

    % Initialize:
    x0=zeros(2*f*g,N);% true trajectory [v ; w]
    % a copy of the ICS for each motif:
    x0(:,1)=[repmat(v0,8,1);repmat(w0,8,1)]; %1 big column vector of all variables for initial time

    %%% FITZHUGH-NAGUMO DRIVING EXTERNAL INPUT
    % A WAVEFORM TO STIMULATE THE EXTERNAL FORCING ON THE SYSTEM
    switch(FNinput)
        % Constant Input:
        case 1
            I=-0.45*ones(1,N);%-1.15

            % Chaos 1,2,3:
        case 2
            I=(0.245/2).*[(square(((2*pi)/1.23)*dt.*(1:N),(.441/1.23)*100)+1)];
        case 3
            I=(0.198/2).*[(square(((2*pi)/1.43)*dt.*(1:N),(.573/1.43)*100)+1)];
        case 4
            I=(0.267/2).*[(square(((2*pi)/2.02)*dt.*(1:N),(.3/2.02)*100)+1)];

            % Standard limit cycle:
        case 5
            I=(0.5/2).*[(square(((2*pi)/5)*dt.*(1:N),(.3/5)*100)+1)];

            % Custom Input:
        case 6
            I=-0.8*[(sin(.1*(1:N).*dt)+1);
                (sin(.1*(1:N).*dt-(2*pi/3))+1);
                (sin(.1*(1:N).*dt-(4*pi/3))+1)]; %-1.410
        otherwise
            error('input selection failure')
    end

    I=repmat(I,3*8,1); % a copy of the inputs for each motif

    % Plot input:
    if doPlots==1
        splt=7734/dt;%1
        eplt=7754/dt;%N
        figure(10000)
        plot(dt.*[splt:eplt],I(1,splt:eplt),'r'); hold on;
        title('FN nodal input')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parameters:
    z=zeros(dq,3*g,N);
    s1=0.7;s2=0.8;s3=10.0; % parameter set
    z(:,1,:)=s1.*ones(3,N);
    z(:,2,:)=s2.*ones(3,N);
    z(:,3,:)=s3.*ones(3,N);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RuKu 4th order integrator:
    for n=1:N-1
        xx=x0(:,n);
        for i=1:nn
            k1=FN3NodeMotifs(xx,z(:,:,n),I(:,n),k,IDcoupling,FNinput); k1=dt*k1;
            k2=FN3NodeMotifs(xx+k1/2,z(:,:,n),I(:,n),k,IDcoupling,FNinput); k2=dt*k2;
            k3=FN3NodeMotifs(xx+k2/2,z(:,:,n),I(:,n),k,IDcoupling,FNinput); k3=dt*k3;
            k4=FN3NodeMotifs(xx+k3,z(:,:,n),I(:,n),k,IDcoupling,FNinput); k4=dt*k4;
            xx=xx+k1/6+k2/3+k3/3+k4/6;
        end;
        x0(:,n+1)=xx;
    end;

    % Create and save the observation data:
    States=[x0(:,stt:edd)];
    if saveStates == 1
        if((FNinput==2||FNinput==3||FNinput==4)&&(IDcoupling==1))
            svname=strcat('ChaosHetIC',num2str(l),'IdK3NodeFNCoupling',num2str(k));
        elseif((FNinput==5)&&(IDcoupling==1))
            svname=strcat('LCHetIC',num2str(l),'IdK3NodeFNCoupling',num2str(k));
        elseif((FNinput==1)&&(IDcoupling==1))
            svname=strcat('CHetIC',num2str(l),'IdK3NodeFNCoupling',num2str(k));
        elseif((FNinput==2||FNinput==3||FNinput==4)&&(IDcoupling==0))
            svname=strcat('ChaosHetIC',num2str(l),'HetK3NodeFNCoupling',num2str(k));
        elseif((FNinput==5)&&(IDcoupling==0))
            svname=strcat('LCHetIC',num2str(l),'HetK3NodeFNCoupling',num2str(k));
        elseif((FNinput==1)&&(IDcoupling==0))
            svname=strcat('CHetIC',num2str(l),'HetK3NodeFNCoupling',num2str(k));
        end

        save(svname,'States');
    end

    %% Plots Loop:
    if doPlots==1
        lbls={'motif 1','motif 1.5','motif 2','motif 3','motif 4','motif 5','motif 6','motif 7'};
        j=1;
        for i=1%1:3:24;  %1,4,7,10,13,16,19,22
            % voltage vars for each motif:
            %     figure(i)
            %     subplot(3,1,1)
            %     plot(dt.*(1:N),x0(i,1:N))
            %     subplot(3,1,2)
            %     plot(dt.*(1:N),x0(i+1,1:N))
            %     subplot(3,1,3)
            %     plot(dt.*(1:N),x0(i+2,1:N));
            splt=7734/dt;%1
            eplt=7754/dt;%N
            % voltage vars on the same plot (check for sync) for each motif:
            figure(2*k-1)%(i+1)
            plot(dt.*(splt:eplt),x0(i,splt:eplt),'go'); hold on;
            title(lbls(j))
            plot(dt.*(splt:eplt),x0(i+1,splt:eplt),'b');
            plot(dt.*(splt:eplt),x0(i+2,splt:eplt),'r--');
            %plot(dt.*(1:N),x0(i+21,1:N),'m');
            ylabel('v')

            % recovery vars on the same plot (check for sync) for each motif:
            figure(2*k)%(i+2)
            plot(dt.*(splt:eplt),x0(i+24,splt:eplt),'ko'); hold on;
            title(lbls(j))
            plot(dt.*(splt:eplt),x0(i+25,splt:eplt),'c');
            plot(dt.*(splt:eplt),x0(i+26,splt:eplt),'m--');
            ylabel('w')

            j = j+1;
        end
    end
   
end;


function [r]=FN3NodeMotifs(x,p,z,k,IDcoupling,FNinput)
global f g dt nn inputIDAvg inputHETAvg demeanInput
I = z;
v = x(1:f*g);
w = x(f*g+1:2*f*g);

% Sigmoid parameters (0 to 1 for input of 0 to 2):
% %%%%%     f=1/2*(tanh((x-1)/(2*1/10))+1);     %%%%%
Kappa=1;h=0;d=1/4;phi=1; %Kappa=1;h=0;d=1/4;phi=1;

% High to Low coupling strength
cspace=(linspace(0.004,0.9999,20));

% Identical Coupling strengths:
if(IDcoupling==1)
    k32=cspace(k);
    k12=k32;k21=k32;k13=k32;k31=k32;k23=k32;
end

% Symmetry breaking coupling strengths:
if(IDcoupling==0)
    rrr=[cspace(k)*0.10,-cspace(k)*0.10,cspace(k)*0.20,-cspace(k)*0.20,cspace(k)*0.30];
    k32=cspace(k);
    k12=k32+rrr(1);k21=k32+rrr(2);k13=k32+rrr(3);k31=k32+rrr(4);k23=k32+rrr(5);
end

% Network Connectivity Motifs:
% Motif 1:
input1m1=k21.*Kappa/2*(tanh((v(2)-h)/(2*d))+phi)+k31.*Kappa/2*(tanh((v(3)-h)/(2*d))+phi);
input2m1=k12.*Kappa/2*(tanh((v(1)-h)/(2*d))+phi)+k32.*Kappa/2*(tanh((v(3)-h)/(2*d))+phi);
input3m1=k13.*Kappa/2*(tanh((v(1)-h)/(2*d))+phi)+k23.*Kappa/2*(tanh((v(2)-h)/(2*d))+phi);
% Motif 1.5:
input1m15=k21.*Kappa/2*(tanh((v(5)-h)/(2*d))+phi);
input2m15=k12.*Kappa/2*(tanh((v(4)-h)/(2*d))+phi)+k32.*Kappa/2*(tanh((v(6)-h)/(2*d))+phi);
input3m15=k13.*Kappa/2*(tanh((v(4)-h)/(2*d))+phi)+k23.*Kappa/2*(tanh((v(5)-h)/(2*d))+phi);
% Motif 2:
input1m2=k21.*Kappa/2*(tanh((v(8)-h)/(2*d))+phi);
input2m2=k12.*Kappa/2*(tanh((v(7)-h)/(2*d))+phi)+k32.*Kappa/2*(tanh((v(9)-h)/(2*d))+phi);
input3m2=k23.*Kappa/2*(tanh((v(8)-h)/(2*d))+phi);
% Motif 3:
input1m3=k21.*Kappa/2*(tanh((v(11)-h)/(2*d))+phi)+k31.*Kappa/2*(tanh((v(12)-h)/(2*d))+phi);
input2m3=k12.*Kappa/2*(tanh((v(10)-h)/(2*d))+phi);
input3m3=k23.*Kappa/2*(tanh((v(11)-h)/(2*d))+phi);
% Motif 4:
input1m4=k21.*Kappa/2*(tanh((v(14)-h)/(2*d))+phi);
input2m4=k12.*Kappa/2*(tanh((v(13)-h)/(2*d))+phi)+k32.*Kappa/2*(tanh((v(15)-h)/(2*d))+phi);
input3m4=0;
% Motif 5:
input1m5=k21.*Kappa/2*(tanh((v(17)-h)/(2*d))+phi);
input2m5=k32.*Kappa/2*(tanh((v(18)-h)/(2*d))+phi);
input3m5=k23.*Kappa/2*(tanh((v(17)-h)/(2*d))+phi);
% Motif 6:
input1m6=k31.*Kappa/2*(tanh((v(21)-h)/(2*d))+phi);
input2m6=k12.*Kappa/2*(tanh((v(19)-h)/(2*d))+phi);
input3m6=k23.*Kappa/2*(tanh((v(20)-h)/(2*d))+phi);
% Motif 7:
input1m7=k21.*Kappa/2*(tanh((v(23)-h)/(2*d))+phi);
input2m7=k32.*Kappa/2*(tanh((v(24)-h)/(2*d))+phi);
input3m7=0;

input=[input1m1;input2m1;input3m1;input1m15;input2m15;input3m15;input1m2;input2m2;input3m2;input1m3;input2m3;input3m3;input1m4;input2m4;input3m4;input1m5;input2m5;input3m5;input1m6;input2m6;input3m6;input1m7;input2m7;input3m7];

% Fitzhugh-Nagumo ODE right hand side (vectorized):
p1=repmat(p(:,1),8,1);
p2=repmat(p(:,2),8,1);
p3=repmat(p(:,3),8,1);
%%%p1=p(:,1);p2=p(:,2);p3=p(:,3);
if FNinput == 1
    % Standard Form Eqns:
    vdot=p3.*(v-v.^3/3+w+I+input);
    wdot=-(v-p1+p2.*w)./p3;
else
    % Chaos form:
    vdot=p3.*(-w+v-(v.^3/3)+I+input);
    wdot=v-(p2.*w)+p1;
end
r=[vdot(:);wdot(:)]; % a column vector with both variables [v;w]
