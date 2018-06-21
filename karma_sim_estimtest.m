% This is a restart of the single-cell estimator code, based off the latest
% version of the statefbtest.m code (from version 17).
% Changes: increased writeint to be large enough to fit 1 bcl at 1000ms,
% loaded processed V data in place of Vd, changed Vd handling in several
% locations, including 1) at the beginning of the for-loop, to read in
% successive pieces of V data, and 2) at the plotting function, so that the
% same piece is not plotted over and over, introduced sensorloc array to
% indicate positions of sensing electrodes.
% To increase conduction speed, increased epsil1. Had to increase current
% energy to compensate for the increase in diffusion (APs have trouble
% reaching threshold if this is not done)

% Older comments: 
% Changed varargin list so that the arguments are (bcl, K, length). The latter
% two are optional. This ordering of arguments is not the same as karma_sim_pyra.m. 
% Use K=0 to indicate no feedback, independently of size of system.
% Should change to include some kind of indexing system for
% desfile names. 
% Added 'desfile' name to list of saved quantities in configinfo. 
% gains P 0, -.02 -.04, I 0, -0.0004, -0.0008
% done: 0,0, 0,4; 2,0, 2,4; 4,4
% todo: 0,8; 2,8; 4,0, 4,8
% Added gain as a parameter
% To convert back to fiber geometry, increase number of cells, reduce
% writeint, increase deltat (for rough testing only), comment out acoeff,
% hprime, and htemp calculations, change epsil1 back to 0.001 from 0,
% increase stimduration from 2.5 to 7, change desfile, comment out A, B1,
% B2, and change 2 save commands to avoid saving acoeff, hprime, B1, B2. In
% plots, restore crange, and subscripts on V and n in legends. 

% YET ANOTHER ERROR in acoeff's: forgot the factor of 1/nB in acoeff2.
% Since nB=1 in Rappel's units, it shouldn't change the answer, but it
% affected previous versions in the old unit system. 

% Older comments: 
% This is a single-cell configuration of the code, used to test state feedback
% (based on either V or n measurements) on either the V or n dynamic
% equations. Comments on gains that work are included below. The default
% setting should be open-loop. 

% Main changes: 1) revised code to use wraparound data structure to avoid memory
% errors, 2) corrected mistakes in Jacobian code. Some results are
% documented in jobtalk.ppt. 

% Warning: I modified acoeff* and hprime to use V, n instead of Vd, nd, to
% simplify the Grammian tests. This should be changed back for later
% versions. 

% Validation of LTV: 325ms, set fbon_time to 1, stimduration to 2.5, use
% 1cellOL_...325_3250.mat. Comment in s, A-related quantities.
% Perturb V and n at 670 and reset s to match. Perturb KIred1 5 steps later.
% Restore acoeff, hprime, and htemp dependence on Vd and nd. 
% Plot V, nonlinear errors, and LTV states. The LTV appears to match the
% nonlinear model initially. The LTV test configuration is saved in
% karma_sim_statefbtest_ltvcheck.m 

% Changed h function to use ch* parameters. This was an omission in previous
% versions that led to the wrong hprime. Fixed acoeff2 definition.

%%%%%%%% Errors upon errors in Rappel version. 
% There is an incomplete conversion to the new "h" that was used from
% initial change (for UNYCES) until 01/06/09. The "ch1,2,3" parameters were initially
% intended for use in the modified tanh approximation to the power-series
% expansion in the Rappel version. I ended up NOT USING the modified tanh;
% my intent eventually was just to use Rappel et al.'s power series for h.
% However, there are at least 2 separate errors in the originally modified
% version of acoeffX. 1) a power series is used for hprime in acoeff1, but
% the tanh form persists in acoeff2. It should be one or the other (both
% power series or both tanh). 2) the chX coefficients were never used in
% the actual computation of the power series version of h, yet they are
% introduced in hprime. These coefficients are nowhere near the correct
% values for the power series, so it makes no sense to insert them in a
% power series derivative for hprime. Reviewing earlier versions would
% help, but it appears that I started to convert the acoeffs to the power
% series version, but didn't complete the conversion, then 
% forgot that the chX were not part of the definition.
% I hardly know what to say -- this makes it hard to trust anything I do. I
% think the damage should be contained to the Grammian computations. Even
% though I did check the new LTV behavior, apparently, checking that
% qualitative behavior (whether the LTV appears to follow the nonlinear
% model for small perturbations) was not nearly a strong enough test. 
% Comments: 1) The conversion to Rappel (between v.13 and 14) was a disaster
% as far as the linearized model was concerned. The fact that I didn't
% comprehend that the conversion was incomplete is frightening. Note that
% another error (epsil1 nonzero) occured in v.14, adding to the errors in
% the Jacobian. I can't account for what I was doing around this time to
% lead to so many errors. I can understand that the conversion was rushed
% ahead of the conference, but after the conference, it appears I forgot to
% check for errors in the archived version of 11/4. The previous version 
% (scaled units, v.13) appears to have had fewer problems, outside of the older sech^2
% error that was corrected in an earlier version. 
% 2) I need better methods for error-checking the linearization sections of
% the code. Any suggestions? 
% 3) Since I'm having so much trouble getting this part right, it doesn't
% seem like a good idea to pursue this approach for more complicated
% models, given that I made so many mistakes in the two-variable model. I hope
% that a numerical Jacobian approach can help avoid some of the mistakes,
% since it will only depend on the nonlinear equations and not a separate
% derivation. This won't help with errors in parameter settings or coding
% logic, though. The only thing I can think of right now is to try to
% hand-compute a few values of the derivatives, and to see if the coded
% versions agree, but I don't know whether this will really help or not. 
% See latest jobtalk.ppt for a discussion of effects of corrections on
% Grammians. 

% To prevent "wraparound" reading of Vd and nd, 
% 1) disable line around 573: if writeintsteps <= size(Vd,2) should be if 0
% 2) revise lines around line 885 (acoeff1, etc.) to use trueindex instead
% of k (except for hprime ref in acoeff1) 

% Rewrote to avoid memory errors for large arrays. Changed function
% argument to bcl.

%%function [Vd,iext,stimindices_new,apedges,rctr,dctr,apds,n,V] = karma_sim_statefbtest(karma_test,varargin)%karma_sim(karma_test,pubyear,stimperiod)
%function [] = karma_sim_statefbtest(karma_test,varargin)%karma_sim(karma_test,pubyear,stimperiod)
%function [] = karma_sim_statefbtest(bcl)%karma_sim(karma_test,pubyear,stimperiod)
%function [] = karma_sim_statefbtest(bcl,Karg)
function [] = karma_sim_statefbtest(bcl,varargin)
% Allow optional input of feedback gain (useful for batch mode) -- Karg =
% varargin{1}
% simulate 2-variable Karma model for 1-D geometry
% use parameter values and functions from either PRL '93 paper or Chaos '94
% paper
% varargin{1} = pubyear (required input)
% varargin{2} = stimperiod (optional)
% if karma_test is nonzero, attempt to reproduce Fig. 3 in '93 or '94 paper
% use values from paper corresponding to 'pubyear': 93 or 94
% if pubyear == 'ext', use 1D external potential version of equations
% outlined in Niels' biodomain model handout
% 6/28/07: changed detection scheme to try to make it work in the case when
% sensor and actuator are in the same cell.
% 7/02/07: added capability for multiple detectors

%clear all;
%close all;

global Vstar Vb Vn Vh b ch0 ch1 ch2 ch3

% load data_80cell_b225_7500_ol/34.mat V n usf
% Vi = V;
% ni = n;
% usfi = usf;

mkdir('data');
% dfl = input('Delete old files? y or n: ','s')
% if dfl=='y'
     delete data/*.mat;
% end

karma_test = 0;
%pubyear = varargin{1};
pubyear = 'ext';
%writeint = 2000; % maximum array length in ms
%writeint = round(2000/bcl)*bcl; % maximum array length in ms
%writeint = round(700/bcl)*bcl; % maximum array length in ms
%writeint = round(1000/bcl)*bcl; % maximum array length in ms, 20 cells
%(was also using 1000 for the single-cell tests, although it wasn't
%necessary)
%writeint = round(3250/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(325/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(500/bcl)*bcl; % maximum array length in ms, 40 cells
%writeint = round(250/bcl)*bcl; % maximum array length in ms, 80 cells
writeint = round(1000/bcl)*bcl; % maximum array length in ms, 80 cells
%writeint = round(7500/bcl)*bcl; % maximum array length in ms
% Confined to be an even number of bcls to help with desired value alignment
% Chosen as an integer number of bcls closest to some arbitrary duration


if pubyear == 93
    deltat = 0.03/16; %ms (deltat not given in paper, just made it small)
elseif pubyear == 94
    %deltat = 0.25; %ms
    deltat = 0.03/16;%0.03/(16*7);%0.25/32; %ms
elseif pubyear == 'ext'
    %    deltat = 0.03/32;%0.03/(16*7);%0.25/32; %ms
    %    dtold=deltat;
    %    deltat = 2e-4;
%    deltat = 0.05/8; %ms
    deltat = 0.05; %ms
    dtold=deltat;
%    deltat = 0.05; %ms
    dtscale = 1;%4; % integer scaling factor
    deltat = 0.05*dtscale; %ms
end
%deltax = .0262; %cm
deltax = .0262*2; %cm

if karma_test
    if pubyear == 93
        fiberlength = 0.216;%2.62; % cm
    elseif pubyear == 94
        fiberlength = 2.62; %.216%6.55;% cm
    elseif pubyear == 'ext'
        disp('Error: External-potential version of equations not set up for Karma test.')
        return;
    end
else
%        fiberlength = 1*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%        fiberlength = 10*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%        fiberlength = 20*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%        fiberlength = 40*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%        fiberlength = 60*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%    fiberlength = 80*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%    fiberlength = 100*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
%    fiberlength = 160*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
if length(varargin) >= 1
%    fiberlength = varargin{1}*deltax
    fiberlength = varargin{2}*deltax
else
    fiberlength = 100*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
end
end
%finaltime = 10; % ms
%finaltime = 12; % ms
%finaltime = 40; % ms
%finaltime = 60; % ms
%finaltime = 80; % ms
%finaltime = 90; % ms
%finaltime = 18.2; % ms
%finaltime = 19.14; % ms
%finaltime = 20; % ms
%finaltime = 20*250; % ms
%finaltime = 2250; % ms
finaltime = 7500;
%finaltime = 4500;
%finaltime = 3250;
%finaltime = 10000;
%finaltime = 2*7500;
%finaltime = 2*7500+4200;
%finaltime = 3*7500;
%finaltime = 5000;
%finaltime = 20000; 
%finaltime = 25000; 
%finaltime = 30000; 
%finaltime = 35000; 
%finaltime = 45000; 
%finaltime = 12000; 
%finaltime = 28000; 
%finaltime = 24200;
%finaltime = 2.5*7500;
%finaltime = 6*250; % ms
%finaltime = 35; % ms
%finaltime = 15; % ms
% derived quantities
numrep = ceil(finaltime/writeint); % number of cycles based on finaltime and writeint
writeintsteps = round(writeint/deltat);
numpart = round(fiberlength/deltax); % number of partitions
numstep = round(finaltime/deltat); % time steps
%finalblocksteps = numstep-writeintsteps; % length of final interval;
finalblocksteps = writeintsteps - (numrep*writeintsteps-numstep); % length of final interval;
%icfile = 'V_b087_time_0_80_sfi'; % initial condition file % This actually doesn't work unless you can store the correct values of rctr and dctr
poster = 1; % use flag to make larger plots for poster
fs = 24; % poster fontsize
lw1 = 3; % poster linewidth 1
lw2 = 1.5; % poster linewidth 2

if pubyear == 93
    tauN = 1;
    cc_ext = 0;
    cc_int = 0;
    if karma_test
        M=30;
        nB=0.525;
    else
        M=30; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        nB=0.507; % increasing makes more "sensitive"
    end
    b = 1;
    %epsil = .009; % increasing too much means no AP formation; some runtime error if pushed below .007
    epsil1 = .009; % increasing too much means no AP formation; some runtime error if pushed below .007
    epsil2 = .009; % increasing too much means no AP formation; some runtime error if pushed below .007
    % Fourier number
    Fom = epsil1*deltat/(deltax)^2
elseif pubyear == 94
    D = 1.1/100 %.009;%cm^2/ms, 1.1 cm^2/s, diffusion coefficient
    tauV = 2.5/100 %.009%2.5; %ms
    tauN = 250/100 %1%250; %ms
    epsil1 = D;
    epsil2 = tauV;
    cc_ext = 0;
    cc_int = 0;
    nB=1;%.507;
    if karma_test
        M=10;
    else
        M=10;
    end
    Re = 1.204;%500000%
    b = 1-exp(-Re);
    % Fourier number
    Fom = D*deltat/(deltax)^2
elseif pubyear == 'ext'
    %    tauN = 1;
    tauN = 250;
%    tauN = 250*1.05; % ofi still stabilizes, but not as well. OL portion maybe looks less stable, which would be consistent with reduced performance.
%    tauN = 250*0.95; % ofi stabilizes better. OL portion may look more stable, which would be consistent
    Dg = 0.0108;
    De = 5*Dg; % De=1/(c*re), inversely proportional to external fluid resistance
    c = 0.009;
    %    cc_ext = max([0 (Dg/(c*(Dg + De)))]);
    %    cc_int = max([0 (De/(c*(Dg + De)))]);
    cc_ext = 1; % uF/cm^2 ? Rappel, et al.
%    cc_ext = 1*1.05; % ofi still stabilizes
%    cc_ext = 1*0.95; % ofi still stabilizes
    cc_int = 1;
    if karma_test
        M=30;
        nB=0.525;
    else
        M=10; % lowering to 10 softens repolarization edges but doesn't amplify alternans
%        M=10*1.05; % ofi seems to stabilize more quickly, although ss-apd value may not be as low. Would it make sense for increasing M to make the system "more stable"?
%        M=10*0.95; % ofi seems to stabilize better than default case. Would it make sense for both increasing and decreasing M to make the system "more stable"?
        %        nB=0.707; %0.680;% increasing makes more "sensitive"
        nB=1; % nB=0.707;
%        nB=1*1.05; %ofi does not stabilize. OL still looks like some sort of alternans, but appears to be transitioning to 2:1 before control is turned on. ofi/5 appears to help restore alternans briefly, but system appears to lapse back toward 2:1 later on. So, ofi/5 may be too weak, but ofi doesn't work either. 
%        nB=1*0.95; %ofi seems to stabilize better? Checked a few cells and tracking looks OK. 
        %Need to check a tracking error measure to be sure. Presumably,
        %judging from the behavior before fbon_time, lowering nB effectively converts the
        %system to a convergent OL trajectory. This is probably the
        %explanation of whether a disturbance in a particular direction
        %"helps" or not; i.e., if a particular direction makes the OL system
        %more stable, then the perturbed CL result could reasonably be better than
        %the unperturbed CL result. 
    end
    %    b = 1;
    %    b = 1/(tauN*(1+2.4)/2.4);
    b = 1/((1+2.4)/2.4); % I think there is a typo at the end of Rappel's paper, where all the terms in g(V,n) should be
%    b = 1.05*1/((1+2.4)/2.4); % ofi and ofi/5 do not stabilize. OL looks like alternans, although it may be transitioning to 2:1 before control is turned on. 
%    b = 0.95*1/((1+2.4)/2.4); % ofi stabilizes more quickly
    % divided by tauN. If so, this should be the correct expression for b.
    % If I use the other value of b implied by the paper (which includes
    % tauN), then the resulting n trajectories look wrong (too steep due to
    % domination of the step function) and I wasn't able to find stimulus
    % values that would allow the cell to fire "normally".
    %epsil = .009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = max([0 Dg*De/(Dg+De)]); %.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.0053; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 9.3*0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.001/4; % cm^2/ms
    %    epsil1 = 0.001/2; % cm^2/ms
    epsil1 = 0.001; % cm^2/ms
    epsil1 = 0.025
%    epsil1 = 0.001*1.1; % cm^2/ms % ofi still stabilizes
%    epsil1 = 0.001*0.9; % cm^2/ms % ofi does not stabilize
%    epsil1 = 0.001*1.05; % cm^2/ms % ofi still stabilizes
%    epsil1 = 0.001*0.95; % cm^2/ms % ofi still stabilizes
%    epsil1 = 0; % cm^2/ms
    %    epsil1 = 0.001*0.75; % cm^2/ms at 225, dur=7, alt is slightly worse but not discordant
    %    epsil1 = 0.001*0.5; % cm^2/ms at 225, dur=7, alt is not discordant (may be slightly less severe than 0.75 case? why?)
    %    epsil1 = 0.001*0.25; % cm^2/ms at 225, dur=7, alt is not discordant (but more spatial variation than 0.5 case, and cell 1 looks more stable). Looks similar for dur=5 and 3, except alt. amp is somewhat reduced
    %    epsil1 = 0.001*0.1; % cm^2/ms at 225, dur=3, alt is starting to look discordant around
    %t=5000ms (end of run), but alt envelope is decaying overall. Running
    %longer shows discordancy between 4000 and 8000, but still decaying.
    %Switching to 215, discordancy shows up around 4000. Not clear if it
    %will grow or decay overall. 205: discordancy starting around 3500, may
    %be expanding. 195: discordancy around 3250. 193: growing discordancy
    %starting around 2000. 190: discordant between 1800 and 2200,
    %2:1 thereafter. 185: concordant to discordant (briefly) back to
    % concordant then to 2:1 around 3200. 175: 2:1 right away
    %(around 1000)
    %    epsil2 = c;
    epsil2 = 5; 
%    epsil2 = 5*1.05; % ofi still stabilizes
%    epsil2 = 5*.95; % ofi still stabilizes
    % Fourier number
    Fom = epsil1*deltat/(deltax)^2
end

% Vstar = 1.5415; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
% Vn = 1; % lowering seems to make alternans go away faster; making too high causes breakup/block
% Vh = 3; % lowering seems to make alternans go away faster, effect of increasing is not clear
Vstar = 1.88; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.88*1.05; % ofi still stabilizes, but performance looks worse
%Vstar = 1.88*0.95; % ofi still stabilizes
Vn = -60; % lowering seems to make alternans go away faster; making too high causes breakup/block
%Vn = -60*1.05; % ofi stabilizes, but not as well
%Vn = -60*0.95; % ofi stabilizes better
Vh = -21; % lowering seems to make alternans go away faster, effect of increasing is not clear. You should comment out Vh since it isn't used anymore
Vb = -85; % mV (offset from Rappel) 
%Vb = -85*1.05; % APD convergence looks "better" but this is misleading, since there is serious (40deg?) phase error (but not total inversion) in the V vs. Vd tracking 
%Vb = -85*0.95; % ofi does not stabilize; neither does ofi/5 (clearly makes things worse). OL does not look like alternans. In fact, it is closer to 1:1, but with a higher "resting" potential, around -77 or -78, which is higher than the specified, perturbed Vb. 
% coefficients in h(V)
%ch1 = 1;
%ch2 = 0;
%ch3 = 2;
%h = 72.83 - .83*V - .04976*V.^2 - (3.52e-4)*V.^3;
%ch1 = 45;
%ch2 = 85;
%ch3 = 56;
%h=((1-tanh((V-Vh)/ch1)).*(V+ch2).^2)/ch3;
% Redefine ch* to match power series coefficients 1/06/09
ch0 = 72.83; 
%ch0 = 72.83*1.05; % ofi does not stabilize. OL does not look like alternans
%ch0 = 72.83*0.95; % ofi probably does not stabilize. apds converge faster, but current looks crazy
%ch0 = 72.83*1.005; % ofi/5 may help. OL looks like alternans. Haven't tried default ofi
%ch0 = 72.83*0.995; % ofi/5 may help. OL looks like alternans. Haven't tried default ofi
ch1 = -0.83;
%ch1 = -0.83*1.05; % ofi does not stabilize. OL does not look like alternans. 
%ch1 = -0.83*0.95; % ofi probably does not stabilize. apds converge faster, but current looks crazy
%ch1 = -0.83*1.005; % ofi/5 may help. OL looks like alternans. Haven't tried default ofi
%ch1 = -0.83*0.995; % ofi/5 may help. OL looks like alternans. Haven't tried default ofi
ch2 = -0.04976;
%ch2 = -0.04976*1.05; % ofi does not stabilize. this is the second most extreme shift I've seen
%ch2 = -0.04976*0.95; % ofi does not stabilize. this is among the most extreme shifts I've seen. Confirmed: OL does not exhibit normal APs, let alone alternans
%ch2 = -0.04976*0.99; % OL does not exhibit normal APs, although unlike the 0.95 case, there are recognizable APs. I'm not sure how to classify the type of arrhythmia. It looks like partial block in some places, but hasn't settled into any recognizable APD rhythm by 15000ms. Adding default ofi control, or gains reduced by factor of 5, doesn't improve anything (current goes crazy)
%ch2 = -0.04976*0.999; % OL looks very similar to unperturbed case (disco alt). 
%ch2 = -0.04976*0.995; % OL does look like some sort of alternans (maybe concordant), but the default ofi controller just goes bonkers (it actually makes things worse). Reducing all gains by factor of 10 is too much (does nothing) and by 2 is too little (still destabilizes system). Reducing by a factor of 5 seems more effective (appears to be some benefit to control, after a while) 
%ch2 = -0.04976*1.005; % OL does look like some sort of alternans. Reducing ofi gains by factor of 5 appears to give some benefit. Didn't test default ofi controller. 
ch3 = -3.52e-4;
%ch3 = -1.05*3.52e-4; % ofi does not stabilize. this is among the most extreme shifts I've seen. OL does not look like alterans. 
%ch3 = -0.95*3.52e-4; % ofi does not stabilize. this is among the most extreme shifts I've seen
%ch3 = -0.995*3.52e-4; % Haven't tried default ofi. ofi/5 may be doing some good. The OL case still looks like alternans. 
%ch3 = -1.005*3.52e-4; % Haven't tried default ofi. ofi/5 may be doing some good. The OL case still looks like alternans. 

deton = 1; % turn on AP detectors
%detloc = [0.1 0.15 1.3 2.5 2.55]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [deltax]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [0.1 0.15 0.2 0.25]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [0.2]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [0.11]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [0.1 0.25 0.40]; % location of measurement electrodes, measured
%from upstream end (x=0) (cm)
detloc = deltax*(1:numpart); % location of measurement electrodes, measured from upstream end (x=0) (cm)
%detloc = [0.5]; % location of measurement electrodes, measured from upstream end (x=0) (cm)
numdet = length(detloc);
detvec = [1:numdet]';
if karma_test
    detthresh = -40*ones(numdet,1);%2*ones(numdet,1); % detection threshold for V
    % -- should be less than stimulus size if sensor and actuator are co-located
else
    detthresh = -40*ones(numdet,1);%2*ones(numdet,1); % detection threshold for V
    % -- should be less than stimulus size if sensor and actuator are co-located
end

fbon = 'sf';%0;% %0%'lyap'%'pyra'%; turn on feedback (0 for none, 'lyap' for Lyapunov, 'pyra' for pyragas)
cgain = 0.15; %0.33;% feedback control gain; 0.1 works well for nB=0.707, bcl=0.95, M=30
fbdet = 1% index of detector used for feedback measurments
stimoutdet = 4; % index of detector used to measure downstream current extraction point

stimloc = deltax;%0.1 % location of stimulating electrode, measured from upstream end (x=0) (cm)
%stimheight = (deltat/dtold)*1.5/5; %1.5/5;%
%stimheight = (deltat/dtold)*1.6; %1.5/5;%
stimheight = (deltat/dtold)*2.5; %1.5/5;%
%stimheight = 1.6; %(deltat/dtold)*0.2;%mV %Rappel et al., 10-20mV? Using Rappel's h(V), the "minimum" stimulus for
% a clean-looking AP seems to be either 0.2 height for 1 ms or 0.1 for 2
% ms, when deltat=0.05/8. Stimheight should go up by a factor of deltat/dtold if
% deltat is increased. For the old approximate h(V), I expect it to be higher, but I haven't
% checked yet. The units of iext as defined in the code should be uF/cm^2
%stimduration = 1/250; %deltat% assumptions: 1) period = 1 corresponds to 250ms, 2) typical stimulus duration is 1-5ms
stimduration = 7;%used for >=10cells; In new units, a stimulus duration of 20 timesteps (old duration) would be equiva
%stimduration = 7;%used for >=10cells; In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
%stimduration = 5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
%stimduration = 3; % 
%stimduration = 2.5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
% 1 cell: A stimduration of 1ms will work at 1.6 height for 1 cell, but need 2.5ms duration to
% get the depolarization edge "all the way up" without changes in slope
% 2 cells: The upstroke in cell 1 takes a dip around -30mV using the 1cell
% settings
% 3 cells: The upstroke in cell 1 takes a dip around -45mV using the 1cell
% settings. If dur is increased to 4ms, the dip moves to -15mV.
% 10 or 20 cells: some dur >5 but <10 ms is probably ideal
% 20 cells: 7ms seems ok (1.6 height)
stimoutloc = 0.2%2*deltax; %2.5;% location where current is extracted
stimstart = deltat; % start time for stimulation (ms)
%if length(varargin) == 1
% The following appears to be the period for
% vtest_periodic_T30_M30_nB707:
%    stimperiod = (deltat+0.94968750000000)*[1 1]; % stimuation period (ms) % with nB=.707 and M=30%
%    stimperiod = [0.91 0.91] % stimuation period (ms) % with nB=.707 and M=30%
%    stimperiod = [0.85 0.85] % stimuation period (ms) % with nB=.707 and M=30%
%    stimperiod = [0.87 0.87] % stimuation period (ms) % with nB=.707 and M=30%
%    stimperiod = 0.97*[1 1] % stimuation period (ms) % with nB=.707 and M=30%
% Warning: Bifurcation behavior depends heavily on the number of cells, as well as
% the values of epsil1, stimheight and stimduration.
%stimperiod = 325*[1 1]; % ms, Rappel et al. Decaying alternans for 1 cell (height,dur) (1.6,1)
%stimperiod = 300*[1 1]; % ms, Rappel et al. Decaying alternans for 1 cell (height,dur) (1.6,1)
%stimperiod = 260*[1 1]; % ms, Rappel et al. Decaying alternans for 1 cell (height,dur) (1.6,1)
%stimperiod = 275*[1 1]; % ms, Rappel et al. Decaying alternans for 1 cell (height,dur) (1.6,1)
%stimperiod = 250*[1 1]; % ms, Rappel et al. Decaying alternans for 1 cell (height,dur) (1.6,1)
%stimperiod = 200*[1 1]; % ms, Rappel et al. Goes from almost no alternans,
%to a sharp expansion, to 2:1 around 4300ms for 1 cell (1.6,1)
%stimperiod = 225*[1 1]; % ms, Rappel et al. Goes from expanding alt to 2:1
stimperiod = bcl*[1 1];
%stimperiod = 215*[1 1]; % ms
%stimperiod = 205*[1 1]; % ms
%stimperiod = 195*[1 1]; % ms
%stimperiod = 175*[1 1]; % ms
%stimperiod = 185*[1 1]; % ms
%stimperiod = 190*[1 1]; % ms
%stimperiod = 193*[1 1]; % ms
%around 3000ms (1.6,1)
%stimperiod = 190*[1 1]; % ms,
% For Rappel, the critical point for 1 cell appears to be between 225 and
% 250 ms for stimheight = 1.6, dur=1
%    stimperiod = 0.87*[1 1] % stimuation period (ms) % with nB=.707 and M=30%
%    stimperiod = 1.12*[1 1] % stimuation period (ms) % with nB=.707 and M=30%
%with nB=.707, M=30 for a single cell, "stable" (nonincreasing, nondecreasing) alternans point, if it exists, is between
%   .915 and .92 (tmin = 0.88)
%For a fiber, bcl=0.917 does *not* produce stable alternans: goes from alternans to 2:1 rhythm after a while
%with nB=.707, M=10 for a fiber, "stable" discordant alternans appears
%to occur between bcl=0.8 and 0.85; blocking occurs at least by 0.75
%Karma says (93): tmin = log(nB/(1-nB))
% elseif length(varargin) > 1
%     stimperiod = varargin{2};
% end

stimtransamt = 0.5; % step amount (ms) by which period is increased or decreased during transition

detlocindices = round(detloc/deltax); % cell indices of measurement electrode locations
stimdurationsteps = round(stimduration/deltat);
stimlocindex = round(stimloc/deltax); % cell index of stimulating electrode location
stimoutlocindex = round(stimoutloc/deltax); % cell index of location where current is extracted
stimstartindex = round(stimstart/deltat); % time index for initial stimulus
stimperiodsteps = round(stimperiod/deltat); % timesteps between stimuli
if ~exist('icfile')
    %fbon_time = 10*stimperiodsteps(1); %round(25/deltat);% time index at which control is switched on
    %fbon_time = round(2200/deltat); %round(25/deltat);% time index at which control is switched on
    %fbon_time = round(1320/deltat); %round(25/deltat);% time index at which control is switched on
    % Discordant:
    %fbon_time = round(2275/deltat); %round(25/deltat);% time index at which control is switched on
    %fbon_time = round(3500/deltat); %round(25/deltat);% time index at which control is switched on
    %fbon_time = round(4500/deltat); %round(25/deltat);% time index at which control is switched on
   fbon_time = round(2800/deltat); %round(25/deltat);% time index at which control is switched on
%%    fbon_time = round(7000/deltat); %round(25/deltat);% time index at which control is switched on
%    fbon_time = round(12000/deltat); %round(25/deltat);% time index at which control is switched on
%    fbon_time = round(4000/deltat); %round(25/deltat);% time index at which control is switched on
%    fbon_time = round(10500/deltat); %round(25/deltat);% time index at which control is switched on
    %fbon_time = 1; %round(25/deltat);% time index at which control is switched on
    %fbon_time = round(2500/deltat); %round(25/deltat);% time index at which control is switched on
else
    fbon_time = 1; % time index at which control is switched on
end
%fbon_time = 28*stimperiodsteps(1); % time index at which control is
%switched on (need to wait longer if you want discordant alternans to fully
%develop before the control is turned on)

stimtransstep = round(stimtransamt/deltat); %

if stimperiod(1)>stimperiod(end)
    stepsign = -1;
else
    stepsign = 1;
end

transperiods = stimperiodsteps(1):(stepsign*stimtransstep):stimperiodsteps(end);
transperiodsteps = cumsum(transperiods);
%stimindices = stimstartindex + [0  transperiodsteps transperiodsteps(end)+cumsum(stimperiodsteps(end)*ones(1,300))];
stimindices = stimstartindex + [0  transperiodsteps transperiodsteps(end)+cumsum(stimperiodsteps(end)*ones(1,ceil(numstep/stimperiodsteps(1))))];

stimindices_new = stimindices;

%V(:,1) = zeros(numpart,1); % initial condition
% ***
% V = zeros(numpart,numstep);
V = zeros(numpart,writeintsteps);
% if karma_test
%     V(1,1) = 2; % initiate pulse that will travel around ring.
%     % If we want to change the initiation point, then the periodic BCs must
%     % change also. Right now, a pulse is started at the leftmost end, after
%     % which the ends of the fiber are tied together: V_dn = V(end,...),
%     % and Vup = Vdn
% end
% %V(4,1) = 2; % cell 4 is approx .1 cm from upstream end, for deltax = 0.0262 cm
% %n(:,1) = zeros(numpart,1); % initial condition
% n = zeros(numpart,numstep);
n = zeros(numpart,writeintsteps);
% ***
%  load xfzero_M10_BCL1_5;
% load  xf_ninthBCL_M10_BCL1_5;
% Xf=x;
%  V(:,1) = Xf(1:numpart); % initial condition
%  n(:,1) = Xf((numpart+1):end); % initial condition
%  V(1,1) = 2.5; % initial condition
%  n(1,1) = 0.5; % initial condition
V(:,1) = -85; % initial condition
n(:,1) = 0; % initial condition
% %load diff0045_0_20_ic
%load diff0045_0_18_2_ic
% load diff0045_0_18_2_ic
% V(:,1) = Vinit;
% n(:,1) = ninit;
% Load ICs from SFI run
if exist('icfile')
    eval(['load ' icfile ' Vinitnew ninitnew iextinitnew'])
    V(:,1) = Vinitnew(:,2);
    n(:,1) = ninitnew(:,2);
end

% *** iext=zeros(numpart,numstep); % external current
%temp=zeros(numpart,numstep); % external current
% *** iint=zeros(numpart,numstep); % internal current
if pubyear == 'ext'
    if fbon
        stimindices_short = stimindices(stimindices <= (fbon_time + 2*stimperiodsteps(1)));
    else
        stimindices_short = stimindices(stimindices <= numstep);
    end
    %     if cc_ext
    %         iext(stimlocindex,stimindices_short)= -stimheight/deltat/cc_ext;
    %         if epsil1
    %             iext(stimoutlocindex,stimindices_short) = -iext(stimlocindex,stimindices_short);
    %         end
    %     end
end

% generate or load desired V profile
%Vd = zeros(numpart,numstep);
%nd = zeros(numpart,numstep);
%iextd = zeros(numpart,numstep); % desired external current
iextd = zeros(numpart,writeintsteps); % desired external current
ctr = 1;
if cc_ext
    while stimindices(ctr) < writeintsteps%numstep
        iextd(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -stimheight/deltat/cc_ext; % desired external current
        if 0%epsil1 & ~strcmp(fbon,'lyap') % more than one cell
            iextd(stimoutlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -iextd(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps));
        end
        %    Vd(detlocindices(2),stimindices(ctr):(stimindices(ctr)+69)) = 3.68;
        %    Vd(detlocindices(4),stimindices(ctr):(stimindices(ctr)+69)) = 3.68;
        if ~fbon %| strcmp(fbon,'sf')
            iext(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -stimheight/deltat/cc_ext; % desired external current
            if 0%epsil1 % more than one cell
                iext(stimoutlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -iext(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps));
            end
        end
        ctr = ctr+1;
    end
end
% if fbon & ~strcmp(fbon,'sf')
%     ctr = 1;
%     if cc_ext
%         while stimindices_short(ctr) < fbon_time + stimperiodsteps(1)
%             iext(stimlocindex,stimindices_short(ctr):(stimindices_short(ctr)+stimdurationsteps)) = -stimheight/deltat/cc_ext; % desired external current
%             if 0%epsil1 & ~strcmp(fbon,'lyap') % more than one cell
%                 iext(stimoutlocindex,stimindices_short(ctr):(stimindices_short(ctr)+stimdurationsteps)) = -iext(stimlocindex,stimindices_short(ctr):(stimindices_short(ctr)+stimdurationsteps));
%             end
%             ctr = ctr+1;
%         end
%     end
% end

%iext(:,1:6) = zeros(numpart,6);

%load('jon_desval_oncell_b091.mat');
%load('jon_desval_2cell_b091.mat');
%load('jon_desval_2cell_b091.mat');
%load('desval_twocell_b087_d0045.mat');
%load('desval_3cell_b087_d0053.mat');
%load('desval_3cell_b087_d0090.mat');
% load('desval_3cell_b087_d0090_shortspike.mat');
% if epsil1==0.0053
% load('desval_3cell_b087_d0053_shortspike.mat');
% end
if 1%epsil1==0.009
    %load('desval_4cell_b087_d0090_shortspike.mat');
    %load('desval_40cell_b087_d0090_shortspike.mat');
    %eval('load desval_40cell_b087_d0090_shortspike_finaltime60.mat Vd nd');
    %desfile = 'desval_20cell_b087_d0090_shortspike_finaltime100.mat';
    %desfile = 'desval_10cell_b225_d001_finaltime40';
    %desfile = 'desval_10cell_b225_d001_finaltime40_no1data'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = 'desval_10cell_b225_d001_finaltime40_allsample'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = 'desval_2cell_b225_d001_finaltime40_allsample_shortspike'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = 'desval_2cell_b225_d001_finaltime40_cutout_shortspike5'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = 'desval_40cell_b087_d0090_shortspike.mat';
    %desfile = 'desval_2cell_b097_d0090_shortspike_finaltime30.mat';
%    desfile = 'desval_10cell_b225_d001_finaltime40_cutout'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = 'desval_10cell_b190_d0001_finaltime40_cutout'; % attempt to use only fiber data for more accurate depolarization slopes
    %desfile = '1cellOL_testA_b225';
%    desfile = '1cellOL_testA_b325';
%    desfile = '1cellOL_testA_b325_3250';
    %desfile = '1cellOL_testA_b250';
%    desfile = 'desval_1cell_b225_d001_finaltime40'; 
%    desfile = 'desval_80cell_b225_d001_finaltime5_cutout'; % attempt to use only fiber data for more accurate depolarization slopes
%    desfile = 'desval_20cell_b165_finaltime5_glitch'; % using _glitch or interp profile seems to produce similar results. Glitches may be too fast or small to induce tracking. 
%    desfile = 'desval_20cell_b165_finaltime5_lininterp'; %
%    desfile = 'desval_40cell_b165_finaltime5_lininterp'; %
%    desfile = 'desval_40cell_b165_finaltime5_lininterp_dapd119'; %
%    desfile = 'desval_60cell_b165_finaltime5_lininterp'; %
%    desfile = 'desval_60cell_b165_finaltime5_lininterp_dapd119'; %
%    desfile = 'desval_60cell_b165_finaltime5_lininterp_dapd115'; %
%    desfile = 'desval_60cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_20cell_b180_finaltime5_lininterp_dapd125'; %
%    desfile = 'desval_20cell_b180_finaltime5_lininterp_dapd110'; %
%    desfile = 'desval_40cell_b180_finaltime5_lininterp_dapd110'; %
%    desfile = 'desval_60cell_b180_finaltime5_lininterp_dapd110'; %
%    desfile = 'desval_80cell_b180_finaltime5_lininterp_dapd110'; %*
%    desfile = 'desval_80cell_b180_finaltime5_lininterp_dapd120'; %
%    desfile = 'desval_80cell_b180_finaltime5_lininterp_dapd105'; %*
%    desfile = 'desval_100cell_b180_finaltime5_lininterp_dapd110'; %*
%    desfile = 'desval_120cell_b180_finaltime5_lininterp_dapd110'; %*
%    desfile = 'desval_140cell_b180_finaltime5_dapd110'; %*
%    desfile = 'desval_160cell_b180_finaltime5_dapd110'; %*
%    desfile = 'desval_20cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_40cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_20cell_b195_finaltime5_lininterp'; %
%    desfile = 'desval_40cell_b195_finaltime5_lininterp'; %
%    desfile = 'desval_60cell_b195_finaltime5_lininterp'; %
%    desfile = 'desval_80cell_b195_finaltime5_lininterp'; %
%    desfile = 'desval_100cell_b195_finaltime5_lininterp'; %
%    desfile = 'desval_120cell_b195_finaltime5'; %
%    desfile = 'desval_140cell_b195_finaltime5'; %
%    desfile = 'desval_160cell_b195_finaltime5'; %
%    desfile = 'desval_20cell_b210_finaltime5_lininterp'; %
%    desfile = 'desval_40cell_b210_finaltime5_lininterp'; %
%    desfile = 'desval_60cell_b210_finaltime5_lininterp'; %
%    desfile = 'desval_80cell_b210_finaltime5_lininterp'; %
%    desfile = 'desval_100cell_b210_d001_finaltime5';
%    desfile = 'desval_120cell_b210_finaltime5'; %
%    desfile = 'desval_140cell_b210_finaltime5'; %
%    desfile = 'desval_160cell_b210_finaltime5'; %
%    desfile = 'desval_20cell_b225_finaltime5_lininterp';
%    desfile = 'desval_40cell_b225_finaltime5_lininterp';
%    desfile = 'desval_60cell_b225_finaltime5_lininterp';
%    desfile = 'desval_80cell_b225_finaltime5_lininterp';
%    desfile = 'desval_100cell_b225_d001_finaltime5'; 
%    desfile = 'desval_120cell_b225_finaltime5_lininterp';
%    desfile = 'desval_140cell_b225_finaltime5_lininterp';
%    desfile = 'desval_160cell_b225_finaltime5_lininterp';
%    desfile = 'desval_20cell_b240_finaltime5_lininterp';
%    desfile = 'desval_40cell_b240_finaltime5_lininterp';
%    desfile = 'desval_60cell_b240_finaltime5_lininterp';
%    desfile = 'desval_80cell_b240_finaltime5_lininterp';
%    desfile = 'desval_100cell_b240_finaltime5_lininterp';
%    desfile = 'desval_120cell_b240_finaltime5_lininterp';
%    desfile = 'desval_140cell_b240_finaltime5_lininterp';
%    desfile = 'desval_160cell_b240_finaltime5_lininterp';
%    desfile = 'desval_80cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_100cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_120cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_140cell_b165_finaltime5_lininterp_dapd100'; %
%    desfile = 'desval_160cell_b165_finaltime5_lininterp_dapd100'; %
%    eval(['load ' desfile ' Vd nd dapd']);
%    eval(['load ' desfile ' Vd dapd']);
%    Vd = zeros(numpart,writeintsteps+1);
%    nd = zeros(size(Vd)); 
%    eval(['load ' desfile ' Vd nd']);
    desfile = 'markdata_proc_50_70'; 
    Vinterp = 0
    eval(['load ' desfile ' Vinterp']);
    Vd = zeros(numpart,writeintsteps+1);
    % assume that for 12/12/07 data, there are 2 PF electrodes (nos. 1 and 3) located at
    % 3mm and 14mm from the proximal end of the sample.
    sensorlocindices = round([0.3 1.4]/deltax); 
    sensorindices = [1 3]; % indices of the sensors relative to the data file
    Vd(sensorlocindices,:) = Vinterp(1:(writeintsteps+1),sensorindices)'; 
    nd = zeros(size(Vd)); 
    if ~exist('dapd')
        dapd = -1;
    end

    if exist('dtscale');
        Vd = Vd(:,1:dtscale:end); 
        nd = nd(:,1:dtscale:end); 
    end
end

% Expand Vd sizes if necessary
if size(Vd,1) < numpart
    %Vdt = zeros(numpart,numstep);
    %ndt = zeros(numpart,numstep);
    %Vdt(1:size(Vd,1),:) = Vd;
    %ndt(1:size(nd,1),:) = nd;
    Vdt = zeros(numpart,writeintsteps+1);
    ndt = zeros(numpart,writeintsteps+1);
    Vdt(1:size(Vd,1),:) = Vd(:,1:writeintsteps+1);
    ndt(1:size(nd,1),:) = nd(:,1:writeintsteps+1);
    Vd = Vdt;
    nd = ndt;
    clear Vdt ndt
end

% comment out condition if you want Vd reading not to wrap around
if writeintsteps <= size(Vd,2)%numstep <= size(Vd,2)
    %    Vd = Vd(1:end,1:numstep);
    %    nd = nd(1:end,1:numstep);
    if ~exist('icfile')
        %        Vd = Vd(1:numpart,1:numstep);
        %        nd = nd(1:numpart,1:numstep);
        % This should be the default case; des traj are shortened to
        % writeinsteps+1 length
        Vd = Vd(1:numpart,1:writeintsteps+1);
        nd = nd(1:numpart,1:writeintsteps+1);
    else
        % 10 is chosen arbitrarily as a stimulus index that occurs after the
        % default fbon_time. If Vd is not sufficiently long, this won't work
        Vd = Vd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
        nd = nd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
    end
else
    %    Vd = Vd(1:end,:);
    %    nd = nd(1:end,:);
    if ~exist('icfile')
        Vd = Vd(1:numpart,:);
        nd = nd(1:numpart,:);
    else
        % 10 is chosen arbitrarily as a stimulus index that occurs after the
        % default fbon_time. If Vd is not sufficiently long, this won't work
        Vd = Vd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
        nd = nd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
    end
end
% figure
% hold on;
% plot(Vd')
% plot(nd','g')
% plot(iextd','c')
% figure
% hold
% plot((1:length(Vinterp))*deltat,Vinterp(:,sensorindices),'b')
% for ii = 1:numrep
%     plot([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)]*deltat, Vinterp([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)],sensorindices)','g--')
% end

if size(Vd,2) < writeintsteps %numstep
    % Interpolate -desired quantities if necesary
    % correct shift induced by initial times
    Vdnew = interp1(dtold*(1:size(Vd,2))-dtold+deltat,Vd',deltat*(1:numstep),'linear',0);
    ndnew = interp1(dtold*(1:size(nd,2))-dtold+deltat,nd',deltat*(1:numstep),'linear',0);

    figure
    plot(dtold*(1:size(Vd,2))-dtold+deltat, Vd')
    hold
    plot(deltat*(1:numstep),Vdnew,'g.')
    plot(deltat*(1:numstep),iextd,'c')

    figure
    plot(dtold*(1:size(nd,2))-dtold+deltat, nd')
    hold
    plot(deltat*(1:numstep),ndnew,'g.')
    plot(deltat*(1:numstep),iextd,'c')

%    Vd = Vdnew';
%    nd = ndnew';
%    Vd = Vdnew; % may need to change back for newer versions of matlab
%    nd = ndnew;
    Vd = [Vdnew zeros(numpart,1)]; % this was added to fix an out-of-bounds error that occurs when writeindex == finaltime 
    nd = [ndnew zeros(numpart,1)];
end

% state feedback term
usf = zeros(1,writeintsteps);
if exist('icfile')
    iext(:,1) = iextinitnew(:,2);
    usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
end
%%%%%%%%%%%%%%%%%%
% V(:,1) = Vi(:,end);
% n(:,1) = ni(:,end);
% usf(1) = usfi(end);
%%%%%%%%%%%%%%%%%%
% nerror = nd(:,1:numstep)-n;
% Verror = Vd(:,1:numstep)-V;
nerror = nd(:,1:writeintsteps)-n;
Verror = Vd(:,1:writeintsteps)-V;
% % Vexx=zeros(numpart,numstep);
% % ferror=zeros(numpart,numstep);
% % gerror=zeros(numpart,numstep);
% % dW1=zeros(1,numstep);
% % dW2=zeros(1,numstep);
% % dW3=zeros(1,numstep);
% % dW4=zeros(1,numstep);
% % W1=zeros(1,numstep);
% % W2=zeros(1,numstep);
% % W=zeros(1,numstep);
% fb=zeros(1,numstep);
% %c1=1; c2=1; % cost function weights
% c1=ones(numpart,1);
% c2=ones(numpart,1); % cost function weights
% %c1(stimlocindex)= 10e-10;
% % c2(stimlocindex)=10e-5;
% z=zeros(numpart,numstep);
% cc=zeros(1,numstep);
% b1=zeros(1,numstep);
% %ac=c1*deltax*(deltat^2)*(cc_ext)^2;
% ac=c1(stimlocindex)*deltax*(deltat^2)*(cc_ext)^2;
% % state feedback quantities:
% A=zeros(2*numpart,2*numpart);
% B1=zeros(2*numpart,1);
% B2=zeros(2*numpart,1);
% %B1(1,1) = cc_ext*deltat;
% %B2(2,1) = cc_ext*deltat;
% B1(1,1) = -cc_ext*deltat; % changed on 7/29/08 to match convention in karma_sim_1cell.m
% B2(2,1) = -cc_ext*deltat;
% %K1=ones(1,2*numpart);
% %K2=ones(1,2*numpart);
% %KIred1=zeros(1,4);
KIred1=zeros(1,2*numpart);
% Ared=zeros(numpart,numpart);
% B1red=zeros(numpart,1);
% B2red=zeros(numpart,1);
% B1red=B1(1:2);
% B2red=B2(1:2);
% Kred=ones(1,numpart);
% C1 = [1 0 0 0]; % measure upstream
% C2 = [0 1 0 0]; % measure downstream
% % state feedback with integrator
% AIred = zeros(2*numpart,2*numpart);
% AIred(1:2,1:2) = eye(2);
% BIred1 = [B1red; B1red];
% BIred2 = [B2red; B2red];
% AI = zeros(4*numpart,4*numpart);
% AI(1:4,1:4) = eye(4);
% BI1 = [B1; B1];
% BI2 = [B2; B2];
%rankms = zeros(1,numstep);
%rankctr = zeros(1,numstep);
%*** un=zeros(numpart,numstep); % optional feedback term added to n-dynamics

% eiga = zeros(2*numpart,numstep);
% eigata = zeros(2*numpart,numstep); % check eigenvalues of A'A (LTV stability criterion)
% eigared = zeros(numpart,numstep);
% eigacl1 = zeros(2*numpart,numstep);
% eigatacl = zeros(2*numpart,numstep); % check eigenvalues of A'A for CL (LTV stability criterion)
% eigaclred = zeros(numpart,numstep);
% eigacl2 = zeros(2*numpart,numstep);
% eigaclred = zeros(2*numpart,numstep);
% eigctrb1 = zeros(2*numpart,numstep);
% eigctrb2 = zeros(2*numpart,numstep);
% eigctrb1red = zeros(numpart,numstep);
% eigctrb2red = zeros(numpart,numstep);
% eigctrbired1 = zeros(2*numpart,numstep);
% eigctrbired2 = zeros(2*numpart,numstep);
% eigctrbi1 = zeros(4*numpart,numstep);
% eigctrbi2 = zeros(4*numpart,numstep);
% eigctrbiredbar1 = zeros(3,numstep);
% eigaired = zeros(2*numpart,numstep);
% eigairedcl = zeros(2*numpart,numstep);
%
%eigobsv1 = zeros(2*numpart,numstep);
%eigobsv2 = zeros(2*numpart,numstep);

% kgainswr = zeros(numpart,numstep);
% kgains = zeros(numpart,numstep);
% load kgains_2cell_b087_Ared_pp0_999;
% load kgains_2cell_b087_Ared_pp0_995;
% load kgains_2cell_b087_Ared_pp0_99;
% evectaclwr = zeros(numpart,numstep);
% evectacl = zeros(numpart,numstep);

%*** usf = zeros(1,numstep);
% if exist('icfile')
%     iext(:,1) = iextinitnew(:,2);
%     usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
% end
capprox = 1; % Gain used for approximation of Heaviside function 1(V-Vn) as 0.5*(1+tanh(capprox*(V-Vn)))
% hprime = zeros(numpart,numstep);
% acoeff1 = zeros(numpart,numstep);
% acoeff2 = zeros(numpart,numstep);
% acoeff3 = zeros(numpart,numstep);
% %s=zeros(2,numstep);

% svector = zeros(2*numpart,numstep);
% %svector(:,1) = (1e-12)*ones(2*numpart,1);
% %svector(1:2,1) = (1e-2)*ones(numpart,1);
% svector(1:numpart,1) = (1e-2)*ones(numpart,1);
% %svector(3:4,1) = (1e-12)*ones(numpart,1);
stimlocindexold = stimlocindex;

% boundary conditions
%bc='d';
bc='n';
%bc='p';
if karma_test
    bc='p';
end
% if bc == 'n', Neumann (flux)
% if bc == 'd', Dirichlet
% if bc == 'p', periodic (ring geometry)
if bc == 'n'
    %    Vx_up = zeros(numstep,1);
    %    Vx_dn = zeros(numstep,1);
    Vx_up = zeros(writeintsteps,1);
    Vx_dn = zeros(writeintsteps,1);
elseif bc == 'd'
    %    Vup = zeros(numstep,1);
    Vup = zeros(writeintsteps,1);
    %Vup(2) = 2;
    Vup(1:400) = 2*ones(400,1);
    %Vup(2:6) = 2*ones(5,1);
    %    Vdn = zeros(numstep,1);
    Vdn = zeros(writeintsteps,1);
elseif bc == 'p'
    %    Vup = zeros(numstep,1);
    %    Vdn = zeros(numstep,1);
    Vup = zeros(writeintsteps,1);
    Vdn = zeros(writeintsteps,1);
end

dctr = zeros(numdet,1);%0; % rising edge counter
rctr = zeros(numdet,1);%0; % falling edge counter
apedges = cell(numdet,1);%[]; % start and end times of AP's
edgediff = cell(numdet,1);%[]; % diff of apedges
apstartindices = cell(numdet,1);%[]; % start indices of AP's
apendindices = cell(numdet,1);%[]; % end indices of AP's
apds = cell(numdet,1);%[]; % AP durations (units = number of timesteps)

% store previous and next values to carry across new looping structure
Vprev = zeros(numpart,1);
nprev = zeros(numpart,1);
Verrorprev = zeros(numpart,1);
nerrorprev = zeros(numpart,1);
usfprev = 0;
Vnext = zeros(numpart,1);
nnext = zeros(numpart,1);
Verrornext = zeros(numpart,1);
nerrornext = zeros(numpart,1);
usfnext= 0;
Vdnnext = 0; % only needed for periodic BC
Vupnext = 0;
%noiseamp = 0.5; 
noiseamp = 0; 
noiseprev = noiseamp*randn(numpart,1); 
% ff=figure;
% hold on;

tic
for ii=1:numrep % number of data writing cycles
    if ii > 1
        %        V = zeros(numpart,numstep);
        V = zeros(numpart,writeintsteps);
        if karma_test
            V(1,1) = 2; % initiate pulse that will travel around ring.
            % If we want to change the initiation point, then the periodic BCs must
            % change also. Right now, a pulse is started at the leftmost end, after
            % which the ends of the fiber are tied together: V_dn = V(end,...),
            % and Vup = Vdn
        end
        %        n = zeros(numpart,numstep);
        %        nerror = nd(:,1:numstep)-n;
        %        Verror = Vd(:,1:numstep)-V;
        n = zeros(numpart,writeintsteps);
    
        % new code to read in successive portions of Vd for use with
        % estimator
        Vd(sensorlocindices,:) = Vinterp([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)],sensorindices)'; 

        nerror = nd(:,1:writeintsteps)-n;
        Verror = Vd(:,1:writeintsteps)-V;
        % Reinstate initial values
        V(:,1) = Vnext;
        n(:,1) = nnext;
        Verror(:,1) = Verrornext;
        nerror(:,1) = nerrornext;
        if  bc == 'p'
            Vup(1) = Vupnext;
            Vdn(1) = Vdnnext;
        end
        usf = zeros(1,writeintsteps);
        if exist('icfile')
            iext(:,1) = iextinitnew(:,2);
            usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
        end
    end % if ii>1
    % iext=zeros(numpart,numstep); % external current
    % iint=zeros(numpart,numstep); % internal current
    % un=zeros(numpart,numstep); % optional feedback term added to n-dynamics
    % usf = zeros(1,numstep);
    iext=zeros(numpart,writeintsteps); % external current
    iint=zeros(numpart,writeintsteps); % internal current
    un=zeros(numpart,writeintsteps); % optional feedback term added to n-dynamics
%     %hprime = zeros(numpart,numstep);
%     %acoeff1 = zeros(numpart,numstep);
%     %acoeff2 = zeros(numpart,numstep);
%     %acoeff3 = zeros(numpart,numstep);
%     hprime = zeros(numpart,writeintsteps);
%     acoeff1 = zeros(numpart,writeintsteps);
%     acoeff2 = zeros(numpart,writeintsteps);
%     acoeff3 = zeros(numpart,writeintsteps);
%     %s=zeros(2,numstep);

    if ii==numrep
        wi = finalblocksteps;
    else
        wi = writeintsteps;
    end
    for k=1:wi %k = 1:numstep-1
        trueindex = (ii-1)*writeintsteps + k;
        if strcmp(fbon,'sf')
            if 0
                           %        hprime(:,k) = (1-(sech(Vd(:,k)-Vh)).^2).*((Vd(:,k).^2)/2)+(1-tanh(Vd(:,k)-Vh)).*Vd(:,k);
%                 %            hprime(:,k) = -((sech(Vd(:,k)-Vh)).^2).*((Vd(:,k).^2)/2)+(1-tanh(Vd(:,k)-Vh)).*Vd(:,k);
%%%%%%%% Errors upon errors in Rappel version. 
%See beginning of file for explanation
% There is an incomplete conversion to the new "h" that was used from
% initial change until 01/06/09. The "ch1,2,3" parameters were initially
% intended for use in the new tanh approximation to the power-series
% expansion in the Rappel version. 
% Original version, with untested corrections to hprime and acoeff2, using
% desired trajectories and relative time indices: 
%                             hprime(:,k) = ch1 + 2*ch2*Vd(:,k) + 3*ch3*Vd(:,k).^2; % Rappel version
%                             htemp = ch0 + ch1*Vd(:,k) + ch2*Vd(:,k).^2 + ch3*Vd(:,k).^3; 
%                             acoeff1(:,k) = (-1+(Vstar - (nd(:,k)/nB).^M).*hprime(:,k))/epsil2;
%% wrong (in UNYCES poster)                             acoeff2(:,k) = -0.5*(M*(nd(:,k)/nB).^(M-1)).*(1-tanh(Vd(:,k)-Vh)).*(Vd(:,k).^2)/epsil2;
% corrected                            acoeff2(:,k) = -(M*(nd(:,k)/nB).^(M-1)).*htemp/epsil2;
%                             acoeff3(:,k) = 0.5*(1/b)*(capprox*sech(capprox*(Vd(:,k)-Vn)).^2)/tauN;
%
% % Change to use true time indices
% %                             hprime(:,k) = ch1 + 2*ch2*Vd(:,trueindex) + 3*ch3*Vd(:,trueindex).^2; % Rappel version
% %                             acoeff1(:,k) = (-1+(Vstar - (nd(:,trueindex)/nB).^M).*hprime(:,k))/epsil2;
% % wrong (in UNYCES poster)                             acoeff2(:,k) = -0.5*(M*(nd(:,trueindex)/nB).^(M-1)).*(1-tanh(Vd(:,trueindex)-Vh)).*(Vd(:,trueindex).^2)/epsil2;
% %                           htemp = ch0 + ch1*Vd(:,trueindex) + ch2*Vd(:,trueindex).^2 + ch3*Vd(:,trueindex).^3; 
% % corrected                             acoeff2(:,k) =  -(M*(nd(:,trueindex)/nB).^(M-1)).*htemp/epsil2;
% %                             acoeff3(:,k) = 0.5*(1/b)*(capprox*sech(capprox*(Vd(:,trueindex)-Vn)).^2)/tauN;
% Rewrite to accommodate non-des trajectories (easier for Grammian tests)
                            hprime(:,k) = ch1 + 2*ch2*V(:,k) + 3*ch3*V(:,k).^2; % Rappel version
                            htemp = ch0 + ch1*V(:,k) + ch2*V(:,k).^2 + ch3*V(:,k).^3; 
% % check whether h's match
% outcheck = -V(:,k) + Vb + (Vstar - (n(:,k)/nB).^M).*htemp;
% outact = karma_f(V(:,k),n(:,k),nB,M);
% figure(ff)
% plot(trueindex,outcheck,'gs');
% plot(trueindex,outact,'bx');
% % The values appear to match, indicating that the h function (or at least the implied f function) used here is
% % the same as the one below
                            acoeff1(:,k) = (-1+(Vstar - (n(:,k)/nB).^M).*hprime(:,k))/epsil2;
% wrong (in UNYCES poster)                            acoeff2(:,k) = -0.5*(M*(n(:,k)/nB).^(M-1)).*(1-tanh(V(:,k)-Vh)).*(V(:,k).^2)/epsil2;
%                            acoeff2(:,k) = -(M*(n(:,k)/nB).^(M-1)).*htemp/epsil2;
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
                if numpart ==1
                    Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,k)+1];
                    A12 = [deltat*acoeff2(1,k)];
                    %        A21 = [0];
                    Auc = [1-deltat/tauN];
                    %%%%%%%%%
                    %%
                    % Alternative: use "true" A matrix
                    A21 = [deltat*acoeff3(1,k)];
                    A = [ Ac   A12
                        A21  Auc];
                end
                %            eiga(:,k) = eig(A);
                % %            eigata(:,k) = eig(A'*A);
                %  %           eigared(:,k) = eig(Ared);
                %             %         eigctrb1(:,k) = eig(ctrb(A,B1));
                %             %         eigctrb2(:,k) = eig(ctrb(A,B2));
                %             %         eigctrb1red(:,k) = eig(ctrb(Ared,B1red));
                %             %         eigctrb2red(:,k) = eig(ctrb(Ared,B2red));
                %             %         eigctrbired1(:,k) = eig(ctrb(AIred,BIred1));
                %             %         eigctrbired2(:,k) = eig(ctrb(AIred,BIred2));
                %             %         eigctrbi1(:,k) = eig(ctrb(AI,BI1));
                %             %         eigctrbi2(:,k) = eig(ctrb(AI,BI2));
            end
            if 0
                [AIredb,BIred1b,CIredb,MIredcomp,ranknum] = ctrbf(AIred,BIred1,[1 0 0 0],1e-7);
                eigctrbiredbar1(:,k) = eig(ctrb(AIredb(2:4,2:4),BIred1b(2:4))); % if using ctrbf
                eigaired(:,k) = eig(AIredb);
                if min(abs(eigctrbiredbar1(:,k))) > 1e-8
                    %    Kt=acker(AIredb(1:3,1:3),BIred1b(1:3),0.9*ones(1,3));
                    %    Kt=acker(AIredb(1:3,1:3),BIred1b(1:3),0.78*ones(1,3));
                    %    Kt=acker(AIredb(1:3,1:3),BIred1b(1:3),0.8*ones(1,3));
                    Kt=acker(AIredb(2:4,2:4),BIred1b(2:4),[0.9 0.9228 0.8838]); % if using ctrbf %***************
                    %            Kt=acker(AIredb(2:4,2:4),BIred1b(2:4),[0.9 0.87 0.8838]); % if using ctrbf
                    %            Kt=acker(AIredb(2:4,2:4),BIred1b(2:4),[0.9 1.1 1]); % if using ctrbf
                    %    [Kt,St,Et] = DLQR(AIredb(1:3,1:3),BIred1b(1:3),eye(3),1);
                    %    Kb=[Kt 1];
                    Kb=[1 Kt]; % if using ctrbf
                    KIred1=Kb*inv(MIredcomp); % This is wrong
                else
                    KIred1=zeros(1,4);
                end

                eigairedcl(:,k) = eig(AIredb-BIred1b*KIred1);
            end

            if 0%min(abs(eigctrb1(:,k))) > 1e-14 %& max(abs(eiga(:,k))) > 1.1
                %            K=acker(A,B1,[0.5; -0.5; 0.1; 0]);
                %            K1=acker(A,B1,[0.5; -0.5; 0.1; 0]);
                %            K=acker(A,B1,zeros(4,1));
                %                        K1=acker(A,B1,0.9*ones(4,1));
                %                        K1=acker(A,B1,0.999*ones(4,1));
                K1=acker(A,B1,[0.9991 0.9991 0.87 0.87]);
                %                        K1=acker(A,B1,0.95*eiga(:,k)');
            else
                K1=zeros(1,4);
            end
            if 0% min(abs(eigctrb2(:,k))) > 1e-8
                %            K=acker(A,B1,[0.5; -0.5; 0.1; 0]);
                K2=acker(A,B2,[0.5; -0.5; 0.1; 0]);
                %            K=acker(A,B1,zeros(4,1));
                %            K=acker(A,B1,0.9*ones(4,1));
            else
                K2=zeros(1,4);
            end

            if 0
                [Ab,B1b,Cb,Ma,ranknum] = ctrbf(A,B1,[1 0 0 0],1e-3); % gives avg 2.01, min 0, max 3 at smaller timestep.
                %        [Ab,B1b,Cb,Ma,ranknum] = ctrbf(A,B1,[1 0 0 0],5e-2); % too
                %        restrictive at smaller timestep

                % The following uses as many controllable states as possible
                numctrl = sum(ranknum);
                rankctr(:,k) = numctrl;
                if numctrl > 0
                    startac=size(A,2)-numctrl+1;
                    %       Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),0.9*ones(1,numctrl));
                    %       Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),0.9*ones(1,numctrl));
                    %if max(abs(eig(Ab(startac:end,startac:end)))) > 1
                    eab = eig(Ab(startac:end,startac:end));
                    Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),0.95*eig(Ab(startac:end,startac:end)));
                    %            Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),0.95*eig(Ab(startac:end,startac:end)));
                    %            Kbp =acker(Ab(startac:end,startac:end),B1(1:2),[eab(1) 0.98*eab(2)]);
                    %            Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),0.95*eig(Ab(startac:end,startac:end)));
                    %else
                    %       Kbp =acker(Ab(startac:end,startac:end),B1b(startac:end),1.05*eig(Ab(startac:end,startac:end)));
                    %end
                    Kbf=[zeros(1,size(A,2)-numctrl) Kbp];
                else
                    Kbf=zeros(1,4);
                end
                %         % This just assumes that the (approximate) n-subspace is always
                %         % uncontrollable, which is not always true
                %         Kbp =acker(Ab(3:end,3:end),B1b(3:end),[0.9 0.9]);
                %         Kbf=[zeros(1,2) Kbp];
                % %        Kred=Kbf*inv(Ma); % Wrong! The transformation matrix has the
                % %        opposite convention from 232 reader and progrep writeup. Replace
                % %        inv(M) (=M') with M.
                Kred=Kbf*Ma;

                %*               K1 = Kred;
                %                if V(stimlocindex,k) >=0 & V(stimlocindex,k) < 0.1
                %                if k > 1 & (V(stimlocindex,k) > V(stimlocindex,k-1))
                %                if k > 1 & (V(stimlocindex,k) < V(stimlocindex,k-1)) & V(stimlocindex,k) >=1.6 & V(stimlocindex,k) < 1.8 % for 0.91
                ea = eig(Ared);
                %                K1r = acker(Ared,B1red,0.95*eig(Ared)); % reducing to 0.9*eig(Ared) doesn't seem to help much. Way too harsh at smaller timestep.
                %%        K1r = acker(Ared,B1red,0.999*eig(Ared)); % for bcl091, with no windowing, wasn't able to find anything in the 0.99...*eig range that produced a desirable effect, either not much happens or the phase shift error starts in
                %        K1r = acker(Ared,B1red,[ea(1) 0.95]);
                %
                if 0 &1& numctrl == 2 & isreal(ea) %& k > 1 &(V(stimlocindex,k) < V(stimlocindex,k-1)) % adding directional constraints doesn't seem to help, probably because all the action appears to be happening on the downstroke
                    eas = sort(ea); % ascending order
                    %             %       K1r = acker(Ared,B1red,[ea(1) 0.95*ea(2)]); % this helps, but induces severe distortion of V1.
                    %             %
                    %             %        K1r = acker(Ared,B1red,[ea(1) 0.98*ea(2)]); %  works but induces
                    %             %        180deg phase shift. Largest reduction in alterrnans amplitude that I've seen so far.
                    %             %
                    %             %        K1r = acker(Ared,B1red,[ea(1) 0.99*ea(2)]);  %destabilizing?
                    %             %        K1r = acker(Ared,B1red,[ea(1) 0.999*ea(2)]);  % doesn't do much
                    %             %        K1r = acker(Ared,B1red,[ea(1) 0.9999*ea(2)]);  %doesn't do anything
                    %             %        K1r = acker(Ared,B1red,[eas(1) 0.98*eas(2)]); % destabilizing
                    %             %        K1r = acker(Ared,B1red,[eas(1) 0.999*eas(2)]); % ineffective
                    %             %        K1r = acker(Ared,B1red,[eas(1) 0.99*eas(2)]); % destabilizing
                    %             %        K1r = acker(Ared,B1red,[eas(1) 0.9*eas(2)]); % cancels entire APs
                    %             K1r = acker(Ab(startac:end,startac:end),B1red,[eab(1) 0.98*eab(2)]); % almost works? no phase shift? but I forgot the reverse transformation. Doesn't do anything useful at smaller timestep
                    %             K1r = acker(Ab(startac:end,startac:end),B1red,[eab(1) 0.99*eab(2)]); % almost works? no phase shift? but I forgot the reverse transformation. Doesn't do anything useful at smaller timestep
                    %             K1r = acker(Ab(startac:end,startac:end),B1red,[eab(1) 0.97*eab(2)]); % almost works? no phase shift? but I forgot the reverse transformation
                    %             K1r = acker(Ab(startac:end,startac:end),B1red,[eab(1) 0.96*eab(2)]); % Helps apds in 2nd cell, but induces 180deg phase error
                    %       K1r = acker(Ared,B1red,[0.98*eas(1) eas(2)]); % works OK, but not as well as the "wrong" ones above. Note: default deltat vs. 2e-4 give radically different results. The latter look very bad (control is too harsh)
                    %       K1r = acker(Ared,B1red,[0.985*eas(1) eas(2)]); % Loses control near the end. Works OK, but not as well as the "wrong" ones above
                    %       K1r = acker(Ared,B1red,[0.99*eas(1) eas(2)]); % Loses control near the end. Works OK, but not as well as the "wrong" ones above
                    %       K1r = acker(Ared,B1red,[0.99*ea(1) ea(2)]); % No effect at default deltat. At smaller deltat, it may speed up the destabilization or be the same as the open-loop case. Hard to tell.
                    %       K1r = acker(Ared,B1red,[eas(1) eas(1)]); % Destabilizing
                    %       K1r = acker(Ared,B1red,[eas(2) eas(2)]); % Hardly any effect
                    %       K1r = acker(Ared,B1red,[eas(1) (eas(1)+eas(2))/2]); % Destabilizing
                    %       K1r = acker(Ared,B1red,[eas(1) (0.01*eas(1)+0.99*eas(2))]); % No effect (at either time step)
                    %       K1r = acker(Ared,B1red,[eas(2) (eas(1)+eas(2))/2]); % No effect
                    %       K1r = acker(Ared,B1red,[eas(2) 0.95]); % Not much of an effect
                    %             % %            K1r = acker(Ared,B1red,[0.9 0.9]);
                    %            K1r=[0 0];
                else
                    K1r=[0 0];
                end

            end

            if 1
                %%%% Attempt at pseudo-single-cell control

                %        K1r1 = acker(Ared(1,1),B1red(1),0.95*eig(Ared(1,1))); % didn't
                %        work (with the vtest quantities)
                %            K1r1 = acker(Ared(1,1),B1red(1),0.95); % This sort-of works (with the vtest quantities), but not as well as performing single-cell control with the single-cell model and attaching the single cell to others, through diffusion. Why the difference?
                % Correction: At larger deltat, this works about as well as the quasi-single-cell
                % control, if (1) stimperiod = 0.91, and (2) the jon_desval_oncell...
                % quantities are used. Under the same conditions, leaving out the
                % diffusion term (as below) doesn't work as well.
                % The above method only works at the larger deltat. Gives bizarre results at the smaller deltat.
                %        K1r1 = acker(deltat*c1(1)+1,deltat*cc_ext,0.95); % Doesn't work (with the vtest quantities) (diffusion term has been removed, so that's not the problem).
                %            K1r1 = acker(Ared(1,1),B1red(1),1.01); %
                %            if k < 12/deltat
                %                K1r = [K1r1 0];
                %                            K1r = [-23 0];
                %            K1r = [-23 -23]; % This works also for the epsil1 0045, bcl
                %            087, M10 case -- gets rid of discordant alternans in both
                %            cells, but "less is more", the -10 -10 case achieves the same
                %            with less control activity and less distortion of the cell 1
                %            APs.
                %            K1r = [-10 -10];
                %            K1r = [0 0];
                %            K1r = [0 0 0];
                % epsil1 = 0.0053 results
                %            K1r = [-23 0 0]; % works
                %            K1r = [-10 -10 0]; % works
                %            K1r = [-7 -7 -7]; % doesn't work
                %            K1r = [-10 0 -10]; % doesn't work
                %            K1r = [-7 0 7]; % works OK
                %            K1r = [0 0 -7]; % kind of works for cells 2 and 3
                %            K1r = [0 0 7]; % doesn't work at all (2:1 is induced)
                %            K1r = [-1 -1 -1]; % doesn't work at all
                %            K1r = [-10 -1 -1]; % works, but not better than feedback based on first or first and 2nd cells.
                %            K1r = [0 -10 0]; % kind of works for cells 2 and 3, messes up cell 1
                K1r = zeros(1,numpart);
                % epsil1 = 0.009 results
                %            K1r(1) = -23; % works OK for 0.009, but unfortunately, the non-resampled
                %            Vd (for 0.0053), seems to give better results
                %            K1r(1:2) = [-10 -10]; %works OK for 0.009(better than 1-cell control with "correct" Vd for 0.009), but some distortion of cell 1
                %            APs is induced
                K1r(1:3) = [-5 -5 -5]; %works OK for 0.009(better than 1-cell control with "correct" Vd for 0.009), but some distortion of cell 1
                %            K1r(1:3) = [-2 -2 -2]; %works OK for 0.009(better than 1-cell control with "correct" Vd for 0.009), but some distortion of cell 1
                %            APs is induced
                %            K1r(1:3) = [-3 -7 -11]; %helps speed up convergence in cells 2 and 3, but spatial attenuation is still present, along with
                % more noticeable distortion of cell 1 APs.

                %            else
                %                K1r = [0 K1r1]; % take cell 1 gain and wrap it around cell 2's error: cancels APD error in cell 2, but cell 1 AP gets fairly distorted (started, cancelled, restarted), although not obviously out of phase.
                %            end
                %            K1 = [K1r 0 0];
                %%%% end   pseudo-single-cell control
            end
            %        kgainswr(:,k) = K1r';
            %        kgains(:,k) = acker(Ared,B1red,[eab(1) 0.99*eab(2)])';

            %        K1 = [K1r 0 0 0];
            K1 = [K1r zeros(1,numpart)];
            if trueindex > fbon_time %***k > fbon_time %& abs(V(stimlocindex,k)) < 4.5
                if trueindex == 1 %***k==1 %|| k < 2*fbon_time % this will not be executed anymore since fbon_time >=1
                    %               usf(k) = -Kred*[Verror(:,k)]; % would prefer to use this

                    %usf(k) = -Kred*[Verror(:,k); nerror(:,k)];
                    usf(k) = -K1*[Verror(:,k); nerror(:,k)];
                else
                    %***                sext=[Verror(:,k); Verror(:,k)-Verror(:,k-1)];
                    sext=[Verror(:,k); Verror(:,k)-Verrorprev];
                    %                          usf(k) = usf(k-1)+KIred1*sext;
                    %                                          KIred1 = K1;
                    %                          KIred1 = [zeros(1,numpart) K1r]; % Note: the first half of the gain vector (the part that multiplies Verror) is the integral gain, and the other half (that multiplies Verror_k-Verror_k-1) is the proportional gain
                    %                           KIred1 = [-1 0 0 K1r]; % I-action added to cell 1 only. Helps cell 1 to converge faster
                    %                           KIred1 = [0 -1 0 K1r]; % I-action added to cell 2 only. Helps convergence of cells 2 and 3, but massive distortion of cell 1 APs. Need better tuning and/or antiwindup
                    %                           KIred1 = [0 0 -1 K1r]; % I-action added to cell 3 only. System blowup (current saturates at lower value). Need better tuning and/or antiwindup
                    %                          KIred1 = [-1 -0.3 0 K1r]; % I-action added to cells 1 and 2 only. Gain on cell 2 is too small, no effect.
                    %                KIred1 = [-1 -0.7 0 K1r]; % I-action added to cells 1 and 2 only. Doesn't seem to be possible to easily improve convergence in cells 2/3 without distorting cell 1 APs, using this approach (0.3 is too small, 0.7 too large, I don't think 0.5 was that good either).
                    %                KIred1 = [-1 0 -0.5 K1r]; % I-action added to cells 1 and 3 only. 0.3 doesn't do anything. 0.5 induces lots of cell 1 AP distortion without helping much with cell 2/3 convergence.
                    %                KIred1 = [-1 0 -0.5 K1r]; % I-action added to cells 1 and 3 only. 0.3 doesn't do anything. 0.5 induces lots of cell 1 AP distortion without helping much with cell 2/3 convergence.
                    %                KIred1 = [-1 -0.5 -0.5 0 0 0]; % OK, seems to work better than PI cases, but still see distortion of cell 1 APs (especially early in the simulation), along with spatial attenuation of control effectiventess
                    %                KIred1 = [-0.5 -1 -0.5 0 0 0]; % too much cell 1 distortion. Doesn't actually improve convergence in cell 2 relative to -1 -0.5 -0.5 case
                    %                KIred1 = [-1 -0.5 -0.5 -0.5 -0.5 0]; % Looks very similar to the case where the P gains are zero.
                    %                KIred1 = [-1 -0.5 -0.5 -0.5 -0.5 -0.5]; % Looks very similar to the case where the P gains are zero.
                    %                KIred1 = [-1 -0.5 0.5 0 0 0]; % effects of I gains appear to cancel for cells 2 and 3.
                    %                KIred1 = [-1 0.5 0.5 0 0 0]; % destabilizing.
                    %                KIred1 = [-1 -1 -1 -2 -2 -2]; % OK, but cell 1 AP's are
                    %                distorted
                    %                KIred1 = [-2 -2 -2 -1 -1 -1]; %  Looks about the same as previous case.
                    %                KIred1 = [-1 -0.5 -0.5 -5 -5 -5]; % Converges somewhat faster than corresponding P-only case, but the steady-state APD for cell 1 is much different (cells 2 and 3 are not strongly affected).
                    %                KIred1 = zeros(1,6); % Converges somewhat faster than corresponding P-only case, but the steady-state APD for cell 1 is much different (cells 2 and 3 are not strongly affected).
                    %                KIred1 = [0 -0.5 -0.5 -5 -5 -5]; % blows up
                    %                KIred1 = [-0.1 -0.5 -0.5 -5 -5 -5]; % cell 1 APs are wacky
                    %                KIred1 = [-0.5 -0.5 -0.5 -5 -5 -5]; % cell 1 APs are wacky
                    %                KIred1 = [-1.5 -0.5 -0.5 -5 -5 -5]; % cell 1 APs are less weird
                    %                KIred1 = [-3 -0.5 -0.5 -5 -5 -5]; %As the size of the first entry is increased, the relative distortion of cell 1 APs appears to decrease, but with no apparent benefit in the convergence times of the other cells. In fact the convergence of the other cells slows down as the first entry is increased in size.
                    %                KIred1 = [-3 -1 -0.5 -5 -5 -5]; % Improves 2/3 convergence relative to previous case, but still not as fast as some other cases with more distorted APs, such as [-1 -.5 -.5 -5 -5 -5]
                    %                KIred1 = [-3 -1.5 -0.5 -5 -5 -5]; % cell 1 AP distortion reappears. Looks similar to [-1 -.5 -.5 -5 -5 -5] case, but convergence times are actually somewhat worse. So, simply tuning up the gains doesn't help. The ratio between KI1 and KI2 seems to affect whether or not the cell 1 APs become distorted (need to keep KI2 << KI1 to avoid distortion).

                    %                if ~isempty(find(stimindices==k))
                    %                    tempsi = stimindices(find(stimindices==k));
                    %                end
                    %                 if k >= tempsi && k < tempsi + 500
                    %                     KIred1 = [0 -0.1 0 -5 -5 -5]; %
                    %                 else
                    %                     KIred1 = [0 0 0 -5 -5 -5]; %
                    %                 end
                    %                     KIred1 = [0 0 0 -20 0 0]; %
                    %                     KIred1 = [0 -0.01 0 -5 -5 -5]; % Works OK with fac=0.999
                    %                     KIred1 = [0 -0.05 0 -5 -5 -5]; % With fac=0.99, doesn't work as well as 0.999 / [0 -0.01 0 -5 -5 -5] combination. So the decay rate may matter more than the value of the gain.
                    %                     KIred1 = [0 -0.03 0 -5 -5 -5]; % Keeping 0.999, this
                    %                     works better than the -0.01 case, but at the expense
                    %                     of making larger negative dents in the cell 1 APs
                    %                KIred1 = [0 0 -0.01 -5 -5 -5]; % Works OK with fac=0.999, a bit better than putting the I-gain on cell 2
                    %                     KIred1 = [0 -0.01 -0.01 -5 -5 -5]; % Works OK with fac=0.999, causes more distortion in V1 than just using either cell 2 or cell 3 alone for the I-gain
                    %                     KIred1 = [0 0 -0.005 -5 -5 -5]; % Works OK with fac=0.999, a bit better than putting the I-gain on cell 2
                    %                     KIred1 = [-0.01 0 -0.01 -5 -5 -5]; % Adding gains back to cell 1 appears to reduce undershoot of V1 somewhat, but worsens convergence of cell 1/2
                    KIred1 = zeros(1,2*numpart);
                    %                     KIred1(3:4)=place(Ared,B1(1:2),0.999*eig(Ared));
                    %                     KIred1(3)=place(Ared(1,1),B1(1),0.999*Ared(1,1));
                    %                     KIred1(3:4)=place(Ared,B1(1:2),0.995*eig(Ared));
                    %                     KIred1(3:4) = min([KIred1(3:4); zeros(1,2)]); % only retain negative gains? (This didn't fix the problem)
                    %                     KIred1(3:4) = Kstore; % use time-averaged gains from previous pole-placement run
                    %                     KIred1(3:4) = [-0.2834   -0.0998];
                    %                     KIred1(3:4) = [-1 0];
                    %                     kgains(:,k) = KIred1(3:4)';
                    %                     eigacl1(:,k) = eig(A-B1*[KIred1(3:4) 0 0]);
                    %                     eigatacl(:,k) = eig((A-B1*[KIred1(3:4) 0 0])'*(A-B1*[KIred1(3:4) 0 0]));
                    %                     eigaclred(:,k) = eig(Ared-B1(1:2)*KIred1(3:4));
                    %                     KIred1 = [0 0 0 -0.01 -5 -5 -5 0];
                    %                     KIred1 = [0 0 0 -0.01 0 -5 -5 -5 0 0];
                    %                     KIred1(1:4) = [0 0 0 -0.01]; % 4-cell SFI controller applied to 40-cell system with stimulus at cell 1
                    %                     KIred1((numpart+1):(numpart+4))=[-5 -5 -5 0]; % 4-cell SFI controller applied to 40-cell system with stimulus at cell 1
                    %                     KIred1(2:5) = [0 0 0 -0.01]; % 4-cell SFI controller
                    %                      applied to 40-cell system with stimulus at cell 2
                    %                      KIred1((numpart+2):(numpart+5))=[-5 -5 -5 0]; %
                    %                      4-cell SFI controller applied to 40-cell system with
                    %                      stimulus at cell 2
                    %                     KIred1(1:5) = [0 0 0 0 -0.01]; % *** default 5-cell setting
                    %                     KIred1(6) = [-0.01]; %
                    %                     KIred1(7) = [-0.01]; %
                    %                     KIred1(8) = [-0.01]; %
                    %                     KIred1(8) = [0.01]; %
                    %                     KIred1(8) = [-0.005]; %
                    %                      KIred1(4) = [-0.01]; %
                    %                      KIred1(7) = [-0.01]; %
                    %                     KIred1((numpart+1):(numpart+20))= -0.1; % too low
                    %                     KIred1((numpart+1):(numpart+20))= -0.5; % too high
                    %                     KIred1((numpart+1):(numpart+20))= -0.15; % still too high?
                    %                     KIred1((numpart+1):end)= -0.1; % too high
                    %                     KIred1((numpart+1):end)= -0.05; % destabilizing
                    %                     KIred1((numpart+1):end)= -0.01; % almost no effect
                    %                     KIred1(numpart+8)= -5; % Some stabilizing effect, but
                    %                     the control current (and perturbations to the cell 1
                    %                     potential) are way too high. The control current is
                    %                     almost always high and saturated. Note that this is a
                    %                     much more severe effect than that produced by putting
                    %                     similar P-gains on cells 1:3.
                    %                     KIred1(numpart+8)= -1; % Strangely, this appears to be destabilizing, although the control currents are less severe than in the -5 case.
                    %                     KIred1(numpart+8)= 1; % This is also destabilizing.
                    %                     KIred1(numpart+8)= -0.5; % Improves long-term behavior a little, but V in cell 8 is 180deg out of phase with Vd (phase tracking error is maximized).
                    %                     KIred1([numpart+1 numpart+8])=[-5 -0.5]; %
                    %                     KIred1(numpart+10)= -0.5; % In addition to the default, Kp = -5 -5 -5, Ki(5) = 0.01, doesn't seem to do much?
                    %                     KIred1(10) = [0.01]; % Destabilizes the system, when used with Kp = -5 -5 -5, Ki(5) = 0.01, and Kp(10)=-10
                    %                     KIred1(10) = [0.005]; % Maybe some improvement in convergence to cells 10-11, when used with Kp = -5 -5 -5, Ki(5) = 0.01, and Kp(10)=-10
                    % Traveling gain test:
                    % speed = 0.7; % speed (cm/time unit) with which additional feedback terms travel down the fiber
                    % startsweep = 2*fbon_time;
                    % if k == startsweep
                    % initoffset = 6; % cell index where extra term starts
                    % inittime = startsweep;
                    % locindex = initoffset;
                    %                     KIred1(numpart+locindex)= -0.5;
                    %                     KIred1(locindex) = [0.005];
                    % elseif k > startsweep
                    %     locindex = initoffset + round(speed*deltat*(k-inittime)/deltax);
                    %     if locindex > numpart
                    %         disp(num2str(locindex))
                    %         locindex = initoffset;
                    %         inittime = k;
                    %     end
                    %                     KIred1(numpart+locindex)= -0.5;
                    %                     KIred1(locindex) = [0.005];
                    % end
                    % Single cell reproductions
                    %                KIred1((numpart+1):(numpart+1))= -20;
                    %                KIred1((numpart+1):(numpart+1))= -1;
                    %                KIred1((numpart+1):(numpart+1))= -0.5; % Rappel
%                                    KIred1((numpart+1):(numpart+1))= -0.1; % Rappel
                    %                KIred1((numpart+1):(numpart+1))= -0.01; % Rappel (too low)
                    % The critical proportional gain is between -0.05 and -0.01
%                                    KIred1((numpart+1):(numpart+1))= -0.05; % Rappel; also works for 2 cells
%                                    KIred1((numpart+1):(numpart+1))= -0.03; % Rappel; sort of works for 1 cell but not clear if it converges completely
                    % Single-cell control with n-based feedback applied to V dynamics 
%                                    KIred1((numpart+1):(numpart+1))= -5; % -4 is too small and -10 is too large (phase-shift error in V) 
%                                     -5 can remove alternans but still produces phase shift error
                                    
                    % Two-cell reproductions
                    %                KIred1((numpart+1):(numpart+2))= [-5 -5];
                    %                KIred1((numpart+1):(numpart+2))= [-2 -2];
                    %                KIred1((numpart+1):(numpart+2))= [-1 -1];
                    %                KIred1((numpart+1):(numpart+2))= [-0.05 -0.05]; % Rappel; works for 2 cells
                    %                if k > fbon_time + 3*stimperiodsteps(1)
                    %                KIred1(2) = [-0.005]; % bad transient
                    %                KIred1(2) = [-0.001]; % Rappel; works for 2 cells but is a bit strong (undershoot in V's)
                    %                KIred1(2) = [-0.0005]; % Rappel; works for 2 cells but is a bit strong (undershoot in V's)
                    %                KIred1(2) = [-0.0001]; % Rappel; works for 2 cells (not much undershoot in V)
                    %                end
                    %                KIred1(2) = [-0.01]; % bad transient
                    %                KIred1(2) = [-0.005]; % bad transient
                    %                KIred1(2) = [-0.01]; % too small
                    %                KIred1(2) = [-0.003]; % bad transient
                    %                KIred1(2) = [-0.002]; % bad transient
                    % 10 cell version (Rappel)
                    %                KIred1((numpart+1):(numpart+1))= -0.1; % Rappel
                    %                KIred1((numpart+1):(numpart+1))= -0.05; % not much effect on alternans (may be destabilizing), big current blobs
                    %                KIred1((numpart+1):(numpart+2))= [-0.05 -0.05]; % converts to stable alt?
                    %                KIred1((numpart+1):(numpart+2))= [-0.1 -0.1]; % Stabilizes alternans (slowly)
                    %                KIred1((numpart+1):(numpart+1))= -0.2; % Stabilizes
                    %                alternans, but -0.1 -0.1 works a bit better
                    %                KIred1((numpart+1):(numpart+10))= -0.05*ones(1,10); % Stabilizes alt (may be too agressive)
                    %                KIred1(10) = [-0.0001]; % When added to the above (all p gains at -0.05), increases rate of convergence slightly.
                    % Discordant (190, start at 2275)
                    %                KIred1((numpart+1):(numpart+2))= [-0.1 -0.1]; % converts to discordant alt (prevents or delays 2:1); stabilizes cells 1:2
                    %                KIred1((numpart+1):(numpart+4))= [-0.1 -0.1 -0.1 -0.1]; %
                    %                converts to discordant alt (prevents or delays 2:1);
                    %                stabilizes cells 1:3? Large undershoot in V1.
                    %                KIred1(10) = [-0.0001]; KIred1((numpart+1):(numpart+2))= [-0.1 -0.1];
                    % The above combo may work OK, if you run it longer. Looks
                    % a little "calmer" in APs than p-only case. However there
                    % is some undershoot in V1.
                    %                KIred1((numpart+1):(numpart+10))= -0.05*ones(1,10); % Converts to nondecaying discordant alt, with large undershoots in V1,2
%                                  KIred1((numpart+1):(numpart+10))= -0.02*ones(1,10); % Converts to expanding discordant alt, which crunches to a lower level around 4000 then starts growing back
%                                  KIred1((numpart+1):(numpart+20))= -0.04*ones(1,20); % 
%                                  KIred1((numpart+1):(numpart+20))= -0.02*ones(1,20); %
%                                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); % 
%                                  KIred1((numpart+1):(numpart+40))= -0.01*ones(1,40); % 
%                                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); % 
%                                  KIred1((numpart+1):(numpart+80))= -0.02*ones(1,80); % Too strong (destabilizing) 
%                                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80); % Too weak (doesn't stabilize enough)
%                                  KIred1((numpart+1):(numpart+80))= -0.01*ones(1,80); % Too strong (slowly destabilizing)
%                                  KIred1((numpart+1):(numpart+80))= -0.007*ones(1,80); % Too strong (quickly destabilizing)
%                                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80); % Too strong (quickly destabilizing)
%                                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80); % Too strong (quickly destabilizing)
%                                  KIred1((numpart+1):(numpart+80))= -0.0055*ones(1,80); % Too weak (doesn't stabilize enough)
%                                  KIred1((numpart+1):(numpart+80))= -0.004*ones(1,80); % Too strong (slowly destabilizing)
%                                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80); % Too strong (slowly destabilizing)
%                                  KIred1((numpart+41):(numpart+80))= 0.003*ones(1,40); % 
%                                  KIred1((numpart+1):(numpart+40))= -0.005*ones(1,40); % 
%                                  KIred1((numpart+1):(numpart+80))= -0.002*ones(1,80); %
%                                 KIred1((numpart+1):(numpart+80))= 0.001*ones(1,80); % Increases alt amp but doesn't destabilize
%                                  KIred1(numpart+80)= 0.1; % 
%if trueindex < 8000/deltat
%                                  KIred1(numpart+80)= 0.1; % 
%end
% Somewhere between -0.005 and -0.006, the current becomes "unmoored" from
% the id train; the total current ends up looking like a wide sinusoid that
% doesn't track the id curve at all. Given that there isn't much room
% between the two gain values, it doesn't look hopeful for a uniform
% stabilizing solution, using Kp only, at least not one that is robust to
% gain shifts. In a weird twist, both -0.0055 and -0.005 appear to be too
% weak, but the former appears to be less stabilizing, so I suppose we
% should go the other way. But -0.004 gives results similar to -0.006 (both
% come unmoored from id). -0.003: doesn't cause unmooring, but can't tell
% if it's stabilizing. If it is, it takes much longer than 15000ms. Adding
% back Ki(40) = -0.0002 converts back to an unmoored (unstable) solution. 
% Removing Ki, putting the left half of Kp at -0.005 and the right at
% -0.003 is destabilizing, but in an uneven way. L -0.005 and R +0.003
% doesn't seem to do any better than -0.003 everywhere, but in contrast it
% appears to convert disco alt to concordant alt. Not sure what to make of
% that. I can't tell if -0.002 everywhere is better or worse than -0.003.
% Adding back Ki(40) = -0.0002 leaves things fairly ambiguous, but
% different looking than Kp -0.002 alone. Tuning Ki up to -0.0003 returns
% to a destabilizing solution. Trying positive gains everywhere (Kp=0.0001
% then 0.001) shows a trend toward enhancing the alternans (by increasing
% its amplitude), but interestingly it doesn't automatically destabilize
% the fiber. This may give some insight to the disco-to-concordant case
% seen using a mix of positive and negative gains? Niels may say that a
% time-varying approach is warranted (covert to concordant, then
% eliminate?), but this already looks like it would take longer than just
% using the default 40-cell controller. Question: can feedback from the
% most distal cell help to improve the performance of the Kp(1:40) = -0.02
% controller? Adding Kp(80) = + or - 0.02 doesn't seem to have much of an
% effect. Same for -0.05. -0.1 is large enough to cause an obvious
% modification to iext. The convergence time seems to be somewhat worse.
% +0.1 seems to induce faster convergence for proximal cells (by about
% 8000ms) but slower convergence for distal cells, compared to no feedback
% from cell 80. Not sure about the dividing line (would need to check
% contour plots or more clearly labeled APD plots). Turning the cell 80
% feedback (0.1) off after 8000ms doesn't seem to produce a clear benefit over
% the 40-cell-only controller. The proximal cells still seem to converge
% faster, but the 40-cell-only controller still seems to yield smaller alternans
% amplitudes, overall, even if some convergence times are longer. I'm not 
% sure why negative feedback from cell 80 would tend to prolong convergence
% times without offering any clear reduction in alternans amplitude. Should
% try running the +/- scenarios beyond 15000 ms to see whether there is any
% clear winner in the longer run. Running longer (to 18000ms) doesn't help.
% I haven't been able to beat the 40-cell SFI controller using + or -
% P-gains from cell 80, or by moving the I-gain away from cell 40, although
% I haven't tried too many variations on the latter. 


%                                  KIred1((numpart+1):(numpart+40))= -0.025*ones(1,40); % Not clearly better or worse than -0.02
%                                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); % 
%                    KIred1(numpart+1) = -0.1; % seems to be stabilzing the left end of the 80-cell fiber, but slowly
%                    KIred1(numpart+1) = -0.5; % seems to be stabilzing the left end of the 80-cell fiber, but slowly
%                    KIred1(numpart+1) = -1; % seems to be stabilzing the left end of the 80-cell fiber, but slowly
%                    KIred1(numpart+1) = -5; % seems to be stabilzing the left end of the 80-cell fiber, but slowly
%                                  KIred1(numpart + [1:4:40])= -4*0.02*ones(1,10); % Does not stabilize the 80-cell fiber in conjunction with Ki(40) = -0.0002
%                                  KIred1(numpart + [1:4:40])= -6*0.02*ones(1,10); % Appears to stabilize the 80-cell fiber in conjunction with Ki(40) = -0.0002
                    
                    %                KIred1(10) = [-0.0001]; % When added to the above (all p gains at-0.02) converts to nondecaying discordant alt (large undershoot in V1,2) .
                    %                KIred1(10) = [-0.0005]; % When added to the above (all p gains at-0.02) too strong (way too much current).
                    %                KIred1(10) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
%                                    KIred1(10) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
%                                    KIred1(20) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
%                                    KIred1(40) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
%                                    KIred1(80) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
                    %                KIred1(10) = [-0.0002]; % When added to the above (all p gains at-0.02) seems to work, although there's a bad transient
                    %                KIred1(10) = [-0.00015]; % When added to the above (all p gains at-0.02) seems to work (but a bit slow). No bad transient, but prolonged undershoots in at least V1 due to lack of fast convergence.
                    %                KIred1(10) = [-0.00017]; % When added to the above (all p gains at-0.02) seems to work (but a bit slow). No bad transient, but prolonged undershoots in at least V1 due to lack of fast convergence.
                    %                KIred1(10) = [-0.00019]; % When added to the above (all p gains at-0.02) seems to work (but a bit slow). Somewhat similar to 0.0002 case, but while the early apd spikes are reduced, convergence doesn't start until around 5000.
%                                    KIred1(40) = [-0.0003]; % When added to Kp(1:40)=-0.02, works better than -0.0002 on the 80-cell system
%                                    KIred1(40) = [-0.0004]; % Better than -0.0002 when added to Kp(1:40)=-0.02, for 100 cells
%                                    KIred1(40) = [-0.0006]; % Converges better (but transient is quite terrible) than -0.0002, -0.0004 when added to Kp(1:40)=-0.02, for 100 cells (deltat=0.2)
%                                    KIred1(40) = [-0.0005]; % Converges better (but there is an intermediately bad transient only revealed at deltat=0.05, not 0.2) than -0.0002, -0.0004 when added to Kp(1:40)=-0.02, for 100 cells 
%                                    KIred1(50) = [-0.0002]; % Doesn't work as well as same gain at cell 40 
%                                    KIred1(60) = [-0.0002]; % Doesn't work as well as same gain at cell 40 
                    % single-cell LTV perturbation test
                    %if trueindex == round(670/deltat) + 5 %k == round(670/deltat) + 5
                    %    KIred1 = [0 -0.1];
                    %else
                    %    KIred1 = [0 0];
                    %end
                    % Lyapunov check
                    %*                KIred1((numpart+1):(numpart+1))= -0.1; % Rappel
                    %                KIred1((numpart+1):(numpart+1))= -0.05; % Rappel
                    %                KIred1((numpart+1):(numpart+1))= -0.02; % Rappel
                    %                KIred1((numpart+1):(numpart+1))= -0.015; % Rappel
                    % Single-site feedback applied to longer fiber (100 cells) 
%                    KIred1((numpart+1):(numpart+1))= -5; % improves upon OL, but not better than default SF (Kp(1:40) = -0.02) 
%                    KIred1((numpart+1):(numpart+1))= -10; % too high (weird behavior)
%                    KIred1((numpart+1):(numpart+1))= -7; % very modest improvement over Kp(1) = -5
%                    KIred1((numpart+1):(numpart+1))= -8; % improves upon OL, but current looks weird 
%                    KIred1((numpart+20):(numpart+20))= -7; % try moving the sensor to another location (destabilizing) 
%                    KIred1((numpart+20):(numpart+20))= -1; % try moving the sensor to another location (not an improvement over Kp(1) = -7)
%                    KIred1((numpart+40):(numpart+40))= -1; % try moving the sensor to another location
%                    For 100 cells, Kp(1) = -7 gives the best performance out of the tested values (-5, -7, -8, -10) 
% Test output feedback (3 electrodes) 100 cells 
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.02*ones(1,3); % May be helping, but not a lot 
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.04*ones(1,3); % slight improvement over 0.02
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -1*ones(1,3); % better convergence than -0.04, but transient is moderately bad                     
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.8*ones(1,3); % very similar to -1, transient still somewhat bad
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.6*ones(1,3); % very similar to -1, transient still somewhat bad
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.4*ones(1,3); % This one seems workable at deltat=0.2. Try adding Ki. Moderately bad transient for deltat=0.05
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.3*ones(1,3); % This one seems workable at deltat=0.2. Try adding Ki. Moderately bad transient for deltat=0.05 (in conjunction with Ki(40)=-0.0004
%                    KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.25*ones(1,3); % This one is designated as the "best" output-feedback-integral result (with Ki(40)=-0.0004) at deltat=0.05, that doesn't produce a weird transient. However, it may not be the "best" standalone output feedback result, since a) I didn't do a careful search for the best 3-element Kp at the smaller timestep, and b) adding Ki has a very large effect on the convergence time (and presumably the transient), so it may be possible to increase Kp further if Ki=0. 
%                    KIred1(40) = [-0.0004]; % Improves on Kp-only of -0.4, for 100 cells
%                    KIred1([(numpart+1) (numpart+20) (numpart+40)])= -1.05*0.25*ones(1,3); 
%                    KIred1(40) = 1.05*[-0.0004]; 
%                     KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.95*0.25*ones(1,3); 
%                     KIred1(40) = 0.95*[-0.0004]; 
%                    KIred1([(numpart+1) (numpart+20) (numpart+40)])= -0.25*ones(1,3)/5; % 
%                    KIred1(40) = [-0.0004]/5; % Improves on Kp-only of -0.4, for 100 cells
%                     KIred1(40) = [-0.0006]; % Not significantly better than -0.0004
% Discordant alt (bcl=185, 100 cells)
%                    KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); % by itself, controls temporarily. Raising/lowering didn't seem to improve
%                    KIred1((numpart+1):(numpart+20))= -0.02*ones(1,20); % untethered? 
%                    KIred1((numpart+1):(numpart+20))= -0.01*ones(1,20); % untethered? 
%                    KIred1((numpart+1):(numpart+20))= -0.002*ones(1,20); % remains 2:1
%                    KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); % untethered
%                    KIred1((numpart+1):(numpart+1))= -5; % I think it sort of works on cell 1 but everything else (including current) looks terrible
%                    KIred1((numpart+1):(numpart+1))= 0.1; % remains 2:1
%                    KIred1((numpart+1):(numpart+1))= -0.1; % weird, but not 2:1
%                    KIred1((numpart+1):(numpart+1))= -1; % some proximal oscillate wildly, distal stay at 2:1
%                    KIred1((numpart+1):(numpart+1))= -10; % some proximal oscillate wildly, some distal stay at 2:1, but the current looks so bad, I wouldn't increase the gain further
                    %                    KIred1(40) = [-0.0002]; % tried -0.0001, -0.0004 
% Discordant alt (bcl=210, 100 cells)
%                    KIred1((numpart+1):(numpart+20))= -0.01*ones(1,20); % proximal cells do something, distal quickly return to 2:1 (deltat 0.2)
%                    KIred1((numpart+1):(numpart+1))= -5; % distal end remains blocked. cell 1 apds look OK, but there is an oscillation that appears to widen with distance from the controller (deltat 0.2)
%                    KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); % restores discordant solution, but eventually returns to 2:1 (deltat 0.2)
%*                    KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); % delays return to 2:1 slighly over -0.02 (deltat 0.2)
%                    KIred1((numpart+1):(numpart+60))= -0.02*ones(1,60); % worse than Kp(1:40) = -0.02 or -0.04 (current becomes untethered) (deltat 0.2)
%                    KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40); % delays return to 2:1 slighly over -0.04 alone. In combo with Ki(40) = -0.001, the transient is worse than -0.04, not sure about long term (deltat 0.2)
%                    KIred1(40) = [-0.0002]; % In combination with Kp(1:40) = -0.02, delays 2:1 over Kp(1:40) = -0.04 alone. Didn't retest at smaller deltat. 
%                    KIred1(40) = [-0.0004]; % In combination with Kp(1:40) = -0.02, outperforms -0.0002 (does not become unmoored at smaller deltat)
%                    KIred1(40) = [-0.0006]; % In combination with Kp(1:40) = -0.02, outperforms -0.0004 (seems to remain discordant, but with an expanding envelope, up to 30000ms). Didn't retest at smaller deltat. 
%                    KIred1(40) = [-0.0008]; % In combination with Kp(1:40) = -0.02, outperforms -0.0006. Didn't retest at smaller deltat. 
%                    KIred1(40) = [-0.001]; % In combination with Kp(1:40) = -0.02, outperforms -0.0008, but current is unmoored at smaller deltat
%                    KIred1(40) = [-0.002]; % In combination with Kp(1:40) = -0.02, looks very weird (not in a good way); high current. Didn't retest at smaller deltat. 
%                    KIred1(40) = [-0.0015]; % In combination with Kp(1:40) = -0.02, outperforms -0.001, but current is high (unmoored at smaller deltat)
%                    KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); % For bcl=225, in combination with Ki(40) = -0.0004, actually seems to be better (convergence and energy) for 100 cell case than Kp(1:40) = -0.02
%                    KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); % For bcl=225, in combination with Ki(40) = -0.0004, converges faster for 100 cell case than Kp(1:40) = -0.02, but uses more energy
%                    KIred1(40) = [-0.0003]; % For bcl=225, in combination with Kp(1:40) = -0.02, actually seems to be better for 100 cell case than -0.0004 
%                    KIred1(40) = [-0.0002]; % 
%                    KIred1(40) = [-0.0006]; % 
%                    KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);
%                    KIred1((numpart+1):(numpart+50))= -0.03*ones(1,50); % Not helpful. Current becomes unmoored in conjunction with Ki(40) = -0.0003
%                    KIred1((numpart+1):(numpart+50))= -0.02*ones(1,50); % Not helpful. Current becomes unmoored in conjunction with Ki(40) = -0.0003
%                    KIred1((numpart+1):(numpart+50))= -0.01*ones(1,50); % Not helpful. Current becomes unmoored, when used by itself. 
%                    KIred1((numpart+1):(numpart+50))= -0.005*ones(1,50); % Probably not helpful. 2:1 block appears to occur sooner, when used by itself. 
%                    KIred1((numpart+1):(numpart+45))= -0.005*ones(1,45); % Probably not helpful. 2:1 block appears to occur sooner, when used by itself. 
%*                    KIred1(40) = [-0.0003]; % -0.02 (1:40) combined with -0.0003 or -0.0004 at 40 (which is fine for 225) doesn't work (becomes unmoored) if initiated earlier (2800 vs. 10000) in 100 cells, 210ms. Adjusting fac to 0.9999 doesn't seem to help; increasing to 0.99999 makes the system unmoored
% Some of the above refer to a 10000ms start time? For 2800 ms, Kp(1:40) =-0.04, Ki(40) = -0.0003, fac=0.9999 appears to stabilize very slowly. Unsat is noticeably better than sat. 
%165, 20 cells
%                    KIred1((numpart+1):(numpart+10))= -0.01*ones(1,10); %  
%                    KIred1((numpart+1):(numpart+20))= -0.005*ones(1,20); % 165, 20 cells, OK for a bit then 2:1 resumes 
%                    KIred1((numpart+1):(numpart+20))= -0.007*ones(1,20); % wild oscillations in cells 1 and 2, doesn't break block in other cells
%                    KIred1((numpart+1):(numpart+20))= -0.004*ones(1,20); % not quite as good as -.005
%                    KIred1((numpart+1):(numpart+1))= -2*ones(1,1); % only works to control cells 1 and 2 (others have expanding oscillations or block)
%                  KIred1((numpart+1):(numpart+2))= -1*ones(1,2); % mildly better than -2 applied to cell 1 only
%                  KIred1((numpart+1):(numpart+2))= -2*ones(1,2); %mildly better than -1 applied to cells 1:2 only
%                  KIred1((numpart+1):(numpart+2))= -4*ones(1,2); %not a significant improvement over -2 applied to cells 1:2 only
%                  KIred1((numpart+1):(numpart+3))= -2*ones(1,3); %possibly a mild improvement in cells 3:4 over -2 applied to cells 1:2 only
%                  KIred1((numpart+1):(numpart+4))= -2*ones(1,4); %an improvement over -2 applied to cells 1:3, but current is very large
%                  KIred1((numpart+1):(numpart+4))= -1*ones(1,4); %slower-acting than applying -1 to cells 1:4, but reduces current somewhat
%                  KIred1((numpart+1):(numpart+5))= -1*ones(1,5); %similar to -2 applied to cells 1:4
%                  KIred1((numpart+1):(numpart+7))= -0.7*ones(1,7); %transfers more cells from block to large oscillation compared with -1 applied to 1:5. Current is very large though. 
%                  KIred1((numpart+1):(numpart+10))= -0.5*ones(1,10); %somewhat better than -1 applied to 1:7.
%                  KIred1((numpart+1):(numpart+12))= -0.4*ones(1,12); %transfers more cells from block to large oscillation compared with -0.5 applied to 1:10. Current is very large though.
%                  KIred1((numpart+1):(numpart+14))= -0.3*ones(1,14); %transfers more cells from block to large oscillation compared with -0.4 applied to 1:12. Current is very large though.
%                  KIred1((numpart+1):(numpart+16))= -0.2*ones(1,16); %transfers more cells from block to large oscillation compared with -0.3 applied to 1:14. Current is very large though.
%                  KIred1((numpart+1):(numpart+18))= -0.1*ones(1,18); %transfers more cells from block to large oscillation compared with -0.2 applied to 1:16. Current is very large though.
%*                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); %transfers more cells from block to large oscillation compared with -0.1 applied to 1:18, although transfer is inconsistent. Current is very large though.
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); %Combined with Ki(20) = -0.001 for 20-cell fiber, sort of works but phase shift is becoming noticeable (not 180 deg though) as gains are reduced. Current is very large though. Current is also arguably unmoored, although it contains periodic inward spikes (just time-shifted ahead of id spikes)
%                  KIred1((numpart+1):(numpart+20))= -0.02*ones(1,20); %Combined with Ki(20) = -0.001 for 20-cell fiber, sort of works but phase shift is becoming noticeable (not 180 deg though) as gains are reduced. Current is very large though. Current is also arguably unmoored, although it contains periodic inward spikes (just time-shifted ahead of id spikes)
%                  KIred1((numpart+1):(numpart+20))= -0.2*ones(1,20); %May be somewhat worse (higher current and fewer cells transferred) than -0.1 applied to 1:20
                  %                    KIred1(10) = [-0.0004];
%                    KIred1(5) = [-0.0001]; % at fac=0.999, no major difference over -1 feedback over cells 1:5 alone
%                    KIred1(5) = [-0.0003]; % at fac=0.999, no major difference over -1 feedback over cells 1:5 alone
%                    KIred1(5) = [-0.001]; % at fac=0.999, no major difference over -1 feedback over cells 1:5 alone
%                    KIred1(20) = [-0.0001]; % wrt Kp=-0.1 on 1:20, with fac=0.999, helps somewhat but doesn't really converge. Very similar to the fac=0.9999 case with -0.001, but with lower current. Best phase-tracking of Vd's I've seen, though, even though it doesn't converge. 
%*                   KIred1(20) = [-0.001]; % combined with -0.1 applied to cells 1:20, eliminates alternans at the expense of large current. The APs are overdriven in proximal cells but still recognizable. There is noticeable phase shift error but not as bad as for lower P-gains such as -0.05 and -0.02 on 1:20
%                    KIred1(20) = [-0.001]; % combined with -0.1 applied to cells 1:20, with fac=0.9999 instead of 0.999, seems to induce somewhat worse phase error and overdriving compared to 0.999. Current is also arguably unmoored, although it contains periodic inward spikes (just time-shifted ahead of id spikes)
%                    KIred1(10) = [-0.001]; % combined with -0.1 applied to cells 1:20, very similar to KIred1(20) = -0.001, just takes longer to converge. 
%                    KIred1(20) = [-0.0005]; % combined with -0.1 applied to cells 1:20, with fac=0.9999 instead of 0.999, seems to induce somewhat worse phase error and overdriving compared to 0.999. Current is also arguably unmoored, although it contains periodic inward spikes (just time-shifted ahead of id spikes). Basically similar behavior to -0.001, just with slower convergence. 
%                    KIred1(20) = [-0.0001]; % combined with -0.1 applied to cells 1:20, with fac=0.9999 instead of 0.999, is not an improvement over -0.0005. It doesn't really converge or recover from the unmooring.  Hence, it seems 0.9999 may only be helpful under certain conditions (longer fibers?) but not others. 
%165, 40 cells
% 20-cell controller (Kp(1:20) = -0.1, Ki(20) = -0.001, fac=0.999) doesn't pull all cells out of block.
%                  KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30); %Superior to Kp(1:20) = -0.1 with Ki(20)=-0.001 at fac=0.999, but ultimately doesn't stabilize system (apds crunch then expand).
%                  KIred1((numpart+1):(numpart+40))= -0.1*ones(1,40); %Inferior to Kp(1:30) = -0.1, with Ki(20)=-0.001 at fac=0.999, but ultimately doesn't stabilize system (apds crunch then expand repeatedly).
%*                  KIred1((numpart+1):(numpart+30))= -0.2*ones(1,30); %Superior to Kp(1:30) = -0.1, with Ki(20)=-0.001 at fac=0.999. May confer limited stability to the system, but there re repeated contractions and expansions. 
%                  KIred1((numpart+1):(numpart+30))= -0.25*ones(1,30); %Inferior to Kp(1:30) = -0.2, with Ki(20)=-0.001 at fac=0.999, but ultimately doesn't stabilize system (apds crunch then expand). 
%                  KIred1(40) = [-0.001]; % combined with -0.1 applied to cells 1:20, with fac=0.999, does not improve over Ki(20) = -0.001
%                  KIred1(20) = [-0.001]; % Partially effective in combination with Kp(1:30) = -0.2. fac=0.9999 is inferior (clearly not stable) compared with 0.999
%*                  KIred1(20) = [-0.002]; % combined with -0.2 applied to cells 1:30, with fac=0.999, has similar performance to Ki(20) = -0.001. Envelope of discordant alternans may be somewhat smaller than that of -0.001.
%                  KIred1(30) = [-0.001]; % combined with -0.2 applied to cells 1:30, with fac=0.999, does not improve over Ki(20) = -0.001. Expands and contracts, but expansions are much worse than the Ki(20) = -0.001 or -0.002 cases. 
%165, 40 cells, revised profile (at dapd=119) 
% Kp(1:30) = -0.2, Ki(20) = -0.002, fac=0.999 performs much better than with original desired profiles, and
% succesfully converts to 1:1 pattern. 
%                  KIred1((numpart+1):(numpart+30))= -0.2*ones(1,30); 
%*                  KIred1((numpart+1):(numpart+30))= -0.2*ones(1,30); 
%                  KIred1((numpart+1):(numpart+20))= -0.25*ones(1,20); % Not strong enough combined with Ki(20) = -0.001
%                  KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30); % too weak in combination with Ki(20) = -0.002, fac=0.999. Stabilizes very slowly. Distal end still shows noticeable alternans. 
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); % too weak
%                  KIred1(30) = [-0.002]; % with Kp(1:30) =  -0.2, and fac=0.999, system explodes
%                  KIred1(30) = [-0.001]; % with Kp(1:30) =  -0.1, and fac=0.999 (similar pattern to 20-cell controller), slowly loses control but uses less current than Kp(1:30) =  -0.2, Ki(20) = -0.002. With Kp(1:30) = -0.2, does an even worse job
%                 KIred1(20) = [-0.001]; % with Kp(1:30) = -0.2, and fac = 0.999, stabilizes very slowly
%*                KIred1(20) = [-0.002]; % with Kp(1:20) = -0.25, and fac = 0.999, system blows up. Works in combo with Kp(1:30)=-0.2

%165, 60 cells
%30-cell controller (Kp(1:30) = -0.2, Ki(20) = -0.002, fac=0.999), lapses in and out of block, then blows membrane potentials up to -infinity
%*                  KIred1((numpart+1):(numpart+40))= -0.1*ones(1,40); %By itself, better than to Kp(1:30) = -0.2, with Ki(20)=-0.002 at fac=0.999, in the sense that the system doesn't blow up, but ultimately doesn't stabilize system (apds crunch then expand in an erratic fashion).
%                  KIred1(30) = [-0.001]; % when combined with Kp(1:40) = -0.1, with fac=0.999, causes blowup
%*                  KIred1(30) = [-0.0005]; % when combined with Kp(1:40) = -0.1, with fac=0.999, outperforms Kp(1:40) = -0.1 by itself (less variability), but the system still lapses in and out of block
%                  KIred1(40) = [-0.0005]; % when combined with Kp(1:40) = -0.1, with fac=0.999, is not clearly better than Ki(30) = -0.0005. The system still lapses in and out of block, and may have more lapses than the previous case. 
%                  KIred1(30) = [-0.00025]; KIred1(40) = [-0.00025]; % when combined with Kp(1:40) = -0.1, with fac=0.999, looks about the same as Ki(40) = -0.0005. 
%                  KIred1(30) = [-0.0002]; % when combined with Kp(1:40) = -0.1, with fac=0.999, is not clearly better than Ki(30) = -0.0005. The system still lapses in and out of block, and may have more lapses than the previous case. 
%165, 60 cells, revised profile
%                  KIred1((numpart+1):(numpart+40))= -0.1*ones(1,40); % Don't know if this is better than the Kp(1:30) = -0.2, Ki(20) = -0.0015, fac=0.999 combo 
%                  KIred1((numpart+1):(numpart+40))= -0.15*ones(1,40); % In combination with Ki(30) = -0.0005, fac=0.999, system blows up. Also blows up in combination with Ki(20) = -0.0015
%                  KIred1(30) = [-0.0005]; % when combined with Kp(1:40) = -0.1, with fac=0.999 and new desired profile, still doesn't stabilize, but the intervals between lapses are longer? 
%                  KIred1(30) = [-0.0007]; % when combined with Kp(1:40) = -0.1, with fac=0.999, seems to increase frequency of lapses relative to -0.0005
%*                  KIred1((numpart+1):(numpart+30))= -0.2*ones(1,30); 
%                  KIred1((numpart+1):(numpart+30))= -0.22*ones(1,30); % May be slightly worse than -0.2, not sure
%                  KIred1(20) = [-0.002]; % with Kp(1:30) = -0.2 and fac = 0.999, holds off block for a while then loses its grip (return to block)
%                  KIred1(20) = [-0.003]; % with Kp(1:30) = -0.2 and fac = 0.999, system blows up 
%                  KIred1(20) = [-0.0025]; % with Kp(1:30) = -0.2 and fac = 0.999, works less well than -0.002
%*                  KIred1(20) = [-0.0015]; % with Kp(1:30) = -0.2 and fac = 0.999, hard to tell if better or worse than -0.002
%165, 60 cells, dapd=115
%Doesn't seem like a good set of desired profiles, overall. 
%                  KIred1((numpart+1):(numpart+30))= -0.2*ones(1,30); 
%                  KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30); % too weak by itself
%                  KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30); 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); 
%165, 60 cells, dapd=100
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.0001]; % OK but doesn't stabilize whole fiber 
%                  KIred1((numpart+1):(numpart+20))= -0.15*ones(1,20);  KIred1(20) = [-0.0001]; % Too high
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(20) = [-0.0001]; % OK but lets go after a while 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(20) = [-0.0005]; % OK but lets go after a while 
%*                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0005]; % works
%165, 80 cells, dapd=100
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0005]; % holds for a bit then lets go
%                  KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30) = [-0.0005]; % worse than Kp(1:30)=-0.07, Ki(30) = -0.0005
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0007]; % Not sure of comparison to Kp(1:30)=-0.07, Ki(30) = -0.0005. The minimum amplitude is higher, but the blowup happens somewhat later? 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.001]; % Seems worse than both Kp(1:30)=-0.07, Ki(30) = -0.0005 or -0.0007. Lapses are more frequent. 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0005]; % with fac=0.9999, blows up 
%                  KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30);  KIred1(30) = [-0.0004]; % with fac=0.9999, not clear if it's any better or worse than previous mixed cases, but has more frequent lapses. 
%                  KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30);  KIred1(30) = [-0.0006]; % with fac=0.9999, holds longer than foregoing one, but still lapses
%                  KIred1((numpart+1):(numpart+20))= -0.02*ones(1,20);  KIred1(30) = [-0.0006]; % with fac=0.9999, clearly worse than foregoing
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20);  KIred1(30) = [-0.0006]; % with fac=0.9999, blows up
%                  KIred1((numpart+1):(numpart+40))= -0.003*ones(1,40);  KIred1(30) = [-0.0006]; % with fac=0.999, holds for a while but with complete phase inversion of APs
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002];% with fac=0.999, lapses. Not better than Kp(1:30)=-0.07, Ki(30) = -0.0005
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [-0.0002];% with fac=0.999, very similar to Kp(1:30)=-0.07, Ki(30) = -0.0007, but maybe slightly better
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [-0.0002];  KIred1(60) = [-0.0002];% with fac=0.999, getting worse again (lapsing)
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [-0.0002];  KIred1(60) = [0.0002];% with fac=0.999, very similar to Kp(1:30)=-0.07, Ki(30) = -0.0007, but not as good as Ki([30 40 50])
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [0.0002];% with fac=0.999, not an improvement over Kp(1:30)=-0.07, Ki(30) = -0.0007
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(20) = [-0.0002]; KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [-0.0002];% with fac=0.999, not an improvement over Kp(1:30)=-0.07, Ki(30) = -0.0007
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(20) = [0.0002]; KIred1(30) = [-0.0002]; KIred1(40) = [-0.0002]; KIred1(50) = [-0.0002];% with fac=0.999, not an improvement over Kp(1:30)=-0.07, Ki(30) = -0.0007. Positive I-gains aren't helping so far
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:35) = [-0.0002];% with fac=0.999, not an improvement over Kp(1:30)=-0.07, Ki(30) = -0.0007. Current looks high, untethered and result lapses in and out of block
%*                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0002];% with fac=0.999, this may be the best one so far, at least through 10000ms. Not stable, but appears to delay blowup. With fac=0.9995, appears to transition to block sooner. With fac=0.995, transitions to block much sooner
%                  KIred1((numpart+1):(numpart+30))= -0.08*ones(1,30);  KIred1(30:32) = [-0.0002];% with fac=0.999, increasing P-gains makes it worse. 
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30);  KIred1(30:32) = [-0.0002];% with fac=0.999, decreasing P-gains also makes it worse. With fac=0.9999, seems to induce block sooner than OL. 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(30:32) = [-0.0002];% with fac=0.999, not bad? Better than Kp(1:30)=-0.07, Ki(30) = -0.0007. but not as good as Kp(1:30)=-0.07, Ki(30:32) = -0.0007 (doesn't delay block for as long?) 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(30:33) = [-0.0002];% with fac=0.999, goes to block sooner than Kp(1:30)=-0.07, Ki(30) = -0.0007 or Kp(1:30)=-0.07, Ki(30:32) = -0.0007 or Kp(1:80) = -0.005, Ki(30:32) = -0.0002
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.0005];% with fac=0.999, blows up around 7000ms
%*                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.999, maybe not quite as good as Kp(1:30)=-0.07, Ki(30) = -0.0007. With fac=0.9995, get a weird 4:4 pattern, where some APs are in phase and many aren't; envelope seems to be expanding though. What I thought was 4:4 is actually partially blocking (hence n:m); the average difference between repolarizations is somewhere between 165 and 2*165ms. fac=0.9997 looks worse.
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000008];% with fac=0.999, worse than Kp(1:80)=-0.005, Ki(1:80) = -0.000005 
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.999, maybe slightly worse than Kp(1:80)=-0.005, Ki(1:80) = -0.000005 and slightly better than Kp(1:80)=-0.005, Ki(1:80) = -0.000008
%                  KIred1((numpart+1):(numpart+80))= -0.008*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.999, maybe worse than Kp(1:80)=-0.006 (transitions sooner?) 
%                  KIred1((numpart+1):(numpart+80))= -0.008*ones(1,80);  % not stabilizing by itself; transitions around 6000ms 
%                  KIred1((numpart+1):(numpart+80))= -0.01*ones(1,80);  % not stabilizing by itself; worse than -0.008; transitions more sharply before 6000ms?  
%                  KIred1((numpart+1):(numpart+80))= -0.015*ones(1,80);  % not stabilizing by itself; worse than -0.01; transitions more sharply around 5500ms?  
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  % This is the same "amount" of gain distributed over a smaller area. Not stabilizing by itself; transitions in and out of block in a seemingly random fashion. Arguably worse (or at least less predictable) than Kp(1:80) =  -0.015
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  % not stabilizing by itself; the transition still starts around 6000ms, but is less complete; the full transition appears to be around 7500ms
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  % APDs rise and fall; not sure if system remains out of block past 10000ms
%                  KIred1((numpart+1):(numpart+80))= -0.001*ones(1,80);  % no apparent effect on 2:1 block (it is in block the whole time)
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(30:32) = [-0.0002];% With fac=0.999, adding I-gain induces more periods of approximate 1:1, but still lapses. The lapsed periods have oscillation, with somewhat reduced APD, though.
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(30:32) = [-0.0003];%With fac=0.999, not sure. Increased I-gain seems to keep it out of things that look like they could lead to block longer, but amplitudes of some of the lapses are worse. Changing fac=0.9995 hardly seems like an improvement; seems more random / less stable.
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); % a brief period without block
%                  KIred1((numpart+1):(numpart+30))= -0.09*ones(1,30); % no periods of convergence; some cells do large oscillations instead of block
%                  KIred1((numpart+1):(numpart+30))= -0.12*ones(1,30); % no periods of convergence; some cells do large oscillations instead of block (one additional cell is pulled into large oscillations compared with -0.09)
%                  KIred1((numpart+1):(numpart+30))= -0.17*ones(1,30); % no periods of convergence; some cells do large oscillations instead of block (one additional cell is pulled into large oscillations compared with -0.12, but at least one of the upstream cells returns to block at a reduced amplitude?)
%                  KIred1((numpart+1):(numpart+30))= -0.22*ones(1,30); % no periods of convergence; some cells do large oscillations instead of block (no additional cells are pulled into large oscillations compared with -0.17)
%                  KIred1((numpart+1):(numpart+35))= -0.12*ones(1,35); % Not bad; there is a crunch by all cells down to 1:1 but envelope is expanding by 10000ms. Comparable to KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0002];% with fac=0.999, at least out to 10000ms
%                  KIred1((numpart+1):(numpart+40))= -0.12*ones(1,40); % Worse than Kp(1:35) = -0.12; escape to block starts around 7000ms
%                  KIred1((numpart+1):(numpart+40))= -0.1*ones(1,40); % Almost the same as Kp(1:40) = -0.12; escape to block starts around 7000ms
%                  KIred1((numpart+1):(numpart+40))= -0.08*ones(1,40); % Almost the same as Kp(1:40) = -0.1; escape to block starts around 7000ms
%                  KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40); % Lapsing in and out of block
%                  KIred1((numpart+1):(numpart+35))= -0.15*ones(1,35); % Transition to block happens sooner than Kp(1:35) = -0.12 (maybe around 9500ms)
%                  KIred1((numpart+1):(numpart+35))= -0.1*ones(1,35); % Transition to block happens much sooner than with Kp(1:35) = -0.12 (around 5500ms)
%                  KIred1((numpart+1):(numpart+35))= -0.2*ones(1,35); % Transition to block happens sooner than Kp(1:35) = -0.12 (around 9000ms) 
%                  KIred1((numpart+1):(numpart+35))= -0.2*ones(1,35); % Transition to block happens sooner than Kp(1:35) = -0.12 (around 9000ms) 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000007];% with fac=0.9995, generally looks worse than the Ki(1:80) = -0.000005 variant, although unlike the former, this one has a brief period of 1:1 before lapsing into block 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.00001];% with fac=0.9995, generally looks worse than the Ki(1:80) = -0.000005 variant, although unlike the former, this one has a brief period of 1:1 before lapsing into block. Not sure of comparison to Ki(1:80) = -0.000005. The oscillations seem larger in this one, post-lapse
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000003];% with fac=0.9995, generally looks worse than the Ki(1:80) = -0.000005 and -0.000007 variants, although unlike the former, this one has a brief period of 1:1 before lapsing into block. The oscillations may be larger than those of the -0.000007 and -0.000001 cases. 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995, generally looks worse than the Ki(1:80) = -0.000005 and -0.000007 variants, although unlike the former, this one has a brief period of 1:1 before lapsing into block. The oscillations may be larger than those of the -0.000007 and -0.000001 cases. 
%                  KIred1((numpart+1):(numpart+20))=[ -0.0573   -0.0818   -0.0452   -0.0201   -0.0906   -0.0571   -0.0891   -0.0037   -0.0895   -0.0325   -0.0676   -0.0239   -0.0502   -0.0291   -0.0380   -0.0779   -0.0520   -0.0965   -0.0612   -0.0970]; KIred1(30) = [-0.0004];
                  % -rand(1,20)/10; 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30)*(1+sin(4*pi*trueindex*deltat/165))/2;  KIred1(30:32) = [-0.0002]*(1+sin(4*pi*trueindex*deltat/165))/2;% with fac=0.999, no benefit
%                  KIred1((numpart+1):(numpart+30))= 0.05*ones(1,30); KIred1(1:80) = [-0.000005];% Positive gain (with or without negative Ki) just seems to reinforce (or at least not disrupt) 2:1 block 
%                  KIred1((numpart+1):(numpart+30))= 0.05*ones(1,30); KIred1(1:80) = [0.000005];% produces sub-threshold APs (below -100) 
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(1:80) = [0.000005];% Looks similar to many negative-Kp / negative-Ki runs; brief convergence to 1:1 followed by lapse back into 2:1 block
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(1:80) = -[0.000005];% converges sooner than positive-Ki run; brief convergence to 1:1 followed by lapse back into 2:1 block. The transition to block is also sooner, so the total time spent in the initial crunch period looks about the same
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0002]; %with sext2=[Verror(:,k)-noisenow; sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))]; %with fac=0.999, appears to produce a wandering solution (intermediate block but not as bad as 2:1?) 
%KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30:32) = [-0.0002]; % with sext2 as above, very similar (maybe slightly smaller amplitude) to the above
%KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30); % with sext2 as above, almost no effect on 2:1 block. Hence any benefit is due to the I-gains.
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0003]; %with sext2=[Verror(:,k)-noisenow; sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))]; %with fac=0.999, appears to produce a wandering solution (intermediate block but not as bad as 2:1?). May keep the average amplitude lower than the -0.0002 case. 1:1 is temporarily achieved but is wildly out of phase. 
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.00005]; %with sext2=[sign(Verror(:,k)-noisenow); sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))], fac=0.999, has no effect on 2:1 block
%                  KIred1((numpart+1):(numpart+80))= -0.0044*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) intermediate blocking pattern until 10000ms, followed by brief 1:1, followed by erratic behavior
%*                  KIred1((numpart+1):(numpart+80))= -0.004*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), stays in an intermediate blocking pattern 
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), spends some time in an intermediate blocking pattern that is different from the previous one; erratic after 13000ms 
%                  KIred1((numpart+1):(numpart+80))= -0.004*ones(1,80);  KIred1(1:80) = [-0.000005]; KIred1(30:32) = [-0.0002];% with fac=0.999, adding extra I-gains is not an improvement over base case. fac=0.9995 seems to be worse, and fac=0.995 seems much worse
%                  KIred1((numpart+1):(numpart+80))= -0.004*ones(1,80);  KIred1(1:80) = [-0.000005]; KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); % with fac=0.9995, not an improvment over Kp(1:30) = -0.07, Ki(30:32) = -0.0002. 
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000005]; KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); % with fac=0.9995, not an improvment over Kp(1:30) = -0.07, Ki(30:32) = -0.0002. -0.006 gain is better than -0.004, though
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000008]; KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);% with fac=0.9995, not an improvment over Kp(1:30) = -0.07, Ki(30:32) = -0.0002. -0.000008 is worse than -0.000005
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000005]; KIred1((numpart+1):(numpart+30))= -0.01*ones(1,30); % with fac=0.9995, just looks like temporary tracking, quickly returning to block 
%                  KIred1((numpart+1):(numpart+80))= -0.008*ones(1,80);  KIred1(1:80) = [-0.000005]; % with fac=0.9995, expanding alternans as of 7000ms
%                  KIred1((numpart+1):(numpart+80))= -0.01*ones(1,80);  KIred1(1:80) = [-0.000005]; % with fac=0.9995, expanding alternans as of 7000ms. Possibly a very slight improvement over -0.008
%                  KIred1((numpart+1):(numpart+80))= -0.02*ones(1,80);  KIred1(1:80) = [-0.000005]; % with fac=0.9995, expanding alternans as of 7000ms. Seems worse than -0.01
%                  KIred1((numpart+1):(numpart+80))= -0.01*ones(1,80);  KIred1(1:80) = [-0.000007]; % with fac=0.9995, expanding alternans as of 7000ms. Seems worse than -0.01
%165, 100 cells, dapd=100 profile  
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0002]; % with fac=0.999, returns to 2:1 block starting around 8000ms. W/o sat, explodes immediately after being turned on (at 2800ms) 
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0004]; % with fac=0.999, shows erratic lapsing behavior throughout. W/o saturation, blows up around 4000ms
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0002]; % with fac=0.999, W/o sat, explodes around 3800ms 
%KIred1((numpart+1):(numpart+30))= -0.04*ones(1,30);  KIred1(30) = [-0.0002]; % with fac=0.999, W/o sat, explodes around 6800ms
%KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30);  KIred1(30) = [-0.0002]; % with fac=0.999, W/o sat, lapses in and out of block. With fac=0.9999, current wanders more, lapses look more erratic
%KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30);  KIred1(30) = [-0.0004]; % with fac=0.999, W/o sat, many cells remain in 2:1 block, but others seem to wander aimlessly and slowly below the 2:1 APD value. No cells remain in 1:1 for any significant length of time. 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.999 (w/sat), appears to return to block around 7500ms
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), is similar to 80-cell case (some kind of intermediate blocking pattern, non-convergent). More erratic after 11000ms
%*                  KIred1((numpart+1):(numpart+90))= -0.0044*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) appears to prevent lapse although it doesn't stabilize
%165, 120 cells, dapd=120 profile  
%KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0004]; % with fac=0.999, shows erratic lapsing behavior throughout. W/o saturation, blows up around 4000ms
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000007];% with fac=0.9995 (w/sat) starts lapsing around 6000ms
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) is better than -0.000007 case
%*                  KIred1((numpart+1):(numpart+90))= -0.0044*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) appears to prevent lapse although it doesn't stabilize
%165, 140 cells, dapd=140 profile  
%*                  KIred1((numpart+1):(numpart+90))= -0.0044*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) appears to prevent lapse although it doesn't stabilize
%165, 160 cells, dapd=160 profile  
%*                  KIred1((numpart+1):(numpart+90))= -0.0044*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) appears to prevent lapse although it doesn't stabilize

% if trueindex > 200000
%     KIred1=zeros(1,2*numpart); 
%                   KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30:32) = [-0.0002];%with fac=0.9995 merging this with Kp(1:80) = -0.005, Ki(1:80)=-0.000005, fac=0.999 looks like the two individual results concatenated; no stabilization 
% end
%180, 20 cells, dapd=125 profile or dapd=110 profile 
% The 20-cell 165ms controller works (fac=0.999 unless otherwise indicated)  
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20);
%                  KIred1(20) = [-0.001];
%180, 40 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.001]; % stabilizes system slowly 
%                  KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30); KIred1(20) = [-0.001]; % worse than above 
%*                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; % works
%180, 60 cells,  dapd=110 profile 
%*                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; % stabilizes slowly 
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); KIred1(30) = [-0.001]; % too strong
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); KIred1(30) = [-0.0005]; % too strong
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.0005]; % may stabilize slowly
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(30) = [-0.0005]; % unstable
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0005]; % more unstable
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.001]; % an improvement over Ki(40) = -0.0005 but still unstable
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.001]; %  still unstable
%                  KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40); KIred1(40) = [-0.001]; % may stabilize but current is very large 
%180, 80 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; % too weak? Lapses in and out of block for 0.999, blows up for 0.9999
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.0005]; % not sure
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.001]; % holds for a while then lets go
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); KIred1(30) = [-0.001]; % seems clearly worse than the above
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(40) = [-0.001]; % seems clearly worse than Kp(1:30)=-0.05, Ki(30) = -0.001
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.002]; % very similar to the -0.001 case
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(40) = [-0.0005]; % inferior to the Ki(30) =-0.001 and -0.002 cases
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); % doesn't appear to be stabilizing by itself
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); % doesn't appear to be stabilizing by itself. Maybe better than Kp(1:30)=-0.05
%                  KIred1((numpart+1):(numpart+40))= -0.06*ones(1,40); % worse than -0.04
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(50) = [-0.0005]; % adding an I-gain is not an improvement
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(40) = [-0.0001]; % adding an I-gain is not an improvement
%                   KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(40) = [-0.00005]; % destabilizing
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(30) = [-0.00005]; % destabilizing
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(30) = [-0.00005]; % destabilizing
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(60) = [-0.00005]; % destabilizing? but not as bad as the same I-gain at 30 
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(60) = [-0.00005]; % destabilizing
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); % doesn't really seem to be stabilizing anything
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.002]; % lapses in and out of block
%                  KIred1((numpart+1):(numpart+20))= -0.2*ones(1,20); KIred1(30) = [-0.002]; % causes explosion
%                  KIred1((numpart+1):(numpart+20))= -0.13*ones(1,20); KIred1(30) = [-0.002]; % possibly worse than -0.1. Not sure
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.002]; %  lapses in and out of block
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.003]; %  lapses in and out of block. Non-blocking periods have somewhat lower amplitude than for -0.002, but the current is more untethered
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.0003]; %  lapses in and out of block. (fac = 0.9999)
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.0006]; %  lapses in and out of block. (fac = 0.9999)
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(35) = [-0.0006]; %  worse than placing at 30. (fac = 0.9999)
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.0006]; %  not clearly better or worse than Kp(1:20) = -0.05. (fac = 0.9999)
%                  KIred1((numpart+1):(numpart+30))= -0.03*ones(1,30); KIred1(30) = [-0.0006]; %  not clearly better or worse than previous (fac = 0.9999)
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30);  KIred1(30) = [-0.0005]; % with fac=0.9995, lapses in and out of block
%                  KIred1((numpart+1):(numpart+80))= -0.004*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), 1:1 for a while, but with slowly increasing envelope
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), 1:1 for a while, but starts lapsing around 7000ms
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope by 7500ms. APs for all similar gain values are almost completely out of phase with desired values. 
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(1:80) = [-0.000007];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Current is larger amplitude than -0.000006 case, and envelope may be increasing faster
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Current is larger, and rate of envelope expansion may be slower than -0.005 case. Expands rapidly toward end of run and starts lapsing around 13000ms
%                  KIred1((numpart+1):(numpart+80))= -0.008*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Current is larger, and rate of envelope expansion may be slower than -0.005 case. Expands rapidly toward end of run and starts lapsing around 11700ms
%                  KIred1((numpart+1):(numpart+80))= -0.006*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Expands rapidly toward end of run and starts lapsing around 11300ms, but lapse looks slightly less severe than previous case
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Erratic after 8400ms
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), 1:1 for a while, but shows signs of expanding envelope. Erratic after 8400ms
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), erratic behavior starting around 6000ms, but almost settles into intermediate block around 10000ms
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000008];% with fac=0.9995 (w/sat), erratic behavior starting around 7000ms, although periods of lapse may be outweighed by partial block 
%                  KIred1((numpart+1):(numpart+80))= -0.002*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), erratic behavior between 7000ms and 11000ms, and starting again just before 15000ms. No clear signs of partial block
%                  KIred1((numpart+1):(numpart+40))= -0.006*ones(1,40);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), quite different compared with spreading P-gain over all 80 cells. Brief lapse around 3000ms, and from 11000ms onward
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9997 (w/sat), increasing fac from 0.9995 to 0.9997 seems to disrupt intermediate block more than it encourages it. Lapsing is dominant between 6000ms and 10000ms, and after 13000ms. With fac=0.9993, discordant alternans but with lapsing between 8000 and 11000ms. With fac=0.999, lapsing and erratic after about 10000ms                  
%*                  KIred1((numpart+1):(numpart+80))= -0.001*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. Transient is pretty bad
%                  KIred1((numpart+1):(numpart+80))= -0.0009*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), very similar to -0.001 case, although transition to intermediate block happens a bit later
%                  KIred1((numpart+1):(numpart+80))= -0.0008*ones(1,80);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), very similar to -0.0009 case, although transition to intermediate block happens a bit later
%180, 80 cells,  dapd=105 profile 
%                  KIred1((numpart+1):(numpart+30))= -0.02*ones(1,30); KIred1(30) = [-0.0006]; %  may be worse than previous (fac = 0.9999). Performance definitely is not improved with dapd=120 trajectories. May be better with 105 trajectories, but current looks weird. Changing back to 0.999 makes current look less weird but tracking doesn't last as long. 
%                  KIred1((numpart+1):(numpart+30))= -0.05*ones(1,30); KIred1(30) = [-0.002]; % still doesn't hold for 105 profile, fac=0.999
%                  KIred1((numpart+1):(numpart+30))= -0.07*ones(1,30); KIred1(30) = [-0.002]; % worse than -0.05 for 105 profile, fac=0.999
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.003]; % may be slightly better for dapd=105 than 110; expanding envelope vs lapsing? fac=0.999
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.004]; % Current is untethered but amplitudes are smaller. Not sure this could be considered to be successful tracking. fac=0.999, dapd=105
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.006]; % Blows up at the end. fac=0.999, dapd=105
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(30) = [-0.0006]; % Probably not the best and not clearly improved by use of dapd=105. Lapses in and out of block (fac = 0.9999)
%180, 100 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+80))= -0.003*ones(1,80);  KIred1(1:80) = [-0.000006];% with fac=0.9995 (w/sat), erratic behavior between 7500 and 10000ms. There is some intermediate block but it doesn't look very regular or stable
%                  KIred1((numpart+1):(numpart+90))= -0.0044*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) erratic after 10000ms
%                  KIred1((numpart+1):(numpart+90))= -0.004*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat) erratic after 6000ms
%                  KIred1((numpart+1):(numpart+90))= -0.003*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), some large-amp partial block, but erratic after 13000ms
%                  KIred1((numpart+1):(numpart+90))= -0.002*ones(1,90);  KIred1(1:80) = [-0.000005];% with fac=0.9995 (w/sat), is mostly erratic throughout
%                  KIred1((numpart+1):(numpart+90))= -0.003*ones(1,90);  KIred1(1:80) = [-0.000008];% with fac=0.9995 (w/sat), is mostly erratic throughout. Increasing I-gain has seemed to make it more erratic. 
%                  KIred1((numpart+1):(numpart+90))= -0.003*ones(1,90);  KIred1(1:80) = [-0.000002];% with fac=0.9995 (w/sat), is mostly erratic throughout. Decreasing I-gain has ambiguous effects: current has smaller amplitude, behavior is still mostly erratic throughout, although lapses might be somewhat less frequent
%                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000002];% with fac=0.9995 (w/sat), this yields an apparently stable intermediate blocking pattern, although the amplitude is larger than comparable results for 165ms
%*                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. The amplitude of the oscillation may be smaller than for -0.000002, but there is a large lapse between 6000 and 8000ms, and the inter-AP spacing may be larger
%                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000003];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. The amplitude of the oscillation may be smaller than for -0.000002. The lapse period is more contained than for -0.000004
%180, 120 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. The amplitude of the oscillation may be smaller than for -0.000002, but there is a large lapse between 6000 and 8000ms, and the inter-AP spacing may be larger
%180, 140 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. The amplitude of the oscillation may be smaller than for -0.000002, but there is a large lapse between 6000 and 8000ms, and the inter-AP spacing may be larger
%180, 160 cells,  dapd=110 profile 
%                  KIred1((numpart+1):(numpart+90))= -0.001*ones(1,90);  KIred1(1:80) = [-0.000004];% with fac=0.9995 (w/sat), this yields a possibly stable intermediate blocking pattern. The amplitude of the oscillation may be smaller than for -0.000002, but there is a large lapse between 6000 and 8000ms, and the inter-AP spacing may be larger

%195, 20 or 40 cells,  dapd~=115 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.001]; % works
%195, 60 cells,  dapd~=115 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.001]; % may stabilize slowly? 
%*                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; % stabilizes slowly 
%195, 80 cells,  dapd~=115 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; % delays 2:1 but doesn't stabilize
%                    KIred1((numpart+1):(numpart+30))= -0.04*ones(1,30);  % doesn't seem to stabilize by itself. Not sure if it's better or worse than Kp(1:20) = -0.1
%                    KIred1((numpart+1):(numpart+30))= -0.04*ones(1,30);  KIred1(30) = [-0.001]; %May hold at discordant alternans. Not sure? 
%                    KIred1((numpart+1):(numpart+30))= -0.04*ones(1,30);  KIred1(40) = [-0.001]; %Worse than putting Ki(30) = -0.001
%                    KIred1((numpart+1):(numpart+30))= -0.04*ones(1,30);  KIred1(40) = [-0.0002]; %Worse than putting Ki(40) = -0.001
%                    KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40);  KIred1(40) = [-0.0002]; %OK but not sure it completely stabilizes
%                    KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  KIred1(40) = [-0.0002]; %OK but not sure it completely stabilizes. Early amplitudes look larger than for Kp(1:40) = -0.02, but long term behavior is better
%*                    KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  KIred1(40) = [-0.0004]; %Appears to stabilize very, very slowly. Works better with fac=0.9999 than 0.999
%195, 100 cells,  dapd~=115 profile 
%                    KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  KIred1(40) = [-0.0004]; %With fac=0.9999 (the 80-cell controller), does not appear to come close to stabilizing. May lapse in and out of block. Haven't tried 0.999. With 0.9995, still spends most of its time in block
%                    KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40);  KIred1(40) = [-0.0004]; %With fac=0.9995, looks almost the same as previous, but with higher current
%                    KIred1((numpart+1):(numpart+60))= -0.0015*ones(1,60);  KIred1(40) = [-0.0004]; %With fac=0.9995, appears to go through expansion and contraction cycles
%                    KIred1((numpart+1):(numpart+60))= -0.0015*ones(1,60);  KIred1(1:60) = [-0.000007]; %With fac=0.9995, horrible transient between 4000 and 7000ms, small amplitude followed by expanding amplitude by 15000ms. It's not actually tracking during the small-amplitude range (180deg phase error) 
%                    KIred1((numpart+1):(numpart+60))= -0.0015*ones(1,60);  KIred1(1:80) = [-0.000006]; %With fac=0.9995, lapses in and out of block
%                    KIred1((numpart+1):(numpart+60))= -0.0015*ones(1,60);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, lapses in and out of block. Similar to -0.000006
%*                    KIred1((numpart+1):(numpart+80))= -0.0005*ones(1,80);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, appears to stabilize about an intermediate blocking rhythm 
%                    KIred1((numpart+1):(numpart+80))= -0.001*ones(1,80);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, lapses in and out of block 
%*                    KIred1((numpart+1):(numpart+90))= -0.0005*ones(1,90);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, also appears to stabilize about an intermediate blocking rhythm 
%                    KIred1((numpart+1):(numpart+90))= -0.0005*ones(1,90);  KIred1(1:90) = [-0.000004]; %With fac=0.9995, looks like lapsing up through 7500; not sure if it will settle down into anything 
% Could try modifying range of I-gains and/or increasing/decreasing fac
%195, 120 cells,  dapd~=115 profile 
%                    KIred1((numpart+1):(numpart+80))= -0.0005*ones(1,80);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, appears to stabilize about an intermediate blocking rhythm 
%195, 140 cells,  dapd~=115 profile 
%                    KIred1((numpart+1):(numpart+80))= -0.0005*ones(1,80);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, appears to stabilize about an intermediate blocking rhythm 
%195, 160 cells,  dapd~=115 profile 
%                    KIred1((numpart+1):(numpart+80))= -0.0005*ones(1,80);  KIred1(1:80) = [-0.000004]; %With fac=0.9995, appears to stabilize about an intermediate blocking rhythm 
%210, 20 cells,  dapd~=125 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.001]; % Works with fac=0.999, although current is probably higher than it needs to be. Could try reducing the gains.
%210, 40 cells,  dapd~=125 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20);
%                  KIred1(20) = [-0.001]; % Works with fac=0.999. Stabilizes slowly. 
%210, 60 cells,  dapd~=125 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(20) = [-0.001]; % With fac=0.999, may stabilize slowly. 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; %  With fac=0.999, stabilizes slowly 
%210, 80 cells,  dapd~=125 profile 
%                  KIred1((numpart+1):(numpart+20))= -0.1*ones(1,20); KIred1(30) = [-0.001]; %  With fac=0.999, stabilizes slowly. Using fac=0.9999 only yields a very minor improvement 
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  KIred1(40) = [-0.0004]; % With fac=0.999, stabilizes slowly 
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40);  KIred1(40) = [-0.0005]; % With fac=0.999, noticeably worse than Ki(40) = -0.0004
%*                  KIred1((numpart+1):(numpart+40))= -0.025*ones(1,40); KIred1(40) = [-0.0004]; % With fac=0.999, may be slightly better than Kp(1:20) = -0.1, Ki(30) = -0.001. The main difference is that the longer-range controller seems to disturb proximal cells more but provides a more even convergence. Could try reducing Kp's further. Turning control on later (12000ms) shows the controller pulls the system out of block and stabilizes (slowly). Transient is bad though. 
%210, 100 cells
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40);  KIred1(40) = [-0.0003]; % With fac=0.9999, may slowly stabilize system, but not very well 
%210, 120 cells
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40);  KIred1(40) = [-0.0003]; % With fac=0.9999, discordant alt that goes to 2:1 block around 10000ms
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40);  KIred1(1:40) = [-0.000008]; % With fac=0.9999, discordant alt that goes to 2:1 block around 6000ms
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40);  KIred1(40:42) = [-0.0002]; % With fac=0.9999, discordant alt, probably not stable 
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, discordant alt, probably not stable. Looks identical to previous up to 7500 
%                  KIred1((numpart+1):(numpart+40))= -0.06*ones(1,40);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, discordant alt, probably less stable than previous
%                  KIred1((numpart+1):(numpart+80))= -0.001*ones(1,80);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, looks like intermediate blocking pattern
%                  KIred1((numpart+1):(numpart+80))= -0.01*ones(1,80);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, discordant alt, haven't run long enough to determine stability, APs out of phase
%*                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, looks like intermediate blocking pattern. With fac=0.9995, lapses in and out of block
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(40:59) = [-0.0000035]; % With fac=0.9999, lapses in and out of block, APs out of phase, initially similar to -0.01, -0.00006 case
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(40:59) = [-0.0000035]; % With fac=0.9999, discordant alt, haven't run long enough to determine stability, APs out of phase
%210, 140 cells
%*                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, looks like intermediate blocking pattern. With fac=0.9995, lapses in and out of block
%210, 160 cells
%                  KIred1((numpart+1):(numpart+80))= -0.005*ones(1,80);  KIred1(40:49) = [-0.00006]; % With fac=0.9999, looks like intermediate blocking pattern. With fac=0.9995, lapses in and out of block
%225, 20 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999 or 0.9999 (the latter is more aggressive), although current is probably higher than it needs to be. Could try reducing the gains.
%225, 40 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999 (didn't try 0.9999), although current is probably higher than it needs to be. Could try reducing the gains. dapd looks too low
%225, 60 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999 (didn't try 0.9999), although current is probably higher than it needs to be. Could try reducing the gains. dapd looks too low, still
%225, 80 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % With with fac=0.999 (didn't try 0.9999), may stabilize slowly. dapd looks too low, still
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0004];  % % With with fac=0.999 (didn't try 0.9999), works
%225, 120 cells
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0004];  % % With with fac=0.999 (didn't try 0.9999), may stabilize slowly
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.0004];  % % With with fac=0.999 (didn't try 0.9999), may stabilize slowly. Performance almost looks same as previous.
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.0006];  % % With with fac=0.999 (didn't try 0.9999), stabilizes slowly. Performance is noticeably better than that with Ki(40) = -0.0004.
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.0008];  % % With with fac=0.999 (didn't try 0.9999), stabilizes slowly. Proximal profiles are quite distorted but phase tracking looks OK.
%225, 140 cells
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.0008];  % % With with fac=0.999, may stabilize slowly. There is some phase lag in tracking, and the current is close to untethered. With fac=0.9999, the current is probably too high (untethered), with maybe 90deg phase shift in tracking. The level of deformation in V's may be tolerable though. 
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(50) = [-0.0008];   % With with fac=0.999, much worse than with Ki(40) = -0.0008. 
%                  KIred1((numpart+1):(numpart+40))= -0.03*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, looks very similar to Ki(40) = -0.0008 (but probably has higher current). 
%                  KIred1((numpart+1):(numpart+50))= -0.025*ones(1,50); KIred1(40) = [-0.001];   % With with fac=0.999, looks worse than previous.  
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, may be better in terms of alternans amplitude than Kp(1:40) = -0.03, Ki(40) = -0.0008 or -0.001, but current is still pretty much untethered, and phase tracking may be somewhat worse (not sure why) 
%                  KIred1((numpart+1):(numpart+40))= -0.01*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, has similar alternans amplitude to Kp(1:40) = -0.03, Ki(40) = -0.0008 or -0.001, but current is completely untethered, and phase tracking error is approaching 180deg 
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, this is the best so far, hence it's better to increase P-gains rather than reduce them (although spreading them over a larger region didn't seem to help). Current is not completely untethered, and phase tracking error is present but not as bad as previous
%*                  KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, appears to stabilize slowly. Although amplitude is somewhat larger than in previous case, the current is noticeably smaller and better-anchored than in previous case
%225, 160 cells
%                  KIred1((numpart+1):(numpart+40))= -0.04*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, may confer "neutral" stability. Current seems large but is not completely untethered
%*                  KIred1((numpart+1):(numpart+40))= -0.05*ones(1,40); KIred1(40) = [-0.001];   % With with fac=0.999, appears to stabilize slowly. Current is noticeably smaller and better-anchored than in previous case
%240, 20 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999, although current is probably higher than it needs to be. Could try reducing the gains. dapd is too low. 
%240, 40 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999, although current is probably higher than it needs to be. Could try reducing the gains. dapd is too low. 
%240, 60 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999, although current is probably higher than it needs to be. Could try reducing the gains. dapd is too low. 
%240, 80 cells
%                  KIred1((numpart+1):(numpart+20))= -0.05*ones(1,20); KIred1(20) = [-0.0005]; % Works with fac=0.999, although current is probably higher than it needs to be. Could try reducing the gains. dapd is too low. 
%240, 100 cells
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0003];  % With fac=0.999 (didn't try 0.9999), stabilizes. dapd is too low.
%240, 120 cells
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0003];  % With fac=0.999 (didn't try 0.9999), stabilizes. dapd is too low.
%240, 140 cells
%                  KIred1((numpart+1):(numpart+40))= -0.02*ones(1,40); KIred1(40) = [-0.0005];  % With fac=0.999 (didn't try 0.9999), stabilizes. dapd is too low.
%240, 160 cells
%                  KIred1((numpart+1):(numpart+40))= -0.025*ones(1,40); KIred1(40) = [-0.0006];  % With fac=0.999 (didn't try 0.9999), stabilizes. dapd is too low.

% Estimator feedback
%KIred1([(numpart+sensorlocindices(1))])= -5;
%KIred1([(numpart+sensorlocindices(2))])= -5;

                  %                 % If gain is an argument
%if length(varargin) >= 1 & KIred1 ~= 0
if length(varargin) >= 1 & varargin{1} ~= 0
%if length(varargin) == 1
%KIred1 = Karg; 
KIred1 = varargin{1}; 
end
%                    fac =0.99;
                    fac =0.999;
%                    fac =0.995; % 
%                    fac =0.9995; % don't see much difference
%                    fac =0.9997; %
%                    fac =0.9993; % 
%                     if trueindex > 200000
%                         fac=0.999;
%                     end
%                    fac =0.9997; % 
%                    fac =0.9999; % 
%                    fac =0.99999; % 
                    %                fac =0.998; % this value doesn't seem to improve the bad transients with -0.01, -0.005
                    %                fac =0.997; % this value doesn't seem to improve the bad transients with -0.01
                    %                usf(k) = usf(k-1)-KIred1*sext;
                    % ***               sext2=[Verror(:,k); Verror(:,k)-fac*Verror(:,k-1)];
%                    sext2=[Verror(:,k); Verror(:,k)-fac*Verrorprev];
                    noisenow = noiseamp*randn(size(Verror(:,k))); 
%                    
                    sext2=[Verror(:,k)-noisenow; Verror(:,k)-noisenow-fac*(Verrorprev-noiseprev)];
%
                    
%                    sext2=[Verror(:,k)-noisenow; sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))]; % in combo with 165,80, KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30:32) = [-0.0002]; fac=0.999, appears to produce a wandering solution (intermediate block but not as bad as 2:1?) 
%                    sext2=[Verror(:,k)-noisenow; sign(Verror(:,k)-noisenow-fac*(Verrorprev-noiseprev))]; % in combo with 165,80, KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30:32) = [-0.0002]; fac=0.999, blows up around 8000ms
%                     sext2=[sign(Verror(:,k)-noisenow); sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))]; % in combo with 165,80, KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30:32) = [-0.0002]; fac=0.999, appears to have no effect on 2:1 block 
%  if mean(abs(Verror(:,k))) < 5%10%< 25 %> 25%> 50 %> 10
%                      sext2=[Verror(:,k)-noisenow; Verror(:,k)-noisenow-fac*(Verrorprev-noiseprev)];
%  else
%                      sext2=[Verror(:,k)-noisenow; sign(Verror(:,k)-noisenow)-sign(fac*(Verrorprev-noiseprev))]; % in combo with 165,80, KIred1((numpart+1):(numpart+30))= -0.1*ones(1,30);  KIred1(30:32) = [-0.0002]; fac=0.999, appears to produce a wandering solution (intermediate block but not as bad as 2:1?) 
%  end
                    noiseprev = noisenow; 
                    
%                    sext2=[nerror(:,k); nerror(:,k)-fac*nerrorprev]; % use this to apply n-based feedback to the V-dynamics

% ***               usf(k) = 1*(fac*usf(k-1)-KIred1*sext2);
                    usf(k) = 1*(fac*usfprev-KIred1*sext2);
                    %                 if k==round(10/deltat)
                    %                     iext(20,k)= 0.01*stimheight/deltat/cc_ext;
                    %                 end
                    %                usf(k) = -K1*[Verror(:,k); nerror(:,k)];

                    %                 if k >= tempsi + 500 && k < tempsi + 550 % % 500 to 1000 deadband is too long; constant offset is removed, but the result is not better than P control
                    %                     %k == tempsi + 500 % still allows constant offset to
                    %                     %persist in usf
                    %                     usf(k) = 0;
                    %                 end
                    %                usat = max([min([usf(k) stimheight/deltat/cc_ext]) -stimheight/deltat/cc_ext]);
                    %                usf(k) = usf(k) + 0.1*(usat-usf(k));
                end
            end
            %eigacl1(:,k) = eig(A-B1*K1);
            %         % Anti-windup scheme (turn off integrator if usf is outside saturation
            %         % limits):
            %         if usf(k) < 0
            %             usf(k) = 0;
            %         elseif usf(k) > stimheight/deltat/cc_ext
            %             usf(k) = stimheight/deltat/cc_ext;
            %         end
            % iextd is defined to be negative when current is injected
            %         if usf(k) < 0 & iextd(stimlocindex,k) == 0 % iext=iextd-usf > 0
            %             usf(k) = 0;
            %         elseif usf(k) < 0 & iextd(stimlocindex,k) == -stimheight/deltat/cc_ext
            %             % usf and iextd have the same sign.  iext=iextd-usf could be +
            %             % or - or zero. Require iextd < usf for iext < 0. Only use usf
            %             % if its magnitude is less than the pulse height.
            %             usf(k) = max([usf(k) -stimheight/deltat/cc_ext]);
            %         elseif usf(k) > 0 & iextd(stimlocindex,k) == 0
            %             % iext=iextd-usf will be negative. Only use usf
            %             % if its magnitude is less than the pulse height.
            %             usf(k) = min([usf(k) stimheight/deltat/cc_ext]); % increase threshold here
            %         elseif usf(k) > 0 & iextd(stimlocindex,k) == -stimheight/deltat/cc_ext
            %             % usf and iextd have opposite signs.  iext=iextd-usf will be negative
            %             % with magnitude greater than the pulse height.
            %             usf(k) = 0;
            %         end
            %
            % if usf(k) > iextd(stimlocindex,k) + stimheight/deltat/cc_ext
            %     usf(k) = iextd(stimlocindex,k) + stimheight/deltat/cc_ext;
            % elseif usf(k) < iextd(stimlocindex,k) - stimheight/deltat/cc_ext
            %     usf(k) = iextd(stimlocindex,k) - stimheight/deltat/cc_ext;
            % end
            if isnan(usf(k))
                usf(k) = 0;
            end

            %         if abs(usf(k)) > stimheight/deltat/cc_ext
            %             usf(k) = stimheight/deltat/cc_ext*sign(usf(k));
            %         end
            %
            iext([stimlocindex],k) = iextd(stimlocindex,k)-usf(k);
            % To move usf to cell 2, comment out the previous line and
            % uncomment the next 3, and alter KIred so it only has nonzero gains on cells 2:end
            %         iext([stimlocindexold],k) = iextd(stimlocindexold,k);
            %         stimlocindex = stimlocindexold+1;
            %         iext([stimlocindex],k) = -usf(k);

            if abs(iext(stimlocindex,k)) > stimheight/deltat/cc_ext
                iext(stimlocindex,k) = stimheight/deltat/cc_ext*sign(iext(stimlocindex,k));
            end

            % implement on LTV model
            %        usf(k) = iextd(stimlocindex,k) - iext(stimlocindex,k);
            %        usf(k) = -K1*svector(:,k);
            %        usf(k) = -Kred*svector(1:2,k);
            %        svector(:,k+1) = A*svector(:,k) + B1*usf(k);
            %       svector(:,k+1) = A*svector(:,k) + B1*usf(k);
            %        s(:,k+1) = (A-B1*KIred1)*s(:,k);
            %                    s(:,trueindex+1) = (A-B1*KIred1)*s(:,trueindex);
  end


        if pubyear ~= 'ext' & ~karma_test & sum(k == stimindices_new) %k == 2132 %799
            % stimulate periodically, at designated location
            %        stimtrue = sum(k==stimindices);
            %        V(stimlocindex,k) = stimtrue*2 + (1-stimtrue)*V(stimlocindex,k);
            V(stimlocindex,k) = stimheight;
            % ??? What's the activation threshold ??? I thought it should be
            % V=Vn=1, but it appears that somewhat lower values work, also?
        end
        %    if deton
        %        vhigh(k) = (V(detlocindex,k) > detthresh);
        %    end
        if numpart > 1
            if bc == 'n'
                %                 V(1,k+1) = V(1,k) + epsil1*deltat*(-Vx_up(k)*deltax - V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                %                 V(numpart,k+1) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k) - V(numpart,k) + Vx_dn(k)*deltax)/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
                Vnext(1) = V(1,k) + epsil1*deltat*(-Vx_up(k)*deltax - V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                Vnext(numpart) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k) - V(numpart,k) + Vx_dn(k)*deltax)/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
            elseif bc == 'd'
                %                 V(1,k+1) = V(1,k) + epsil1*deltat*(Vup(k)-2*V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                %                 V(numpart,k+1) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k)-2*V(numpart,k)+Vdn(k))/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
                Vnext(1) = V(1,k) + epsil1*deltat*(Vup(k)-2*V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                Vnext(numpart) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k)-2*V(numpart,k)+Vdn(k))/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
            elseif bc == 'p'
                %                 V(1,k+1) = V(1,k) + epsil1*deltat*(Vup(k)-2*V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                %                 V(numpart,k+1) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k)-2*V(numpart,k)+Vdn(k))/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                %                     - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
                Vnext(1) = V(1,k) + epsil1*deltat*(Vup(k)-2*V(1,k)+V(2,k))/(deltax)^2 + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
                Vnext(numpart) = V(numpart,k) + epsil1*deltat*(V(numpart-1,k)-2*V(numpart,k)+Vdn(k))/(deltax)^2 + deltat*karma_f(V(numpart,k),n(numpart,k),nB,M)/epsil2 ...
                    - deltat*cc_ext*iext(numpart,k) + deltat*cc_int*iint(numpart,k);
                if k>1
                    Vdnnext = Vnext(numpart);
                    Vupnext = Vdnnext;
                    %                     Vdn(k+1) = V(numpart,k+1);
                    %                     Vup(k+1) = Vdn(k+1);
                    %                     %              Vup(k+1) = V(1,k+1);
                    %                     %              Vdn(k+1) = Vup(k+1);
                end
            end
            Vnext(2:(numpart-1)) = V(2:(numpart-1),k) + epsil1*deltat*(V(1:(numpart-2),k)-2*V(2:(numpart-1),k)+V(3:numpart,k))/(deltax)^2 + deltat*karma_f(V(2:(numpart-1),k),n(2:(numpart-1),k),nB,M)/epsil2 ...
                - deltat*cc_ext*iext(2:(numpart-1),k) + deltat*cc_int*iint(2:(numpart-1),k);
            %             V(2:(numpart-1),k+1) = V(2:(numpart-1),k) + epsil1*deltat*(V(1:(numpart-2),k)-2*V(2:(numpart-1),k)+V(3:numpart,k))/(deltax)^2 + deltat*karma_f(V(2:(numpart-1),k),n(2:(numpart-1),k),nB,M)/epsil2 ...
            %                 - deltat*cc_ext*iext(2:(numpart-1),k) + deltat*cc_int*iint(2:(numpart-1),k);
        else
            Vnext(1) = V(1,k) + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
                - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
            %             V(1,k+1) = V(1,k) + deltat*karma_f(V(1,k),n(1,k),nB,M)/epsil2 ...
            %                 - deltat*cc_ext*iext(1,k) + deltat*cc_int*iint(1,k);
        end

        %         if sum(k == stimindices_new)
        %             V(stimlocindex,k+1) = - deltat*cc_ext*iext(stimlocindex,k) + deltat*cc_int*iint(stimlocindex,k);
        %         end
        %    V(:,k+1) = max([V(:,k+1) -4*ones(size(V(:,k+1)))],[],2); % Must remove
        %    if using Rappel units
        %        V(:,k+1) = max([V(:,k+1) 0*ones(size(V(:,k+1)))],[],2);
        %    n(:,k+1) = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN;
%%            if k > fbon_time
%            if trueindex > fbon_time
%                un(:,k) = -0.005*deltat*Verror(:,k)/tauN; % proportional feedback term applied to n instead of V
%                gains in foregoing line (using Verror for fb): positive
%                gains, such as 0.00075*deltat/tauN, don't appear to work.
%                However, -0.005*deltat/tauN stabilizes, -0.003 stabilizes slowly, and -0.001 doesn't. 
               
%                un(:,k) = 0.75*deltat*nerror(:,k)/tauN; % proportional feedback term applied to n instead of V
%                gains in foregoing line: 1*deltat/tauN and 0.75*deltat/tauN stabilize, 0.5*deltat/tauN stabilizes very
%                slowly, 0.25*deltat/tauN does not stabilize
%%                un(:,k) = nerror(:,k)/tauN; % proportional feedback term applied to n instead of V
%            end


            if trueindex > fbon_time % This set should be consistent with b2 = -0.05
%                un(:,k) = 0.005*Verror(:,k)/tauN; % proportional feedback term applied to n instead of V
%                gains in foregoing line (using Verror for fb): negative
%                gains, such as -0.005/tauN, don't appear to work.
%                However, 0.005*deltat/tauN stabilizes, 0.003 stabilizes slowly, and 0.001 doesn't. 
%                un(:,k) = 0.05*Verror(:,k)/tauN; % very aggressive; no apparent oscillations
%                un(:,k) = 0.01*Verror(:,k)/tauN; % several apparent oscillations
%                un(:,k) = 0.03*Verror(:,k)/tauN; % almost 1 oscillation (debatable)
%                un(sensorlocindices(1),k) = 0.03*Verror(sensorlocindices(1),k)/tauN; % 
%                un(sensorlocindices,k) = 0.03*Verror(sensorlocindices(1),k)/tauN; % 
%                un(:,k) = 0.03*Verror(sensorlocindices(1),k)/tauN; % destabilizing? 
                un(:,k) = 0.01*Verror(sensorlocindices(1),k)/tauN; % destabilizing? 
                
 %               un(:,k) = 0.04*Verror(:,k)/tauN; % almost 1 oscillation (debatable)
               
%                un(:,k) = -0.75*nerror(:,k)/tauN; % proportional feedback term applied to n instead of V
%                gains in foregoing line: -1/tauN and -0.75/tauN stabilize, -0.5/tauN stabilizes very
%                slowly, -0.25/tauN does not stabilize
            end
%        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k); % This one is consistent with b2=1 assumption
        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN - deltat*un(:,k); % This should be consistent with b2 = -0.05


%        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k); % This one is consistent with b2=1 assumption
%%        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + deltat*un(:,k);
        %        n(:,k+1) = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k);
        % Apply a small perturbation to V in the first cell near the time when
        % it crosses V=1
        % Approximate times for V(1,:)
        % time    Vcross
        % 3.4832  0.928 up
        % 3.4834  0.998 up
        % 3.4836  1.07  up
        % 4.174   1.01  dn
        % 4.1742  1.00  dn
        % 4.1744  0.989 dn
        % Approximate times for V(5,:)
        % time    Vcross
        % 4.2638   1.01  dn
        % Approximate times for V(5,:), bcl 0.92
        % time    Vcross
        % 4.511   1.01  dn
        % Approximate times for V(5,:), bcl 0.97
        % time    Vcross
        % 4.7498   1.01  dn
        % Approximate times for V(5,:), bcl 1.12
        % time    Vcross
        % 5.441   1.02  dn
        if trueindex == 1 % ***k==1
            tpert = 5.441;
            xpert = 5*deltax;
        end
        if 0&trueindex==round(tpert/deltat)%***k==round(tpert/deltat)%k==round(4.173/deltat)%%k==round(3.5672/deltat)%k==round(3.484/deltat)%k>=round(4.173/deltat) & k<round(4.173/deltat)+20%
            %        V(1,k+1) = V(1,k+1) + 1e-6;
            %        disp('.')
            %        V(1,k+1) = V(1,k+1) + 1e-5;
            %        V(1,k+1) = V(1,k+1) + 1e-3;
            %        V(round(xpert/deltax),k+1) = V(round(xpert/deltax),k+1) + 1e-2;
            Vnext(round(xpert/deltax)) = Vnext(round(xpert/deltax)) + 1e-2;
            %        V(1,k+1) = V(1,k+1) + 1e-2;
            %        V(1,k+1) = V(1,k+1) + 5e-1;
            %        V(1,k+1) = V(1,k+1) + 1;
            %        V(1,k+1) = V(1,k+1) + 10;
        end

        %    Verror(:,k) = Vd(:,k)-V(:,k);
%         %     if k == round(670/deltat)
%         %         V(1,k+1) = V(1,k+1) + 1e-1;
%         %         n(1,k+1) = n(1,k) + 1e-10;
%         %     end
%             if trueindex == round(670/deltat)
%                 Vnext = Vnext + 1e-1;
%                 nnext = nnext + 1e-2;
%             end
if abs(Vnext) > 500
    Vnext
    break
end
        Verrornext = Vd(:,k+1)-Vnext;
        nerrornext = nd(:,k+1)-nnext;
        %     Verror(:,k+1) = Vd(:,k+1)-V(:,k+1);
        %     nerror(:,k+1) = nd(:,k+1)-n(:,k+1);
%         %     if k==round(670/deltat) % periodically reset the LTV system to the correct IC
%         %         s(:,k+1) = [Verror(1,k+1);nerror(1,k+1)];
%         %     end
%             if trueindex==round(670/deltat) % periodically reset the LTV system to the correct IC
%                 s(:,trueindex+1) = [Verrornext;nerrornext];
%             end

        if deton & (trueindex > 1) %***(k > 1)
            %***        dtest = (V(detlocindices,k) > detthresh) & (V(detlocindices,k-1) <= detthresh);  % rising edge detected
            dtest = (V(detlocindices,k) > detthresh) & (Vprev(detlocindices) <= detthresh);  % rising edge detected
            dwhichdets = find(dtest.*detvec); % indices of detectors that show a rising edge at this timestep
            dctr = dctr + dtest.*ones(numdet,1); % increment counter for activated detectors
            if sum(dtest)
                for jj=dwhichdets'
                    apstartindices{jj}(dctr(jj)) = trueindex; %*** k;
                    apedges{jj} = [apedges{jj} apstartindices{jj}(dctr(jj))];
                end
            end
            % ***        rtest = (V(detlocindices,k) <= detthresh) & (V(detlocindices,k-1) > detthresh); % falling edge detected
            rtest = (V(detlocindices,k) <= detthresh) & (Vprev(detlocindices) > detthresh); % falling edge detected
            rwhichdets = find(rtest.*detvec); % indices of detectors that show a falling edge at this timestep
            rctr = rctr + rtest.*ones(numdet,1); % increment counter for activated detectors
            if sum(rtest)
                for kk=rwhichdets'
                    %***                apendindices{kk}(rctr(kk)) = k-1;
                    apendindices{kk}(rctr(kk)) = trueindex-1;
                    apedges{kk} = [apedges{kk} apendindices{kk}(rctr(kk))];
                    %            end
                    %            for ii=detvec
                    edgediff{kk}=diff(apedges{kk});
                    % *** Warning: the following assumes that the smallest element of
                    % apedges is a rising (depolarization) edge, rather than a falling
                    % edge
                    apds{kk}=edgediff{kk}(1:2:end); % AP durations (units = number of timesteps)
                end
                %        end
                %         if V(detlocindices,k) > detthresh & V(detlocindices,k-1) <= detthresh
                % %if V(detlocindex,k+1) > detthresh & V(detlocindex,k) <= detthresh
                %             % rising edge detected
                %             dctr = dctr+1;
                %             apstartindex(dctr) = k;%k+1;
                %             apedges = [apedges apstartindex(dctr)];
                %         end
                %         if V(detlocindices,k) <= detthresh & V(detlocindices,k-1) > detthresh
                % %if V(detlocindex,k+1) <= detthresh & V(detlocindex,k) > detthresh
                %             % falling edge detected
                %             rctr = rctr+1;
                %             apendindex(rctr) = k-1;%k;
                %             apedges = [apedges apendindex(rctr)];
                %             edgediff=diff(apedges);
                %         % *** Warning: the following assumes that the smallest element of
                %         % apedges is a rising (depolarization) edge, rather than a falling
                %         % edge
                %             apds=edgediff(1:2:end); % AP durations (units = number of timesteps)
                %            if strcmp(fbon,'pyra') % FOR MULTIPLE DETECTORS, THIS IS WRONG
                %            (only want to perturb timing when feedback detector detects a
                %            falling edge, not when any detector detects a
                %            falling edge)
            end % if falling edge
        end % if deton
        if ~mod(trueindex,1000) % *** ~mod(k,1000)
            disp(trueindex) % *** disp(k)
        end
        Vprev = V(:,k);
        nprev = n(:,k);
        Verrorprev = Verror(:,k);
        nerrorprev = nerror(:,k);
        usfprev = usf(k);
        if k < writeintsteps
            % write future values to array
            V(:,k+1) = Vnext;
            n(:,k+1) = nnext;
            Verror(:,k+1) = Verrornext;
            nerror(:,k+1) = nerrornext;
            Vdn(k+1) = Vdnnext;
            Vup(k+1) = Vupnext;
        end
    end % for k time loop
%     %eval(['save data/' num2str(ii) ' hprime acoeff* usf Verror nerror iext V n un']);  % s iint Vx_dn Vx_up Vdn Vup
%     eval(['save data/' num2str(ii) ' usf Verror nerror iext V n']);  % s iint Vx_dn Vx_up Vdn Vup
%     clear hprime acoeff* usf Verror nerror iext V n un
%    eval(['save data/' num2str(ii) ' usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
%    eval(['save data/' num2str(ii) ' usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
%    eval(['save data/' num2str(ii) ' un hprime acoeff* usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
    eval(['save data/' num2str(ii) ' un usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
    clear hprime acoeff* usf Verror nerror iext V n un 
    if ii == numrep
%        save data/configinfo stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec dctr rctr fbon_time KIred1 stimindices iextd    
%        save data/configinfo B1 B2 stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength
%        save data/configinfo B1 B2 stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength
%        save data/configinfo stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength
%        save data/configinfo stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength stimlocindex dapd
        save data/configinfo stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength stimlocindex dapd desfile
    end
end % data write loop
toc
k
% Store time-averages of pole-placement K gains
%Kstore = mean(kgains,2)' % this should be modified to skip ahead past fbon_time
% Originally, including all the zeros in the avg (for pp reduction of
% 0.999) Kstore =  -0.2834   -0.0998
% (not much effect) Excluding the zeros from the avg gives   -0.5015   -0.1766 (no good effect, appears to be destabilizing)
% For pp reduction of 0.995, the old values were   -1.4167   -1.5043.
%Kstore = mean(kgains(:,(fbon_time+1):end),2)'
%save kgains_2cell_b087_Ared_pp0_999 Kstore
%save kgains_2cell_b087_Ared_pp0_995 Kstore
%save kgains_2cell_b087_Ared_pp0_99 Kstore
if 0 % Controllability Grammians
    numdiv = 20; % number of subdivisions of a single AP
    cgrams= cell(1,numdiv);
    cgramns=cell(1,numdiv);
    svdcgs=cell(1,numdiv);
    svdcgns=cell(1,numdiv);
    
    %    init=round(3.483/deltat); % start time of a selected AP
    %    increment = 0.2; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    %    init=round(450.05/deltat); % start time of a selected AP
    init=round(650.05/deltat); % start time of a selected AP
    %    increment = stimperiod(1)/5; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    increment = stimperiod(1)/numdiv; % number of time units between starttimes (divide a bcl into about 5 equal intervals)
    for nn=1:numdiv
        cgram = zeros(2,2);
        cgramn = zeros(2,2);
        Bn = [0;1]; % B matrix for n as input?
        %endstep = numstep-1;
        %       startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        %       endstep = stimperiodsteps(1)+startstep;%round(stimperiodsteps(1)/2)+startstep;%1*stimperiodsteps(1);
        startstep = init+(nn-1)*round(increment/deltat);%round(3.483/deltat);%round(6.9632/deltat);%1; % Put near V=1 crossings? ("controllable from initial time")
        endstep = startstep+5;%+2;%Rappel-units model has very poor scaling, so shouldn't run for very many timesteps to prevent blowup of Grammians
        svdcg=zeros(2,endstep-startstep+1);
        svdcgn=zeros(2,endstep-startstep+1);
        for kk=startstep:endstep
            % original version was wrong
            % Define state transition matrix
            stmktoend = eye(2);
            for jj = kk:endstep-1 %kk:(endstep) % wrong
                %            stmktoend = reshape(allAr(:,jj),2,2)*stmktoend; %A(numstep-1)...A(kk)
                Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
                A12 = [deltat*acoeff2(1,jj)];
                Auc = [1-deltat/tauN];
                A21 = [deltat*acoeff3(1,jj)];
                A = [ Ac   A12
                    A21  Auc]
                stmktoend = A*stmktoend; %A(numstep-1)...A(kk)
            end
            cgram = cgram + stmktoend*B1*B1'*stmktoend';
            cgramn = cgramn + stmktoend*Bn*Bn'*stmktoend';
            svdcg(:,kk) = svd(cgram);
            svdcgn(:,kk) = svd(cgramn);
            if ~mod(kk,1000)
                disp(kk)
            end
            %         if kk==round(endstep/2) % save middle value of grammians
            %             cgrammid = cgram;
            %             cgramnmid = cgramn;
            %         end
        end
        cgrams{nn} = cgram;
        cgramns{nn} = cgramn;
        svdcgs{nn} = svdcg;
        svdcgns{nn} = svdcgn;
    end
    % should probably let this run for several periods. So far, grammian results are
    % consistent with the instantaneous results, but should (1) run for a
    % longer number of bcls (improving the efficiency of the routine would
    % help), and (2) verify the LTV model. In addition, should compute a few
    % steps by hand to confirm that the calculations are sensible.
    %save ctrbgram_poster startstep endstep cgram* svd* stimperiod numpart numstep deltax deltat increment init

    figure
    subplot(2,1,1)
    hold on;
    subplot(2,1,2)
    hold on;
    % for nn=1:5
    %     startstep = init+(nn-1)*round(increment/deltat);
    % %    plot(startstep*deltat,svd(cgrams{nn}),'*')
    % %    plot(startstep*deltat,svd(cgramns{nn}),'*')
    %     plot(startstep*deltat,min(svd(cgrams{nn})),'b*')
    %     plot(startstep*deltat,min(svd(cgramns{nn})),'gs')
    % end
    for nn=1:numdiv
        startstep = init+(nn-1)*round(increment/deltat);
        %    plot(startstep*deltat,svd(cgrams{nn}),'*')
        %    plot(startstep*deltat,svd(cgramns{nn}),'*')
        subplot(2,1,1)
        p1(nn)=plot(startstep*deltat,min(svd(cgrams{nn})),'b*');
        subplot(2,1,2)
        p2(nn)=plot(startstep*deltat,min(svd(cgramns{nn})),'gs');
    end
    subplot(2,1,1)
    yl=ylabel('min singular value');
    p3=plot((init:endstep)*deltat,4e-23*(V(:,init:endstep)+85)/85,'--');
    legend('control V directly')
    ti=title('Minimum singular values of controllability Grammians');
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    subplot(2,1,2)
    yl=ylabel('min singular value');
    legend('control n directly','Location','southeast')
    xl=xlabel('start time for Grammian');
    p3=plot((init:endstep)*deltat,5*(V(:,init:endstep)+85)/85,'--');
    if poster
        set(gca,'FontSize',fs);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
    end

    figure
    subplot(2,1,1)
    hold on;
    subplot(2,1,2)
    hold on;
    % for nn=1:5
    %     startstep = init+(nn-1)*round(increment/deltat);
    % %    plot(startstep*deltat,svd(cgrams{nn}),'*')
    % %    plot(startstep*deltat,svd(cgramns{nn}),'*')
    %     plot(startstep*deltat,min(svd(cgrams{nn})),'b*')
    %     plot(startstep*deltat,min(svd(cgramns{nn})),'gs')
    % end
    for nn=1:numdiv
        startstep = init+(nn-1)*round(increment/deltat);
        %    plot(startstep*deltat,svd(cgrams{nn}),'*')
        %    plot(startstep*deltat,svd(cgramns{nn}),'*')
        subplot(2,1,1)
        p1(nn)=plot(startstep*deltat,max(svd(cgrams{nn})),'b*');
        subplot(2,1,2)
        p2(nn)=plot(startstep*deltat,max(svd(cgramns{nn})),'gs');
    end
    subplot(2,1,1)
    yl=ylabel('max singular value');
    p3=plot((init:endstep)*deltat,4e-23*(V(:,init:endstep)+85)/85,'--');
    legend('control V directly')
    ti=title('maximum singular values of controllability Grammians');
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    subplot(2,1,2)
    yl=ylabel('max singular value');
    legend('control n directly','Location','southeast')
    xl=xlabel('start time for Grammian');
    p3=plot((init:endstep)*deltat,5*(V(:,init:endstep)+85)/85,'--');
    if poster
        set(gca,'FontSize',fs);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
    end



end %if 1

% Try to check whether scaling of Bn matters. Replacing Bn with [0;B1(1)]
% % will change the grammian by a factor of B1(1)^2
% figure
% hold on;
% for nn=1:5
%     startstep = init+(nn-1)*round(increment/deltat);
% %    plot(startstep*deltat,svd(cgrams{nn}),'*')
% %    plot(startstep*deltat,svd(cgramns{nn}),'*')
%     plot(startstep*deltat,min(svd(cgrams{nn})),'b*')
%     plot(startstep*deltat,min(svd(B1(1)^2*cgramns{nn})),'gs')
% end

% figure
% plot((1:endstep)*deltat,svdcg)
% figure
% plot((1:endstep)*deltat,svdcg(1,:)./svdcg(2,:),'r:')
% figure
% plot((1:endstep)*deltat,svdcgn)
% figure
% plot((1:endstep)*deltat,svdcgn(1,:)./svdcgn(2,:),'r:')

% check a selected value manually
%A0=[0.9778 0; 0 0.9998]
%A1=[0.9821 0; 0 0.9998]
%A2=[0.9864 0; 0.0001 0.9998]
%cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e-004 *
%    0.395979625244837                   0
%                    0                   0
% jj=1
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A0 = [ Ac   A12
%             A21  Auc];
% jj=2
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A1 = [ Ac   A12
%             A21  Auc];
%
% jj=3
% Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,jj)+1];
%         A12 = [deltat*acoeff2(1,jj)];
%         A21 = [deltat*acoeff3(1,jj)];
%         A2 = [ Ac   A12
%             A21  Auc];
% cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e-004 *
%    0.395999903558393   0.000017710454326
%    0.000017710454326   0.000000001335500
% cgram with endstep = 3 (new)
%   1.0e-004 *
%    0.395999903558393   0.000017710454326
%    0.000017710454326   0.000000001335500
% cgram with endstep=3 (old):
%   1.0e-004 *
%    0.385323995616494   0.000037221998733
%    0.000037221998733   0.000000004138828
%
% Redo manual test for Rappel units
%A0 = 1e4*[2.251320187928600  -0.000001496755894; 0.000000000000000 0.000099980000000]
%A1 = 1e4*[2.166800616612064  -0.000001438301991; 0.000000000000000   0.000099980000000]
%A2 = 1e4*[2.084614165392807  -0.000001381583702; 0.000000000000000   0.000099980000000]
% cgramcheck1to3=A1*A0*B1*B1'*A0'*A1' + A1*B1*B1'*A1' + B1*B1'
%   1.0e+014 *
%    5.949116066528354                   0
%                    0                   0
% cgram =
%   1.0e+014 *
%    5.949116066528354   0.000000000000000
%    0.000000000000000   0.000000000000000

%if 1%length(varargin) == 1
    % h1=figure;
    % hold on;
    % %plot(1:numpart,V(:,1:1000:end));
    % plot((1:numpart)*deltax,V(:,k+1),'g-o');
    % %xlabel('spatial grid element')
    % xlabel('distance')
    % ylabel('V')
    %
    % h2=figure;
    % hold on;
    % %plot(1:numpart,n(:,1:1000:end));
    % plot((1:numpart)*deltax,n(:,k+1),'g-o');
    % %xlabel('spatial grid element')
    % xlabel('distance')
    % ylabel('n')

    % figure
    % hold
    % p1=plot((1:numstep)*deltat, abs(eiga'));
    % p2=plot((1:numstep)*deltat, abs(eigacl1'),'.');
    % plot((1:numstep)*deltat, abs(0.995*eiga'),'m--');
    % legend([p1(1); p2(1)],'A','A-BK')
    % title('A and A-BK eigenvalue moduli vs. time')

    % figure
    % hold
    % p1=plot((1:numstep)*deltat, abs(eigared'));
    % p2=plot((1:numstep)*deltat, abs(eigaclred'),'.');
    % plot((1:numstep)*deltat, abs(0.99*eigared'),'m--');
    % legend([p1(1); p2(1)],'A','A-BK')
    % %title('A and A-BK eigenvalue moduli vs. time')

    % figure
    % hold
    % p1=plot((1:numstep)*deltat, abs(eigata'));
    % p2=plot((1:numstep)*deltat, abs(eigatacl'),'.');
    % %plot((1:numstep)*deltat, abs(0.995*eiga'),'m--');
    % legend([p1(1); p2(1)],'A','A-BK')
    % title('A''A and (A-BK)''(A-BK) eigenvalue moduli vs. time')
    %
    % figure
    % hold
    % p1=plot((1:numstep)*deltat, abs(eigatacl')-abs(eigata'));
    % %plot((1:numstep)*deltat, abs(0.995*eiga'),'m--');
    % title('(A-BK)''(A-BK)- A''A  eigenvalue moduli vs. time')

    % figure
    % plot((1:numstep)*deltat,kgains);

    %figure
    %plot((1:numstep)*deltat,V','.-')
    %
    % figure
    % plot((1:numstep)*deltat,n','.-')
    %
    % figure
    % plot((1:numstep)*deltat,iext','.-')
    %     % compare linearized and nonlinear models
    %     figure
    %     hold on;
    %     %     plot(dt*(1:2:numstep),V(:,1:2:numstep));
    %     %     plot(dt*(1:2:numstep),n(:,1:2:numstep));
    %     p1=plot(deltat*(1:numstep),Verror(:,1:numstep));
    %     p2=plot(deltat*(1:numstep),nerror(:,1:numstep),'r');
    %     p3=plot(deltat*(1:numstep),s(:,1:numstep),'g--');
    %     legend([p1(1);p2(1);p3(1)],'Verror (nonl)', 'nerror(nonl)', 'LTV states 1 and 2')
    %     title('Linearized and nonlinear states')
% 
% % Enable below to plot n-perturbations instead of V-perturbations    
% if 0 
%     h=figure;
%     clear p p1;
%     for ii=1:numrep % number of data writing cycles
%         subplot(2,1,1)
%         hold on;
%         eval(['load data/' num2str(ii)])
%        p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps),'b');
%         clear hprime acoeff* usf Verror nerror iext V n un
%     end
%     yl=ylabel('perturbation to n_{ t}');
%     ti=title(['State feedback term vs. time']);
%     %    axis([0 finaltime -20 20])
%     axis([0 finaltime -6e-5 6e-5])
%     if poster
%         set(gca,'FontSize',fs);
%         set(p1,'LineWidth',lw1);
%         set(yl,'FontSize',fs);
%         set(ti,'FontSize',fs);
%     end
% end
h=figure;
clear p p1 p2;
for ii=1:numrep % number of data writing cycles
    subplot(2,1,1)
    hold on;
    eval(['load data/' num2str(ii)])
    p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([stimlocindex],1:writeintsteps),'b');
    p2(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iextd([stimlocindex],1:writeintsteps),'r--');
    %         p1(ii)=plot(([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat,iext([stimlocindex],1:10:writeintsteps),'b');
    %         p2(ii)=plot(([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat,iextd([stimlocindex],1:10:writeintsteps),'r--');
    %    p1(ii)=plot((1:numstep)*deltat,iext([stimlocindex],1:numstep),'b');
    %    p2(ii)=plot((1:numstep)*deltat,iextd([stimlocindex],1:numstep),'r--');
    clear hprime acoeff* usf Verror nerror iext V n un
end
legend('i','i^d','Location','West')
%    xl=xlabel('time');
yl=ylabel('current (\muA/cm^2)');
ti=title(['Control current vs. time']);
%    axis([0 finaltime -20 20])
axis([0 finaltime -35 35])
if poster
    set(gca,'FontSize',fs);
    set(p1,'LineWidth',lw1);
    set(p2,'LineWidth',lw2);
    %       set(xl,'FontSize',fs);
    set(yl,'FontSize',fs);
    set(ti,'FontSize',fs);
end

% Not sure how best to obtain the color order, so open and close a
% figure:
if ~exist('psyms')
    hfig = figure
    colororder=get(gca,'ColorOrder');
    close(hfig);
    nrco = size(colororder,1); % number of rows (distinct colors)

    %    psyms = strvcat('g-o','--*','m-s','k-x','r-+','y-h','c-^','g-v','--x','m-+','k-h');
    %    psyms = strvcat('+','o','*','x','s','d','^','v','>','<','p','h');    %same order as in Matlab help file
    psyms = repmat(strvcat('o','*','s','+','x','d','^','v','>','<','p','h'),2,1);
end

%    detvec = (71:80)';
colorcounter = 1;
markerindex = 1;
subplot(2,1,2)
for i=detvec'
    hold on;
    if ~isempty(apds{i})
        if length(apendindices{i}) ~= length(apds{i})
            p(i) = plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
            % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
            % stimulus is below peak AP height, (c) BCL is less than initial APD
        else
            p(i) = plot(deltat*apendindices{i}, deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
            %              axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
            axis([0 finaltime 0 200])%axis([0 10 .12 .14])%
            %                axis([0 finaltime 0 275])%axis([0 10 .12 .14])%
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
end

if 1
    % Define an offset colororder
    colororderoffset = colororder([3:7 1:2],:);

    %    numrep = 20;
    h=figure;
    crange = sensorlocindices %1:numpart; %detlocindices
    %    crange = 1; %detlocindices
    clear p1 p2;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        %    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
        set(gca,'ColorOrder',colororderoffset);
        %    p2=plot(deltat*(1:1:numstep),Vd(detlocindices,1:1:numstep),'--');
        %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Vd(detlocindices,1:writeintsteps),'--');
%        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(crange,1:50:writeintsteps),'--');
        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vinterp([((ii-1)*writeintsteps+1):50:ii*writeintsteps],sensorindices),'--');
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(detlocindices,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange,1:50:writeintsteps));
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(crange,1:50:writeintsteps)-V(crange,1:50:writeintsteps));
%        plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(1,1:50:writeintsteps),'g--');
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    %    yl=ylabel('membrane potential (mV)');
    yl=ylabel({'membrane', 'potential (mV)'});
    ti=title(['V and V^d vs. time']);
    legend([p1(end); p2(end)],'V_i','V^d_i')
    %    legend([p1(end); p2(end)],'V','V^d')
    axis([0 finaltime -100 40])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        %       set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,n(detlocindices,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,n(crange,1:50:writeintsteps));
        set(gca,'ColorOrder',colororderoffset);
        %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,nd(detlocindices,1:writeintsteps),'--');
        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,nd(crange,1:50:writeintsteps),'--');
        %    p1=plot(deltat*(1:1:numstep),n(detlocindices,1:1:numstep));
        %    set(gca,'ColorOrder',colororderoffset);
        %    p2=plot(deltat*(1:1:numstep),nd(detlocindices,1:1:numstep),'--');
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    xl=xlabel('time (ms)');
    %    yl=ylabel('ion channel variable');
    yl=ylabel({'ion channel', 'variable'});
    ti=title(['n and n^d vs. time']);
    legend([p1(end); p2(end)],'n_i','n^d_i')
    %    legend([p1(end); p2(end)],'n','n^d')
    axis([0 finaltime 0 1.2])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    h=figure;
    crange = 1:numpart; %detlocindices
    %    crange = 1; %detlocindices
    clear p1 p2;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        %    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
        set(gca,'ColorOrder',colororderoffset);
        %    p2=plot(deltat*(1:1:numstep),Vd(detlocindices,1:1:numstep),'--');
        %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Vd(detlocindices,1:writeintsteps),'--');
%        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(crange,1:50:writeintsteps),'--');
        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vinterp([((ii-1)*writeintsteps+1):50:ii*writeintsteps],sensorindices),'--');
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(detlocindices,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange,1:50:writeintsteps));
        plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(sensorlocindices,1:50:writeintsteps),'g-');
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(crange,1:50:writeintsteps)-V(crange,1:50:writeintsteps));
%        plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(1,1:50:writeintsteps),'g--');
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    %    yl=ylabel('membrane potential (mV)');
    yl=ylabel({'membrane', 'potential (mV)'});
    ti=title(['V and V^d vs. time']);
    legend([p1(end); p2(end)],'V_i','V^d_i')
    %    legend([p1(end); p2(end)],'V','V^d')
    axis([0 finaltime -100 40])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        %       set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    clear p2;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,n(detlocindices,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,n(crange,1:50:writeintsteps));
        set(gca,'ColorOrder',colororderoffset);
        %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,nd(detlocindices,1:writeintsteps),'--');
        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,nd(crange,1:50:writeintsteps),'--');
        plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,n(sensorlocindices,1:50:writeintsteps),'g-');
        %    p1=plot(deltat*(1:1:numstep),n(detlocindices,1:1:numstep));
        %    set(gca,'ColorOrder',colororderoffset);
        %    p2=plot(deltat*(1:1:numstep),nd(detlocindices,1:1:numstep),'--');
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    xl=xlabel('time (ms)');
    %    yl=ylabel('ion channel variable');
    yl=ylabel({'ion channel', 'variable'});
    ti=title(['n and n^d vs. time']);
    legend([p1(end); p2(end)],'n_i','n^d_i')
    %    legend([p1(end); p2(end)],'n','n^d')
    axis([0 finaltime 0 1.2])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    h3=figure;
    hold on;
    %            eval(['load data_80cell_b225_15000_ol/configinfo'])
    %            writeintsteps = round(writeint/deltat);
    %            fiberlength = deltax*numpart;
    reprange = (numrep-10):(numrep-1)%1:10%;%10:20;%11:20%41:50%
    for ii=reprange % number of data writing cycles
        %            eval(['load data_80cell_b225_15000_ol/' num2str(ii)])
        eval(['load data/' num2str(ii)])
        %    eval(['load data_80cell_b225_7500_ol/' num2str(ii)])
        tr=([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat;
        imagesc(deltax*(1:numpart), tr, V(:,1:10:writeintsteps)')
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    axis xy
    axis([deltax fiberlength ((reprange(1)-1)*writeintsteps+1)*deltat tr(end)])
    colorbar
    yl=ylabel('time (ms)');
    xl=xlabel('distance (cm)');
    %    ti=title('Contour plot of membrane potential, V(x,t) (mV)');
    ti=title('Membrane potential, V(x,t) (mV)');
    if poster
        set(gca,'FontSize',fs);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    h31=figure;
    hold on;
    %             eval(['load data_80cell_b225_15000_ol/configinfo'])
    %             writeintsteps = round(writeint/deltat);
    %             fiberlength = deltax*numpart;
    for ii=reprange % number of data writing cycles
        %            eval(['load data_80cell_b225_15000_ol/' num2str(ii)])
        eval(['load data/' num2str(ii)])
        tr=([((ii-1)*writeintsteps+1):10:ii*writeintsteps])*deltat;
        imagesc(deltax*(1:numpart), tr, n(:,1:10:writeintsteps)')
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    axis xy
    axis([deltax fiberlength ((reprange(1)-1)*writeintsteps+1)*deltat tr(end)])
    colorbar
    yl=ylabel('time (ms)');
    xl=xlabel('distance (cm)');
    ti=title('Contour plot of ion-channel variable, n(x,t)');
    if poster
        set(gca,'FontSize',fs);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    h8=figure;
    hold on;
    for ii=1:numrep % number of data writing cycles
        eval(['load data/' num2str(ii)])
        %    plot((1:numstep)*deltat,un([stimlocindex],1:numstep))
        plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps))
        xlabel('time')
        ylabel('npert')
        title(['Perturbation applied to n in cell ' num2str(stimlocindex)])
        clear hprime acoeff* usf Verror nerror iext V n un
    end

    h8=figure;
    hold on;
    %plot((1:numstep)*deltat,iext([stimlocindex stimoutlocindex],1:numstep))
    % plot((1:numstep)*deltat,iextd([stimlocindex stimoutlocindex],1:numstep),'--')
    % plot((1:numstep)*deltat,iext(:,1:numstep))
    %    plot((1:numstep)*deltat,iext([stimlocindex],1:numstep),'c')
    %    plot((1:numstep)*deltat,iextd([stimlocindex],1:numstep),'--')
    p1=plot((1:numstep)*deltat,iext([stimlocindex],1:numstep),'b');
    p2=plot((1:numstep)*deltat,iextd([stimlocindex],1:numstep),'r--');
    %    plot((1:numstep)*deltat,-fb,'g:')
    %    legend('iext_d','iext','-fb')
    %    legend('iext','iext_d')
    legend('i','i^d')
    xl=xlabel('time');
    %   ylabel('iext')
    yl=ylabel('current');
    %    title(['External current vs. time. Current entry location: ' num2str(stimloc) ', exit location: ' num2str(stimoutloc) ', stimulus height ' num2str(stimheight) ', stimulus duration ' num2str(stimduration)])
    %    title(['External current vs. time. Current entry location: ' num2str(stimloc) ', stimulus height ' num2str(stimheight) ', stimulus duration ' num2str(stimduration)])
    ti=title(['Control current vs. time']);
    axis([0 finaltime -20 20])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end


    h7=figure;
    hold on;
    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
    set(gca,'ColorOrder',colororderoffset);
    p2=plot(deltat*(1:1:numstep),Vd(detlocindices,1:1:numstep),'--');
    yl=ylabel('membrane potential');
    xl=xlabel('time');
    ti=title(['V and V^d vs. time']);
    legend([p1; p2],strvcat([repmat('V_',numpart,1) num2str((1:numpart)')],[repmat('V^d_',numpart,1) num2str((1:numpart)')]))
    %    legend([p1(1); p2(1)],'V','V^d')
    axis([0 finaltime -2 5])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    h6=figure;
    hold on;
    p1=plot(deltat*(1:1:numstep),n(detlocindices,1:1:numstep));
    set(gca,'ColorOrder',colororderoffset);
    p2=plot(deltat*(1:1:numstep),nd(detlocindices,1:1:numstep),'--');
    xl=xlabel('time');
    yl=ylabel('ion channel variable');
    ti=title(['n and n^d vs. time']);
    %    legend([p1(1); p2(1)],'n','n^d')
    legend([p1; p2],strvcat([repmat('n_',numpart,1) num2str((1:numpart)')],[repmat('n^d_',numpart,1) num2str((1:numpart)')]))
    axis([0 finaltime 0 1])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(p2,'LineWidth',lw2);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end


    % Not sure how best to obtain the color order, so open and close a
    % figure:
    if ~exist('psyms')
        hfig = figure
        colororder=get(gca,'ColorOrder');
        close(hfig);
        nrco = size(colororder,1); % number of rows (distinct colors)

        %    psyms = strvcat('g-o','--*','m-s','k-x','r-+','y-h','c-^','g-v','--x','m-+','k-h');
        %    psyms = strvcat('+','o','*','x','s','d','^','v','>','<','p','h');    %same order as in Matlab help file
        psyms = strvcat('o','*','s','+','x','d','^','v','>','<','p','h');
    end

    colorcounter = 1;
    markerindex = 1;
    figure
    for i=detvec'
        hold on;
        if ~isempty(apds{i})
            %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
            %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
            if length(apendindices{i}) ~= length(apds{i})
                %                plot(deltat*apds{i},psyms(i,:))
                p(i) = plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
                % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                % stimulus is below peak AP height, (c) BCL is less than initial APD
            else
                %                plot(deltat*apendindices{i},deltat*apds{i},psyms(i,:))
                p(i) = plot(deltat*apendindices{i}, deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:));
                %                plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                %            ylabel('AP duration')
                %                ylabel(['detector at x=' num2str(detloc(i))]);
                %                axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
                axis([0 finaltime 0 400])%axis([0 10 .12 .14])%
            end
            if i == detvec(1)
                ti=title('AP duration vs. time');
            end
        end
        colorcounter = colorcounter + 1;
        if colorcounter > nrco
            colorcounter = 1;
            markerindex = markerindex + 1;
        end

        xl=xlabel('repolarization time (ms)');
        yl=ylabel('APD (ms)');
    end
    legend(num2str(detvec),'Location','EastOutside')
    if poster
        set(gca,'FontSize',fs);
        set(p,'LineWidth',lw1);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
        if numpart == 1
            legend off;
        end
    end


    colorcounter = 1;
    markerindex = 1;
    figure
    for i=1:11%detvec'
        hold on;
        if ~isempty(apds{i})
            %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
            %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
            if length(apendindices{i}) ~= length(apds{i})
                %                plot(deltat*apds{i},psyms(i,:))
                plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                % stimulus is below peak AP height, (c) BCL is less than initial APD
            else
                %                plot(deltat*apendindices{i},deltat*apds{i},psyms(i,:))
                plot(deltat*apendindices{i}, deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                %                plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                %            ylabel('AP duration')
                %                ylabel(['detector at x=' num2str(detloc(i))]);
                axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
            end
            if i == detvec(1)
                title('AP duration vs. time')
            end
        end
        colorcounter = colorcounter + 1;
        if colorcounter > nrco
            colorcounter = 1;
            markerindex = markerindex + 1;
        end

        xlabel('repolarization time')
    end
    legend(num2str((1:i)'),'Location','EastOutside')


    colorcounter = 1;
    markerindex = 1;
    figure
    for i=31:40%detvec'
        hold on;
        if ~isempty(apds{i})
            %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
            %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
            if length(apendindices{i}) ~= length(apds{i})
                %                plot(deltat*apds{i},psyms(i,:))
                plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                % stimulus is below peak AP height, (c) BCL is less than initial APD
            else
                %                plot(deltat*apendindices{i},deltat*apds{i},psyms(i,:))
                plot(deltat*apendindices{i}, deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                %                plot(deltat*apds{i},'Color',colororder(colorcounter,:),'Marker',psyms(markerindex,:))
                %            ylabel('AP duration')
                %                ylabel(['detector at x=' num2str(detloc(i))]);
                axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
            end
            if i == detvec(1)
                title('AP duration vs. time')
            end
        end
        colorcounter = colorcounter + 1;
        if colorcounter > nrco
            colorcounter = 1;
            markerindex = markerindex + 1;
        end

        xlabel('repolarization time')
    end
    legend(num2str((31:i)'),'Location','EastOutside')

    figure
    hold on;
    %        plot(deltat*(1:numstep),V(detlocindices,:),'linewidth',2);
    plot(deltat*(1:32:numstep),V(detlocindices,1:32:numstep),'linewidth',2);
    ylabel('V')
    xlabel('time')
    %        title(['V vs. time: check the phase'])
    title(['V vs. time'])
    %        legend(num2str(detlocindices'))
    legend(num2str(detlocindices'),'Location','EastOutside')


    h4=figure;
    %mesh(deltat*(1:numstep), deltax*(1:numpart), n)
    %    mesh(deltat*(1:2:numstep), deltax*(1:numpart), n(:,1:2:numstep))
    %    mesh(deltat*(1:32:numstep), deltax*(1:numpart), n(:,1:32:numstep))
    %    mesh(deltat*(1:256:numstep), deltax*(1:numpart), n(:,1:256:numstep))
    %    view(90,-90)
    imagesc(deltax*(1:numpart),deltat*(1:256:numstep), n(:,1:256:numstep)')
    axis xy
    %view(120,50)
    colorbar
    %    xlabel('time')
    %    ylabel('distance')
    ylabel('time')
    xlabel('distance')
    title('Contour plot of slow potassium ion-channel variable, n(x,t)')


    if 0
        h5=figure;
        %plot(deltat*apendindex,deltat*apds,'g-o')
        for i=detvec'
            subplot(numdet,1,i)
            if ~isempty(apds{i})
                %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
                %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
                if length(apendindices{i}) ~= length(apds{i})
                    plot(deltat*apds{i},'g-o')
                    % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                    % stimulus is below peak AP height, (c) BCL is less than initial APD
                else
                    plot(deltat*apendindices{i},deltat*apds{i},'g-o')
                    %            ylabel('AP duration')
                    ylabel(['detector at x=' num2str(detloc(i))]);
                    axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
                end
            end
            if i == detvec(1)
                title('AP duration vs. time')
            end
        end
        xlabel('repolarization time')
    end

    if 1

        cellrange=5;
        figure
        hold on;
        %        plot(deltat*(1:numstep),V(detlocindices,:),'linewidth',2);
        plot(deltat*(1:numstep),V(cellrange,1:numstep),'.-');
        ylabel('V')
        xlabel('time')
        %        title(['V vs. time: check the phase'])
        title(['V vs. time'])
        %        legend(num2str(detlocindices'))
        legend(num2str(cellrange),'Location','EastOutside')
        grid

        figure
        hold on;
        %        plot(deltat*(1:numstep),V(detlocindices,:),'linewidth',2);
        plot(deltat*(1:32:numstep),V(1:8,1:32:numstep),'linewidth',2);
        ylabel('V')
        xlabel('time')
        %        title(['V vs. time: check the phase'])
        title(['V vs. time'])
        %        legend(num2str(detlocindices'))
        legend(num2str(detlocindices'),'Location','EastOutside')


        h7=figure;
        hold on;
        %        plot(deltat*(1:numstep),V(detlocindices,:),'linewidth',2);
        %        plot(deltat*(1:numstep),Vd(detlocindices,1:numstep),'--')
        plot(deltat*(1:32:numstep),V(detlocindices,1:32:numstep),'linewidth',2);
        plot(deltat*(1:32:numstep),Vd(detlocindices,1:32:numstep),'--')
        ylabel('V')
        xlabel('time')
        title(['V and Vd vs. time, time step size ' num2str(deltat) ', spatial grid size ' num2str(deltax)])

        if 0
            figure;
            hold on;
            plot(deltat*(1:32:numstep),V(8,1:32:numstep),'linewidth',2);
            plot(deltat*(1:32:numstep),Vd(8,1:32:numstep),'r--')
            ylabel('V')
            xlabel('time')
            title(['V and Vd in cell 8 vs. time, time step size ' num2str(deltat) ', spatial grid size ' num2str(deltax)])

            figure
            for i = 1:numpart
                subplot(numpart,1,i)
                hold on;
                plot(deltat*(1:numstep),V(i,:),'linewidth',2);
                plot(deltat*(1:numstep),Vd(i,1:numstep),'--')
                ylabel(['cell ' num2str(i)])
                xlabel('time')
                axis([0 finaltime -1 8])
                if i == 1
                    title(['V and Vd vs. time'])
                    legend('V','Vd')
                end
            end
        end

        h6=figure;
        hold on;
        %        plot(deltat*(1:numstep),n(detlocindices,:),'linewidth',2);
        %        plot(deltat*(1:numstep),nd(detlocindices(fbdet),1:numstep),'--')
        p1=plot(deltat*(1:32:numstep),n(detlocindices,1:32:numstep),'linewidth',2);
        p2=plot(deltat*(1:32:numstep),nd(detlocindices(fbdet),1:32:numstep),'--')
        xlabel('time')
        ylabel('n')
        legend([p1(1); p2(1)],'n','nd')
        title(['n and nd vs. time, time step size ' num2str(deltat) ', spatial grid size ' num2str(deltax)])

        if 0
            figure
            for i = 1:numpart
                subplot(numpart,1,i)
                hold on;
                plot(deltat*(1:numstep),n(i,:),'linewidth',2);
                plot(deltat*(1:numstep),nd(i,1:numstep),'--')
                ylabel(['cell ' num2str(i)])
                xlabel('time')
                axis([0 finaltime 0 1.5])
                if i == 1
                    title(['n and nd vs. time'])
                    legend('n','nd')
                end
            end
        end

    end


end


if 0%strcmp(fbon,'sf')

    figure
    hold
    plot(svector')
    plot((Vd-V)', ':')
    title('LTV state vector')

    %     figure
    %     hold
    %     plot(min(abs(eigctrb1)))
    %     plot(min(abs(eigctrb2)),'g--')
    %     title('Minimum of ctrb eigenvalue moduli vs. time')


    %     figure
    %     hold
    %     for k=1:numstep
    %         plot(real(eigacl1(:,k)),imag(eigacl1(:,k)),'x')
    %     end
    %     title('Eigenvalues of A-B1*K1')
    %
    %     figure
    %     hold
    %     for k=1:numstep
    %         plot(real(eigacl2(:,k)),imag(eigacl2(:,k)),'x')
    %     end
    %     title('Eigenvalues of A-B2*K2')
    %
    %     figure
    %     hold
    %     for k=1:numstep
    %         plot(real(eigaclred(:,k)),imag(eigaclred(:,k)),'x')
    %     end
    %     title('Eigenvalues of Ared-B1red*Kred and Auc eigenvalues')

    %     figure
    %     hold
    %     plot(min(abs(eigobsv1)))
    %     plot(min(abs(eigobsv2)),'r--')
    %     title('Minimum of observability eigenvalue moduli vs. time')

    %     figure
    %     hold
    %     for k=1:numstep
    %         plot(real(eiga(:,k)),imag(eiga(:,k)),'x')
    %     end
    %     title('Eigenvalues of A')
    %
    %     figure
    %     hold
    %     plot(abs(eigaired'))
    %     title('AIred eigenvalue moduli vs. time')
    if 0
        figure
        hold
        %    plot(abs(eigairedcl'))
        plot(abs(eigacl1'))
        title('AIredcl eigenvalue moduli vs. time')

        figure
        hold
        plot(abs(eiga'))
        title('A eigenvalue moduli vs. time')

        figure
        hold
        plot(kgainswr')

        figure
        hold
        plot(kgains')
    end
    %
    %     figure
    %     hold
    %     plot(min(abs(eigctrbired1)))
    %     plot(min(abs(eigctrbired2)),'r--')
    %     title('Minimum of SFI reduced system OL eigenvalue moduli vs. time')
    %
    %     figure
    %     hold
    %     plot(min(abs(eigctrbi1)))
    %     plot(min(abs(eigctrbi2)),'c:')
    %     title('Minimum of SFI system OL eigenvalue moduli vs. time')
    if 0
        figure
        subplot(2,2,1)
        mesh(deltat*(1:32:numstep), deltax*(1:numpart), V(:,1:32:numstep))
        view(90,-90)
        %view(120,50)
        colorbar
        xlabel('time')
        ylabel('distance')
        title('Contour plot of membrane potential, V(x,t)')

        subplot(2,2,2)
        for i=detvec'
            hold on;
            if ~isempty(apds{i})
                %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
                %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
                if length(apendindices{i}) ~= length(apds{i})
                    plot(deltat*apds{i},psyms(i,:))
                    % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                    % stimulus is below peak AP height, (c) BCL is less than initial APD
                else
                    plot(deltat*apendindices{i},deltat*apds{i},psyms(i,:))
                    %            ylabel('AP duration')
                    %                ylabel(['detector at x=' num2str(detloc(i))]);
                    axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
                end
                if i == detvec(1)
                    title('AP duration vs. time')
                end
            end
            xlabel('repolarization time')
        end
        subplot(2,2,[3:4])
        hold on;
        plot(deltat*(1:numstep),V(detlocindices,:),'linewidth',2);
        ylabel('V')
        xlabel('time')
        title(['V vs. time: check the phase'])
        legend(num2str(detlocindices'))
    end

    if 0
        figure
        for i = 1:numpart
            subplot(numpart+2,1,i)
            hold on;
            plot(deltat*(1:numstep),V(i,:),'linewidth',2);
            plot(deltat*(1:numstep),Vd(i,1:numstep),'--')
            ylabel(['cell ' num2str(i)])
            xlabel('time')
            axis([0 finaltime -1 8])
            if i == 1
                title(['V and Vd vs. time'])
                legend('V','Vd')
            end
        end
        subplot(numpart+2,1,i+1)
        for i=detvec'
            hold on;
            if ~isempty(apds{i})
                %        if stimlocindex == detlocindices(i) & detthresh(i) < stimheight
                %            plot(deltat*apendindices{i}(2:end),deltat*apds{i},'g-o')
                if length(apendindices{i}) ~= length(apds{i})
                    plot(deltat*apds{i},psyms(i,:))
                    % Easier to skip first AP if (a) sensor/actuator are co-located, (b)
                    % stimulus is below peak AP height, (c) BCL is less than initial APD
                else
                    plot(deltat*apendindices{i},deltat*apds{i},psyms(i,:))
                    %            ylabel('AP duration')
                    %                ylabel(['detector at x=' num2str(detloc(i))]);
                    axis([0 finaltime 0 1])%axis([0 10 .12 .14])%
                end
                if i == detvec(1)
                    title('AP duration vs. time')
                end
            end
            xlabel('repolarization time')
        end
        subplot(numpart+2,1,i+2)
        hold on;
        %plot((1:numstep)*deltat,iext([stimlocindex stimoutlocindex],1:numstep))
        % plot((1:numstep)*deltat,iextd([stimlocindex stimoutlocindex],1:numstep),'--')
        % plot((1:numstep)*deltat,iext(:,1:numstep))
        plot((1:numstep)*deltat,iext([stimlocindex],1:numstep),'c')
        plot((1:numstep)*deltat,iextd([stimlocindex],1:numstep),'--')
        %    plot((1:numstep)*deltat,-fb,'g:')
        %    legend('iext_d','iext','-fb')
        legend('iext','iext_d')
        xlabel('time')
        ylabel('iext')
        %    title(['External current vs. time. Current entry location: ' num2str(stimloc) ', exit location: ' num2str(stimoutloc) ', stimulus height ' num2str(stimheight) ', stimulus duration ' num2str(stimduration)])
        title(['External current vs. time. Current entry location: ' num2str(stimloc) ', stimulus height ' num2str(stimheight) ', stimulus duration ' num2str(stimduration)])
    end

end
% mean(rankctr)
% max(rankctr)
% min(rankctr)
%
% id = -iextd;
% save jon_desval Vd nd id
%
%figure(7)
%plot(deltat*(1:numstep),V(1,1:numstep),'m.')
%plot(deltat*(1:numstep),V(2,1:numstep),'g.')
%grid
% proptime_long = 17.4809-17.4015; % propagation time of last long AP, epsil1=0.005
% proptime_short = 18.3614-18.2715;   % propagation time of last short AP, epsil1=0.005
% comment = 'The propagation times (proptimes) of the last 2 APs (long and short) are visual estimates based on where the V values crossed V=2 on the upstrokes';
% save check_prop_speed_D0050_b087 V n iext stimperiod nB M apds apendindices finaltime deltat deltax numstep numpart proptime* comment
% proptime_long = 17.5265-17.4015; % propagation time of last long AP, epsil1=0.0045
% proptime_short = 18.4032-18.2715;   % propagation time of last short AP, epsil1=0.0045
% comment = 'The propagation times (proptimes) of the last 2 APs (long and short) are visual estimates based on where the V values crossed V=2 on the upstrokes';
% save check_prop_speed_D0045_b087 V n iext stimperiod nB M apds apendindices finaltime deltat deltax numstep numpart proptime* comment

%save check_prop_speed_D0053_b087_3cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_3cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_3cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_4cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_5cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_6cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0053_b087_3cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b087_40cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save check_prop_speed_D0090_b097_2cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vunpert  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart

%save Vunpert10  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em6_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em1_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em2_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert5em1_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em0_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1ep1_1_07d  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em1_1_00u  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em1_3_78p  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vpert1em1_1_07d20  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%
%save Vpert1em1_0_93u V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_1_00u V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_1_07u V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_3_78p V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_1_01d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_1_00d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vpert1em1_0_99d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save V5pert1em1_1_01d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save V5pert1em2_1_01d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save V5pert1em3_1_01d V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert
%save Vunpert10b092  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vunpert10b097  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save Vunpert10b112  V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart
%save V5pert1em2_1_02db112 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart xpert tpert

% Save final condition as new initial condition
% note = 'Try to save a final condition corresponding to a stimulus index, so that the signals will be aligned correctly when the simulation is restarted';
% Vinitold = V(:,stimindices(1));
% ninitold = n(:,stimindices(1));
% iextinitold = iext(:,stimindices(1));
% Vinitnew = V(:,stimindices(92)+[-1 0 1]); % save values to the left and right of chosen stimindex in case I am off by one
% ninitnew = n(:,stimindices(92)+[-1 0 1]);
% iextinitnew = iext(:,stimindices(92)+[-1 0 1]);
% save V_b087_time_0_80_sfi KIred1 desfile stimindices stimperiod stimlocindex stimduration stimheight nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart V n iext Vinit* ninit* iextinit* note
% save V_b087_time_0_90_sfi KIred1 desfile stimindices stimperiod stimlocindex stimduration stimheight nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart V n iext Vinit* ninit* iextinit* note
% save V_b087_time_0_90_sf KIred1 desfile stimindices stimperiod stimlocindex stimduration stimheight nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart V n iext Vinit* ninit* iextinit* note
% save V_b087_time_0_90_OL KIred1 desfile stimindices stimperiod stimlocindex stimduration stimheight nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart V n iext Vinit* ninit* iextinit* note

%save 1cellOL_M15_b250 V n iext stimperiod nB M apds apendindices finaltime deltat deltax numstep numpart apds stimduration stimheight
%save 1cellOL_M15_b225 V n iext stimperiod nB M apds apendindices finaltime deltat deltax numstep numpart apds stimduration stimheight
%save 1cellOL_M10_b225 V n iext stimperiod nB M apds apendindices finaltime deltat deltax numstep numpart apds stimduration stimheight
%save check_prop_speed_D0010_M10_b250_20cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0010_M10_b225_10cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0010_M10_b225_2cell V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0010_M10_b225_2cell_shortspike V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0010_M10_b225_2cell_shortspike5 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0001_M10_b225_10cell_shortspike3 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0001_M10_b195_10cell_shortspike3 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0001_M10_b175_10cell_shortspike3 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight
%save check_prop_speed_D0001_M10_b190_10cell_shortspike3 V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart stimduration stimheight

% Vd = V;
% nd = n;
% save 1cellOL_testA_b225 Vd nd
% save 1cellOL_testA_b250 Vd nd
% save 1cellOL_testA_b225_5000 Vd nd
% save 1cellOL_testA_b325 Vd nd
% eval(['load data/' num2str(1)])
% Vd = V;
% nd = n;
% save 1cellOL_testA_b325_3250 Vd nd

%save V1cell_b225_unstable Verror nerror V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart fbon_time KIred1 stimindices iextd stimduration stimheight
%save V1cell_b225_stable_K0_1 Verror nerror V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart fbon_time KIred1 stimindices iextd stimduration stimheight
%save V1cell_b225_stable_K0_02 Verror nerror V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart fbon_time KIred1 stimindices iextd stimduration stimheight
%save V1cell_b225_stable_K0_015 Verror nerror V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart fbon_time KIred1 stimindices iextd stimduration stimheight
%save V_coderewrite_test Verror nerror V n iext stimperiod nB M epsil1 apds apendindices finaltime deltat deltax numstep numpart fbon_time KIred1 stimindices iextd stimduration stimheight

% pause
% fn= input('Enter a filename: ','s')
% eval(['save ' fn]);

%---------------
function out1=karma_f(V,n,nB,M)

global Vstar Vb Vn Vh b ch0 ch1 ch2 ch3
%h=((1-tanh(V-Vh)).*V.^2)/2;
%h=((1-tanh((V-Vh)/ch1)).*(V+ch2).^2)/ch3;
% plot(VV,((1-1*tanh((VV+21)/45)).*(VV+85).^2)/56 + 0.036,'m--') Leave the
% constant offset (0.036) off for now (hoping it doesn't matter much) for
% simplicity. The Jacobians will still be wrong and need to be readjusted
% regardless of the inclusion of the offset, due to other changes in h, f,
% and g. My original motivation for keeping the structure of h is that the
% Jacobians wouldn't have to be changed too much, but it may be easier just
% to switch to Rappel's version of h(V) (see below), since it is easier to
% differentiate than tanh. In any case, the equations on the poster should
% be rewritten.
%
%h=V.^2-delta*V.^3;
%h = 72.83 - .83*V - .04976*V.^2 - (3.52e-4)*V.^3;
h = ch0 + ch1*V + ch2*V.^2 + ch3*V.^3;
% This function seems to produce similar results to the other one, but is
% more sensitive to stimheight and stimduration. I can get some activity
% with height=0.5, duration = 2 or 3, but with the modified h, I need to
% use more current,height 5 and duration 5. I get blowup with height fixed at 0.5 and
% duration above 5 for the former, and heights above 5 for the latter.
%out1= -V + (Vstar - (n/nB)^M)*h;
%out1= -V + (Vstar - (n/nB).^M).*h;
% Include offset from Rappel et al. paper
%out1= -V - 85 + (Vstar - (n/nB).^M).*h;
out1= -V + Vb + (Vstar - (n/nB).^M).*h;


function out2=karma_g(V,n)
global Vstar Vn Vh b ch0 ch1 ch2 ch3

thet = V > Vn;

out2 = (1/b)*thet - n;

%function out=h2(V) % From Rappel's paper
%out = 72.83 - .83*V - .04976*V.^2 - (3.52e-4)*V.^3;