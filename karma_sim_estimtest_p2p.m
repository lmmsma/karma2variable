% Modified karma_sim_estimtest2.m for use with fixed-point estimation.
% Argument = initial condition.
% Added 'load kseparams' to handle batch processing of gains and numbers of cells.
% Assuming the state vector is of the form [V~1 ... V~N ... n~1 ... n~N]',
% a possible C-matrix is [1 0 ... 0; 0 1 0 ...0; ...; 0 ... 1 0 ... 0],
% an rx2N matrix where ones can occur in the first N columns
% (corresponding to cells that have sensing electrodes), and r is in 1:N.
% The observer gain matrix is then 2Nxr, L = [L11 ... L1r; ...; L2N,1 ... L2N,r].
% For coding reasosns, it's easier to assume that C has ones in all
% possible V-locations, then L is 2NxN, with zeros in the appropriate
% columns (those V's that are not fed back). Zeroing out a row in 1:N will prevent
% feedback from going to the V-dynamics, and doing the same in rows
% (N+1):2N will prevent feedback from going to the n-dynamics.
% Then gain L and number of cells (numpart) can be read in from kseparams.
% If L is not loaded, it is set to zero.
% Added L to the list of variables saved in configinfo.
% Added feedback to n via un(:,k) = L((numpart+1):(2*numpart),:)*Verror(:,k);
% Changed trueindex > fbon_time to trueindex >= fbon_time in 3 places.
% Imported comments from ...estimtest2.m:
% This is a multicell estimator code, based off karma_sim_estimtest_1cell.m
% I'm not using karma_sim_estimtest, since I believe it's farther out of
% date. Warning: This still uses the smooth-theta version of the Karma
% model. The justification for this was to make the simulation consistent with the
% single-cell analytical Jacobians. It's not clear whether this is
% necessary or useful for multi-cell simulations, but if Jacobian
% estimation is implemented for larger systems, it may make sense to retain
% it.
%
% Changed epsil1 to 0.001 and stimduration to 7. Note that the n-update law
% is different from that of _statefbtest, since it is a holdover from
% _1cell that is supposed to be consistent with b2=-0.05.
%
% Commented out Lyapunov-test perturbation, since it could accidentally be
% engaged during a sufficiently long simulation. Commented out hprime and
% acoeff* saving commands for now.

% Made configinfo save command consistent with _statefbtest.
% With the default deltax, the conduction speed in the data appears to be
% 11 cells per 1 ms.
% The initial peak to initial peak speed with M=1 (and others at default)
% is 11 cells over 22 ms.
% Increased M to 30. This only seemed to increase the speed to 11 cells /
% 20 ms.
% Increased epsil1 to 0.002 (watch out for Fourier number) and increased
% stimduration to 10 ms to compensate for increase in epsil1, keeping M=1.
% This seemed to produce only a small reduction: 11 cells / 20 ms.
% Increasing stimduration further to 12 ms produced a much larger
% improvement: 11 cells / 10 ms. However, increasing M to 30 gave a better match
% to the upstroke speed than increasing epsil1 and stimduration. The
% downside to increasing M is that it increases the plateau period, whereas
% the measured APs are more triangular than rectangular.
% Increased epsil1 to 0.002 (Fourier number is getting high) and increased
% stimduration to 18 ms. The peak-to-peak speed is 11 cells / 2.5 ms but
% there is large variation in the delay across different parts of the
% upstroke (the maximum gap is around 5 ms). Increasing stimduration is not
% helpful, since the rise time is essentially the same as the stimduration
% for these scenarios, so the peak time is delayed.
% Increased stimheight to (deltat/dtold)*3.2 and decreased stimduration to
% 7. This reduces the rise time while still (apparently) being sufficient
% to trigger APs.
% Increasing M to 10 gives a better match to the rise time, but the
% conduction speed is not reduced much (from 11 cells / 10 ms to 11 cells /
% 7 ms). In addition, the AP shape is fairly square at this value. Notably,
% the increase in M shifts the behavior away from conduction block, to
% decaying alternans. This is more consistent with the data, which shows
% stable 1:1 behavior at 230 ms.
% Decreased epsil2 to 3. This makes the upstroke faster (maybe somewhat
% faster than the measured upstrokes for the second sensor) and increases
% the speed somewhat (to 11 cells / 5 ms).
% For reasons I don't understand, increasing M to 30 doesn't appear to
% affect the conduction speed (it just makes the APs more square), so I'm
% leaving it at 10.
% To reduce Fourier number, increased deltax to 0.05 and reprocessed Vinter
% to match a 37-cell fiber. Sensors 1 and 2 are now assumed to be located
% at cells 1 and 7. For some reason, this seems to reduce the conduction
% speed to 6 cells / 8 ms. This drop could be due to the increased fiber
% length (I was only using 12 cells at 0.0262cm, not the full fiber).
% Decreased epsil2 to 1, which again gives better agreement with the
% measured rise times and reduces the conduction speed to 6 cells / 5 ms.
% Note that the "max" conduction speed for the original deltax is
% 0.262/0.05 = 0.524 cm/ms, and for the larger deltax it is 1.0 cm/ms. The
% measured speed is 0.3 cm/ms. So the goal is to get the speed up to a
% third of the maximum value, although I don't know whether this is
% possible.
% Increased epsil1 to 0.01. The conduction speed doesn't seem to be very
% sensitive to large increases in epsil1. It's down to 6 cells / 3 ms.
% Decreased stimduration to 2.5 ms and increased stimheight to
% (deltat/dtold)*3.2*2.8. This helps increase the speed to 6 cells / 2.5 ms. Further reductions
% in stimduration (and increases in stimheight) may be warranted.
% Replotting V and Vinter with all data points represented reveals a
% mismatch in firing times between cell 1 and detctor 1. If these are to
% match, I shouldn't include the shallow rise period that precedes the main
% depolarization in the measured AP. It also reveals that the conduction
% speed is almost exactly twice as slow as it should be.
% I corrected the upstroke timing in the data to match better with the
% offset in the simulation.
% Raising M to 30 had almost no effect on CV, so left it at 10.
% Decreased stimduration to 1 ms and increased stimheight to
% (deltat/dtold)*3.2*2.8*2. This appears to shift upstrokes backward a bit
% but the CV is similar to before.
% Reduced timestep to 0.025 and adjusted Vinter accordingly. This seems to
% give APs that look different in various ways (which is not a good sign)
% -- faster upstrokes? But the CV is still too low.
% The upstroke timing only matched for the first 3 detectors. Beyond that,
% the simulated CVs were too slow. Various feedback experiments don't seem
% to improve matters very much (local V-based feedback on the corresponding
% V or n channels), and benefits are only local to the feedback location.
% Not sure why this fiber parameter set is harder to control, since
% it seems like increasing the diffusion should couple the cells more
% closely, not make them independent of one another.
% Tried matching AP heights more closely: Decreasing Vstar to 1.2 gives
% better agreement, but 1.0 is too low since it prevents APs from forming.
% Unfortunately, decreasing Vstar also appears to decrease the CV.
% Decreased epsil2 to 0.08, which seems to give a decent match to CV.
% Since high diffusion may be making control difficult, tried reducing
% epsil1 to 0.015 and compensating for the decrease in CV by decreasing
% epsil2 to 0.04. This didn't seem to be a good idea, since epsil2=0.04 introduced weird
% oscillations near the peak value of the simulated V(1,:). Put back the
% previous values of epsil1 and epsil2.
% It may be that the selected data set is not a good choice, since the APDs
% and other aspects (alternans dies out very quickly) may be the "best"
% match for the data, and feedback can do little to correct the details of
% the AP shape. It would be interesting to see whether feedback could
% stabilize alternans (since this was possible with other parameter sets).
%
% Testing other BCLs: alternans still dies out very quickly at 150ms. There
% is alternans transitioning to 2:1 block at 100ms. Not sure where
% alternans emerges in the pacedown (or where the control starts) on 121702
% -- need to search further.
% One thing to try to determine is what "breaks" the ability to control the
% fiber: is it a particular parameter? Can any of the parameters be
% restored to a more "normal" value and still get approximate tracking?
% Could also try testing Pyragas on it to see if it responds.
% The best match I've found is for epsil1=0.03, epsil2=0.08, with deltax =
% 0.05 and deltat=0.025
%
% Problem: for gains over -15 or so (applied to detector 1 or detectors 1 and 2),
% the APs at the feedback locations show "ringing" or buzzing behavior,
% i.e. very
% fast oscillations, at the end of the upstroke. This may be due to numerical instability;
% since the Fourier number is 0.3, feedback may introduce dynamics that are too fast for
% the system? Or it could be that the selected parameters just produce unrealistic behavior
% (could check if it still happens for a much lower Fourier number).
%
% Tests: Try altering epsil1, epsil2, M, and deltax, with feedback applied
% only to cell 1, and gauge the effect on the cell 2 AP. Should avoid
% using feedback at other locations, since the tracking of the timing of
% the measured upstrokes breaks down quickly if any of the parameters are
% altered, whereas the cell 1 upstroke is by default initially aligned with the data,
% and is artificially fast due to the current injection.
% 1) Default parameters with K(1) = -25. Buzzing occurs at end of cell 1
% upstroke.
% 2) Increase epsil2 to 0.3. The simulated CV is now so much slower than that
% of the data, that I can't see how any of the detectors besides the first
% one could reasonably be used for feedback. The buzzing is still present
% (but significantly reduced in its duration) in cells 1 and 2, and the
% effect on cell 2 is significant: its AP now much more strongly represents
% that of cell 1, with a corresponding shift in the n-trajectory. So
% control is "easier" for a higher value of epsil2.
% 3) Decrease epsil1 to 0.008 (this is supposed to be a reduction by the
% same multiplicative factor as epsil2 was increased previously), while
% putting all other parameters at the Test 1 default values. Similarly to Test 2,
% the CV is greatly reduced and causes large misalignment with the measured upstrokes.
% The buzzing is still present in cell 1, and is reduced compared to Test 1, but is somewhat worse
% than that of Test 2 (so increasing epsil2 has a bigger effect on the
% buzzing than decreasing epsil1). As expected, decreasing epsil1 seems to
% make the system "harder" to control: the effect on cell 2 is reduced
% compared to the Test 1 results.
% 4) Increase M to 37.5. The CV is still approximately right, and the
% buzzing is worse (larger amplitude). It looks like the system has become
% somewhat harder to control compared to Test 1 (since the impact on cell 2 seems
% slightly reduced), but it's hard to distinguish this effect from the
% increased mismatch between AP shapes induced by the increase in M.
% 5) Decrease M to 2.6667. The CV is now too slow, but the buzzing is
% significantly reduced compared to Test 1. The ease of control looks
% similar to that of Test 1. The plateau shape is now in better agreement
% with that of the data, but the APDs are way too short. The AP shape is
% only triangular during the plateau, and transitions to a square (very
% steep) drop during repolarization. Hence, the overall match to the AP
% shape is worse than for M=10. I conclude that M doesn't seem to strongly
% affect the control propagation, and M=10 is probably a good value for
% balancing the errors in AP-shape matching (not too square but the APD is
% a decent match).
%
% New parameters: epsil1=0.1, epsil2=1, deltat=0.004, deltax=0.05.
% Fom=0.1600 (<.5) and deltax < space constant =
% sqrt(epsil1*epsil2)=0.3162. With this value of epsil2, the system appears
% to be controllable again, but the value of epsil1 is 100 times higher
% than it should be, and the "true" space constant should be on the order
% of .025 to .1 cm.
% Is the value of deltat OK? reducing it by 1/2 shows a clear increase in
% steepness of the upstroke and downstroke, but I'm not sure whether the
% changes are significant enough to cause worry. It's hard to judge the
% duration of the upstroke: for the proximal cells, it's short, but has a
% noticeable inflection point (reduction in steepness) before the peak
% value. The upstroke duration could be anywhere from .1 ms to 1 ms. The
% "worst case" number of points is then 0.1/deltat = 2 for deltat=0.05
% (this would be far too few intervals for the default deltat, since we want it to be
% "at least 5" intervals per upstroke), and more than enough (25) for
% deltat = 0.004.
% Note: epsil1=0.1 gives a very rough match to the overall conduction time
% from detector 1 to detector 7. The proximal upstroke is faster than that
% of the data, while the remaining simulated upstrokes are much shallower
% than their measured counterparts. Increasing epsil1=0.2 does't really
% help: the CV is increased too much, but the upstroke slope is still much
% too shallow. I haven't figured out yet whether there is a more
% "reasonable" combination of epsil1 and epsil2 that will allow a lower
% epsil1, better matching of upstroke times, and matching of CV, without
% ruining the effectiveness (in terms of impact on adjacent cells) of
% estimation or control feedback.
%
% Upon advice from Robert, changed data scaling. In the model, changed Vb
% from -85 to -90 mV.
% Changed Vn to 65 to compensate for shift in Vb, but it doesn't seem to
% make much difference (just shortens APD slightly by making downstroke
% happen earlier).
% Changed Vstar to 1.3 to increase pulse height, but probably could have
% left it at 1.2. Unfortunately, decreasing Vstar further could be
% desirable (to decrease the pulse height) but decreasing Vstar
% significantly reduces the CV and upstroke slope.
% Try to develop two parameter sets: one where upstroke timing is achieved
% mostly by adjusting Vstar and the other where it is achieved by adjusting
% epsil2.
%
% Another "good" match to the CV occurs with epsil1=0.05, epsil2=0.2, and
% Vstar=1.5. This value of Vstar makes the APs taller than they should be
% compared to the data, but maybe AP height isn't the most important
% feature to get right? "Sacrificing" Vstar allows for larger values of
% epsil2 and smaller values of epsil1 compared to the previous set I found.
% Maximizing epsil2 (for controllability) and minimizing epsil1 (to give epsil1
% a more realistic value) both appear to be desirable. Decreasing epsil1 to
% 0.04 appears to give worse CV agreement overall (comparing both long and
% short beats to the measured values).
% Yet another good match to the CV occurs with epsil1=0.0225, epsil2=0.2,
% and Vstar=1.88 (Rappel's value), although the simulated pulse height is
% higher (and hence more inaccurate) compared to that of the previous parameter
% set (0.05, 0.2, 1.5). I don't think leaving Vstar at the default value is
% necessarily justified, since other parameters such as Vb and Vn were
% altered to better match the data. However, this set uses a lower
% diffusion coefficient (presumably desirable) while still matching the CV.
% Leaving epsil1 at the default value of 0.001, I wasn't able to find a
% sufficiently low value of epsil2 (holding Vstar=1.88) that would allow CV
% matching. For some value of epsil2 above 0.2 but below 1, the APs don't
% propagate out of cell 1 into any of the other cells. Increasing the
% stimulus height didn't appear to help. For epsil2=1, the simulated CV is
% almost 16 times slower than the measured value. I don't think epsil2 can
% be reduced enough to compensate while still preserving AP propagation.
% The (0.05, 0.2, 1.5) set appears to allow slightly better control
% effectiveness than the (0.0225,0.2,1.88) set, although neither is
% particularly good in that respect.
% Another set, with slightly improved control effectiveness (but worse AP
% height matching) over the previous two is (0.04,0.5,2).
% Another set, with very similar control effectiveness (but worse AP
% height matching) compared to the previous one is (0.033,1,2.8). The
% control effectiveness of (0.033,1,2.8) may be slightly better than that
% of (0.04,0.5,2), since it seems to give similar results with a smaller
% feedback signal. That is the feedback signal produced by
% iext(7,k)=-50*Verror(7,k) seems smaller for the newer set compared to
% the previous one. Although the AP height is a worse match, and the
% plateau error is correspondingly larger, the time-avg OL and CL errors
% look like they could be smaller for (0.033,1,2.8) than (0.04,0.5,2). I'm
% not sure how important it is to get the upstroke speed correct (how right
% does it have to be to get multisite feedback to work).
%
% Commented out Vd resizing commands since they aren't necessary for the
% estimator setup.

% Older comments:
% This is a restart of the single-cell estimator code, based off the latest
% version of the statefbtest_1cell code (from version 17).

% This is almost the same as version 15 of ..._statefbtest.m.
% Differences should include the definition of Bn, which was changed to match the "size" of B1,
% as well as the structure of the n_t update law? .
% YET ANOTHER ERROR in acoeff's: forgot the factor of 1/nB in acoeff2.
% Since nB=1 in Rappel's units, it shouldn't change the answer, but it
% affected previous versions in the old unit system.

% Older comments: 
% To prevent "wraparound" reading of Vd and nd,
% 1) disable line around 573: if writeintsteps <= size(Vd,2) should be if 0
% 2) revise lines around line 885 (acoeff1, etc.) to use trueindex instead
% of k (except for hprime ref in acoeff1)

% Rewrote to avoid memory errors for large arrays. Changed function
% argument to bcl.

%function [] = karma_sim_statefbtest(bcl,varargin)
function xout = karma_sim_statefbtest(x)%p2p
%function altamp = karma_sim_statefbtest(x)%p2p
% simulate 2-variable Karma model for 1-D geometry

global Vstar Vb Vn Vh b ch0 ch1 ch2 ch3 capprox

mkdir('data');
%dfl = input('Delete old files? y or n: ','s')
%if dfl=='y'%p2p
delete data/*.mat;
%end

if exist('kseparams.mat')
    load kseparams.mat % load gain L, number of cells, bcl
else
    disp('Parameter file does not exist')
end

%L = zeros(212,106); 
%L(16,16) = -1;
%L(16,16)=  -0.009278521180909; 
%L(106+16,16)=   0.000721478819091; 
%L(16,16)=  -10; 
%L(106+16,16)=   0.0001; 

if ~exist('bcl') % p2p; if it was not loaded from kseparams
    bcl = 230 %p2p
end

karma_test = 0;
pubyear = 'ext';
dataflag = 0%'zoh'; % Use zero-order hold if set to 'zoh', otherwise use linearly interpolated data
%writeint = round(1000/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(3250/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(325/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(500/bcl)*bcl; % maximum array length in ms, 40 cells
%writeint = round(250/bcl)*bcl; % maximum array length in ms, 80 cells
writeint = bcl; %p2p
%writeint = 2*bcl; %p2p
% Confined to be an even number of bcls to help with desired value alignment
% Chosen as an integer number of bcls closest to some arbitrary duration

if pubyear == 93
    deltat = 0.03/16; %ms (deltat not given in paper, just made it small)
elseif pubyear == 94
    %deltat = 0.25; %ms
    deltat = 0.03/16;%0.03/(16*7);%0.25/32; %ms
elseif pubyear == 'ext'
    deltat = 0.008; %ms
    dtold = deltat;
    dtscale = 1;%4; % integer scaling factor
    deltat = deltat*dtscale; %ms
end
%deltax = .0262; %cm
%deltax = .05; %cm
%deltax = .01; %cm
deltax = .02; %cm
%deltax = .1; %cm

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
    if ~exist('numpart') % p2p; if it was not loaded from kseparams
        fiberlength = 1*deltax;%p2p
    else
        fiberlength = numpart*deltax;
    end
    %        fiberlength = 2*deltax;%p2p
end
%finaltime = 4200; %ms
%finaltime = 1150; 
%finaltime = 19*bcl;%p2p
%finaltime = 1*bcl;%p2p
%finaltime = 2*bcl;%p2p
finaltime = 5*bcl;%p2p

% derived quantities
numrep = ceil(finaltime/writeint); % number of cycles based on finaltime and writeint
writeintsteps = round(writeint/deltat);
if ~exist('numpart') % p2p; if it was not loaded from kseparams
    numpart = round(fiberlength/deltax); % number of partitions
end
numstep = round(finaltime/deltat); % time steps
%finalblocksteps = numstep-writeintsteps; % length of final interval;
finalblocksteps = writeintsteps - (numrep*writeintsteps-numstep); % length of final interval;
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
    %    tauN = 250;
    %    tauN = 200;
    tauN = 170;
    Dg = 0.0108;
    De = 5*Dg; % De=1/(c*re), inversely proportional to external fluid resistance
    c = 0.009;
    %    cc_ext = max([0 (Dg/(c*(Dg + De)))]);
    %    cc_int = max([0 (De/(c*(Dg + De)))]);
    cc_ext = 1; % uF/cm^2 ? Rappel, et al.
    cc_int = 1;
    if karma_test
        M=30;
        nB=0.525;
    else
        M=10; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=37.5; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=.08*10/.3; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=1 % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=30 % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        nB=0.707; %0.680;% increasing makes more "sensitive"
        nB=1;
    end
    %    b = 1;
    %    b = 1/(tauN*(1+2.4)/2.4);
    b = 1/((1+2.4)/2.4); % I think there is a typo at the end of Rappel's paper, where all the terms in g(V,n) should be
    % divided by tauN. If so, this should be the correct expression for b.
    % If I use the other value of b implied by the paper (which includes
    % tauN), then the resulting n trajectories look wrong (too steep due to
    % domination of the step function) and I wasn't able to find stimulus
    % values that would allow the cell to fire "normally".
    %epsil = .009; % increasing too much means no AP formation; some runtime error if pushed below .007
    if numpart > 1
        epsil1 = 0.005; % cm^2/ms
    else
        epsil1 = 0 %p2p
    end
    %    epsil1 = 0.2; % cm^2/ms % with epsil2=1, the propagation time to the distal end is shorter than the measured time (CV too large)
    %    epsil1 = 0.0225; % cm^2/ms % suitable for epsil2=0.2 and Vstar=1.88
    %    epsil1 = 0.001*0.75; % cm^2/ms at 225, dur=7, alt is slightly worse but not discordant
    %    epsil1 = 0.001*0.5; % cm^2/ms at 225, dur=7, alt is not discordant (may be slightly less severe than 0.75 case? why?)
    %    epsil1 = 0.001*0.25; % cm^2/ms at 225, dur=7, alt is not discordant (but more spatial variation than 0.5 case, and cell 1 looks more stable). Looks similar for dur=5 and 3, except alt. amp is somewhat reduced
    epsil2 = 0.7;
    % Fourier number
    Fom = epsil1*deltat/(deltax)^2
end

Vstar = 4.0; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
Vn = -65; % lowering seems to make alternans go away faster; making too high causes breakup/block
Vh = -21; % lowering seems to make alternans go away faster, effect of increasing is not clear
%Vb = -85; % mV (offset from Rappel)
Vb = -90; % mV 
% coefficients in h(V)
%h = 72.83 - .83*V - .04976*V.^2 - (3.52e-4)*V.^3;
%h=((1-tanh((V-Vh)/ch1)).*(V+ch2).^2)/ch3;
% Redefine ch* to match power series coefficients 1/06/09
%   ch0 = 72.83;
%   ch1 = -0.83;
%   ch2 = -0.04976;
ch3 = -3.52e-4;
% for roots(f) = -86, -80, 17.8 (Vstar = 2.8)
%ch0 = 76.8441;
%ch1 = -1.0177;
%ch2 = -0.0523;
% for roots = -89, -80, 17
% ch0 = 77.3083;
% ch1 =  -1.1269;
% ch2 = -0.0537;
% for roots = -90, -85, 17
%ch0 = 80.0510;
%ch1 = -1.2548;
%ch2 = -0.0555;
% for roots = -90, -85, 25
% ch0 = 100.7446;
% ch1 = -0.7814;
% ch2 = -0.0528;
% for roots = -90, -80, 20 (Vstar = ?)
%ch0 =  84.1126;
%ch1 = -0.9662;
%ch2 = -0.0528;
% for roots = -90, -65, 20 (Vstar = 2.5)
% ch0 =  72.7228;
% ch1 = -0.7760;
% ch2 = -0.0493;
% for roots = -90, -80, 20 (Vstar = 2.5)
% ch0 = 92.5228;
% ch1 = -0.8728;
% ch2 = -0.0528;
% for roots = -90, -80, 25 (Vstar = 2.5)
%ch0 = 105.1948;
%ch1 = -0.5736;
%ch2 = -0.0510;
% for roots = -90, -80, 26 (Vstar = 2.5)
%ch0 = 107.7292;
%ch1 = -0.5137;
%ch2 = -0.0507;
% for roots = -90, -80, 24 (Vstar = 2.5)
% ch0 = 102.6604;
% ch1 = -0.6334;
% ch2 = -0.0514;
% for roots = -90, -80, 32 (Vstar = 3.0, n=.8)
% ch0 = 115.2557;
% ch1 = -0.3351;
% ch2 = -0.0496;
% for roots = -90, -87, 27 (Vstar = 3.0, n=.8)
% ch0 = 105.5299;
% ch1 = -0.7282;
% ch2 = -0.0528;
% for roots = -90, -70, 28 (Vstar = 3.6, n=.8)
%  ch0 = 87.8614;
%  ch1 = -0.3543;
%  ch2 = -0.0465;
% for roots = -90, -72, 35 (Vstar = 4.5, n=.8)
% ch0 = 100.3225;
% ch1 = -0.0575;
% ch2 = -0.0447;
% for roots = -90, -80, 35 (Vstar = 5.0, n=.8)
%  ch0 = 107.0990;
%  ch1 = -0.2356;
%  ch2 = -0.0475;
%for roots = -90, -76, 35 (Vstar = 5.0, n=.8)
%ch0 = 102.6638;
%ch1 = -0.1582;
%ch2 = -0.0461;
% for roots = -90, -65, 35 (Vstar = 1.88, n=.6)
%  ch0 = 120.0988;
%  ch1 = 0.3840;
%  ch2 = -0.0422;
% for roots = -90, -75, 35 (Vstar = 1.88, n=.6)
%  ch0 = 131.1868;
%  ch1 = 0.1904;
%  ch2 = -0.0458;
% %for roots = -90, -65, 35 (Vstar = 5.0, n=1.02)
% ch0 = 95.8752;
% ch1 = 0.1149;
% ch2 = -0.0422;
%for roots = -90, -75, 35 (Vstar = 2.8, 0.7)
%   ch0 = 115.6304;
%   ch1 = 0.0176;
%   ch2 = -0.0458;
%for roots = -90, -65, 35 (Vstar = 4, 0.7)
%  ch0 = 94.5775;
%  ch1 = 0.1005;
%  ch2 = -0.0422;
%for roots = -90, -65, 35 (Vstar = 4, 0.7)
ch0 = 107.8831;
ch1 = -0.1319;
ch2 = -0.0465;

deton = 1; % turn on AP detectors
%detloc = [0.1 0.25 0.40]; % location of measurement electrodes, measured
%from upstream end (x=0) (cm)
detloc = deltax*(1:numpart); % location of measurement electrodes, measured from upstream end (x=0) (cm)
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
%fbdet = 1% index of detector used for feedback measurments
stimoutdet = 4; % index of detector used to measure downstream current extraction point

stimloc = deltax;%0.1 % location of stimulating electrode, measured from upstream end (x=0) (cm)
stimheight = (deltat/dtold)*3.2*2.8; %1.5/5;%
%stimheight = 1.6; %(deltat/dtold)*0.2;%mV %Rappel et al., 10-20mV? Using Rappel's h(V), the "minimum" stimulus for
% a clean-looking AP seems to be either 0.2 height for 1 ms or 0.1 for 2
% ms, when deltat=0.05/8. Stimheight should go up by a factor of deltat/dtold if
% deltat is increased. For the old approximate h(V), I expect it to be higher, but I haven't
% checked yet. The units of iext as defined in the code should be uF/cm^2
%stimduration = 1/250; %deltat% assumptions: 1) period = 1 corresponds to 250ms, 2) typical stimulus duration is 1-5ms
stimduration = 1;
%stimoutloc = 0.2%2*deltax; %2.5;% location where current is extracted
if ~exist('stimstart') % p2p; if it was not loaded from kseparams
    stimstart = deltat; % start time for stimulation (ms)
end
stimperiod = bcl*[1 1];
% For Rappel, the critical point for 1 cell appears to be between 225 and
% 250 ms for stimheight = 1.6, dur=1

stimtransamt = 0.5; % step amount (ms) by which period is increased or decreased during transition

detlocindices = round(detloc/deltax); % cell indices of measurement electrode locations
stimdurationsteps = round(stimduration/deltat);
stimlocindex = round(stimloc/deltax); % cell index of stimulating electrode location
%stimoutlocindex = round(stimoutloc/deltax); % cell index of location where current is extracted
stimstartindex = round(stimstart/deltat); % time index for initial stimulus
stimperiodsteps = round(stimperiod/deltat); % timesteps between stimuli
%if ~exist('icfile')
%        fbon_time = round(450/deltat); %round(25/deltat);% time index at which control is switched on
    fbon_time = 1; %p2p
%else
%    fbon_time = 1; % time index at which control is switched on
%end

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

V = zeros(numpart,writeintsteps);
% if karma_test
%     V(1,1) = 2; % initiate pulse that will travel around ring.
%     % If we want to change the initiation point, then the periodic BCs must
%     % change also. Right now, a pulse is started at the leftmost end, after
%     % which the ends of the fiber are tied together: V_dn = V(end,...),
%     % and Vup = Vdn
% end
n = zeros(numpart,writeintsteps);
%V(:,1) = -85; % initial condition
%n(:,1) = 0; % initial condition
V(:,1) = x(1:numpart,:); %p2p
n(:,1) = x((numpart+1):end,:); %p2p
% Load ICs from SFI run
if exist('icfile')
    eval(['load ' icfile ' Vinitnew ninitnew iextinitnew'])
    V(:,1) = Vinitnew(:,2);
    n(:,1) = ninitnew(:,2);
end

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
iextd = zeros(numpart,writeintsteps); % desired external current
ctr = 1;
if cc_ext
    while stimindices(ctr) < writeintsteps%numstep
        iextd(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -stimheight/deltat/cc_ext; % desired external current
        if ~fbon %| strcmp(fbon,'sf')
            iext(stimlocindex,stimindices(ctr):(stimindices(ctr)+stimdurationsteps)) = -stimheight/deltat/cc_ext; % desired external current
        end
        ctr = ctr+1;
    end
end

if 1%epsil1==0.009
    if dataflag ~= 'zoh'
        %    desfile = '121702control1c_proc_50_52_43cells_deltax0p05_deltat0p008';
        %desfile = 'desval_1cell_b230_smooththeta_nsoli';
        %    desfile = 'desval_2cellOL_b230_smooththeta_nsoli';
        %   desfile = 'desval_1cell_b160_2to1_smooththeta_nsoli';
        %    desfile = 'desval_2cellOL_b160_2to1_smooththeta_nsoli';
        if bcl == 230
            desfile = 'desval_106cellOL_b230_smooththeta_nsoli';
        elseif bcl == 200
            desfile = 'desval_106cellOL_b200_smooththeta_nsoli';
        end
%            desfile = 'desval_106cellOL_b200_smooththeta_nsoli20';
        %    eval(['load ' desfile ' Vinter']);
%        eval(['load ' desfile ]);%p2p
        load(desfile); 
        %    Vd = Vinter;
        Vd = Vinter(:,1:(writeintsteps+1));%p2p
        %      Vinter=zeros(numpart, finaltime/deltat + writeintsteps+1); % use if
        %      needed for a long OL simulation (e.g. trying to find long-term
        %      pattern)
    else
        %    desfile = '121702control1c_proc_50_55_43cells_deltax0p05_deltat1';
        %    desfile = '121702control1c_proc_50_55_211cells_deltax0p01_deltat1';
        desfile = '121702control1c_proc_50_55_106cells_deltax0p02_deltat1';
        eval(['load ' desfile ' Vinter dtact']);
    end
    %    Vinter = Vinter(1:numpart,:); %p2p
    % Comment out the previous line and comment in the following 3 lines if the 1:1 or data trajectories have
    % too few cells (e.g., before the fixed point has been found)
    %     Vinter = zeros(numpart,(writeintsteps+1)); %p2p
    %     Vd = Vinter; %p2p
    %     nd = Vinter; %p2p
    
    %    Vd(:,1) = Vinter(:,1); %p2p
    if ~exist('nd')
        nd = zeros(size(Vd));
    end
    if ~exist('dapd')
        dapd = -1;
    end
end

% state feedback term
usf = zeros(1,writeintsteps);
if exist('icfile')
    iext(:,1) = iextinitnew(:,2);
    usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
end
nerror = nd(:,1:writeintsteps)-n;
Verror = Vd(:,1:writeintsteps)-V;
KIred1=zeros(1,2*numpart);

capprox = 1; % Gain used for approximation of Heaviside function 1(V-Vn) as 0.5*(1+tanh(capprox*(V-Vn)))
% hprime = zeros(numpart,numstep);
% acoeff1 = zeros(numpart,numstep);
% acoeff2 = zeros(numpart,numstep);
% acoeff3 = zeros(numpart,numstep);
% %s=zeros(2,numstep);

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

perttime = 19900; %time of perturbation for Lyapunov test
%pertval = -10; % reset value in mV
pertval = 0; % reset value in mV

xout = zeros(2*numpart,1);
% ff=figure;
% hold on;

if dataflag == 'zoh'
    vdctr = 1; % Counter for ZOH on data
end

Vtmax = 120 %mV; approximate span of an AP
ntmax = 1.1; %dimensionless; approximate max of n variable
u1max = 97; %muA/cm^2; max for current injections (based on estimate of max usf value for "smallest" stabilizing gain)
u2max = 0.065; %dimensionless; approximate max of n-perturbations (based on estimate of max n-pert value for "smallest" stabilizing gain)

Smat = diag([1/Vtmax 1/ntmax]); % scaling matrix for error states
Smatinv = inv(Smat);
Tscale = 1; % scaling for time derivatives. For cts time it should have units of time, but for discrete time, it should be dimensionless

tic
for ii=1:numrep % number of data writing cycles
    if ii > 1
        V = zeros(numpart,writeintsteps);
        if karma_test
            V(1,1) = 2; % initiate pulse that will travel around ring.
            % If we want to change the initiation point, then the periodic BCs must
            % change also. Right now, a pulse is started at the leftmost end, after
            % which the ends of the fiber are tied together: V_dn = V(end,...),
            % and Vup = Vdn
        end
        n = zeros(numpart,writeintsteps);
        if dataflag ~= 'zoh'
            Vd = Vinter(:,((ii-1)*writeintsteps+1):(ii*writeintsteps+1));
            nerror = nd(:,1:writeintsteps)-n;
            Verror = Vd(:,1:writeintsteps)-V;
        else
            nerror = zeros(numpart,writeintsteps);
            Verror = zeros(numpart,writeintsteps);
            Vd(:,1) = Vinter(:,vdctr);
        end
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
    %s=zeros(2,numstep);
    
    if ii==numrep
        wi = finalblocksteps;
    else
        wi = writeintsteps;
    end
    for k=1:wi %k = 1:numstep-1
        trueindex = (ii-1)*writeintsteps + k;
        if strcmp(fbon,'sf')
            if 0
                acoeff1(:,k) = (-1+(Vstar - (n(:,k)/nB).^M).*hprime(:,k))/epsil2;
                acoeff2(:,k) = -(M*(1/nB)*(n(:,k)/nB).^(M-1)).*htemp/epsil2;
                acoeff3(:,k) = 0.5*(1/b)*(capprox*sech(capprox*(V(:,k)-Vn)).^2)/tauN;
                if numpart ==1
                    Ac = [(-deltat*epsil1)/deltax^2+deltat*acoeff1(1,k)+1];
                    A12 = [deltat*acoeff2(1,k)];
                    %        A21 = [0];
                    Auc = [1-deltat/tauN];
                    %%%%%%%%%
                    % Alternative: use "true" A matrix
                    A21 = [deltat*acoeff3(1,k)];
                    A = [ Ac   A12
                        A21  Auc];
                end
            end
            
            if trueindex >= fbon_time %***k > fbon_time %& abs(V(stimlocindex,k)) < 4.5
                sext=[Verror(:,k); Verror(:,k)-Verrorprev];
                KIred1 = zeros(1,2*numpart);
                
                fac =0.999;
                %                fac =0.998; % this value doesn't seem to improve the bad transients with -0.01, -0.005
                %                fac =0.997; % this value doesn't seem to improve the bad transients with -0.01
                %                usf(k) = usf(k-1)-KIred1*sext;
                % ***               sext2=[Verror(:,k); Verror(:,k)-fac*Verror(:,k-1)];
                sext2=[Verror(:,k); Verror(:,k)-fac*Verrorprev];
                %                    sext2=[nerror(:,k); nerror(:,k)-fac*nerrorprev]; % use this to apply n-based feedback to the V-dynamics
                % ***               usf(k) = 1*(fac*usf(k-1)-KIred1*sext2);
                usf(k) = 1*(fac*usfprev-KIred1*sext2);
            end
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
            if isnan(usf(k))
                usf(k) = 0;
            end
            
%            if trueindex >= fbon_time
            if trueindex >= fbon_time %& (mod(trueindex,stimperiodsteps(1))<=99) 
                
                if exist('L') %p2p
                    uppergaintimesve = L(1:numpart,:)*Verror(:,k);
                    usf(k) = -uppergaintimesve(1); % usf uses opposite sign convention.
                    % usf being a separate variable made more sense in the boundary-control
                    % study, where all feedback terms were applied to cell 1. For
                    % distributed feedback (i.e. the observer) it makes less sense.
                    if numpart > 1
                        iext(2:numpart,k) = uppergaintimesve(2:numpart);
                    end
                else
                    L = 0;
                end
            end
            iext([stimlocindex],k) = iextd(stimlocindex,k)-usf(k);
            % To move usf to cell 2, comment out the previous line and
            % uncomment the next 3, and alter KIred so it only has nonzero gains on cells 2:end
            %         iext([stimlocindexold],k) = iextd(stimlocindexold,k);
            %         stimlocindex = stimlocindexold+1;
            %         iext([stimlocindex],k) = -usf(k);
            
            if abs(iext(stimlocindex,k)) > stimheight/deltat/cc_ext
                iext(stimlocindex,k) = stimheight/deltat/cc_ext*sign(iext(stimlocindex,k));
                %            usf(k) = 0;
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
        
        
        if pubyear ~= 'ext' & ~karma_test & sum(k == stimindices_new) 
            % stimulate periodically, at designated location
            %        stimtrue = sum(k==stimindices);
            %        V(stimlocindex,k) = stimtrue*2 + (1-stimtrue)*V(stimlocindex,k);
            V(stimlocindex,k) = stimheight;
            % ??? What's the activation threshold ??? I thought it should be
            % V=Vn=1, but it appears that somewhat lower values work, also?
        end
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
        
        if trueindex >= fbon_time % This set should be consistent with b2 = -0.05
            %                un(:,k) = 0.005*Verror(:,k)/tauN; % proportional feedback term applied to n instead of V
            %                gains in foregoing line (using Verror for fb): negative
            %                gains, such as -0.005/tauN, don't appear to work.
            %                However, 0.005*deltat/tauN stabilizes, 0.003 stabilizes slowly, and 0.001 doesn't.
            
%            if exist('L') & size(L,1) > numpart %p2p
            if exist('L') & size(L,1) > numpart %& (mod(trueindex,stimperiodsteps(1)) <=99) %p2p
                un(:,k) = L((numpart+1):(2*numpart),:)*Verror(:,k); % assume that only Verror is available for feedback
            end
        end
        %        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k); % This one is consistent with b2=1 assumption
        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN - deltat*un(:,k); % This should be consistent with b2 = -0.05
        
        if dataflag == 'zoh' & k > 1 & ~mod(k,1000*dtact/deltat)
            vdctr=vdctr+1;
        end
        
        if dataflag ~= 'zoh'
            %         Verrornext = Vd(:,k+1)-Vnext;
            Verrornext = Vinter(:,k+1)-Vnext;
        else
            Verrornext = Vinter(:,vdctr)-Vnext;
        end
        nerrornext = nd(:,k+1)-nnext;
        
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
            end % if falling edge
        end % if deton
        %         if ~mod(trueindex,1000) % *** ~mod(k,1000)%p2p
        %             disp(trueindex) % *** disp(k)
        %         end
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
    %    eval(['save data/' num2str(ii) ' un hprime acoeff* usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
    checkii = ii
    datafname = ['data/' num2str(ii)];
%    eval(['save data/' num2str(ii) ' un usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
    save(datafname,'un', 'usf', 'Verror', 'nerror', 'iext', 'V', 'n', 'apds', 'apendindices', 'dctr', 'rctr', '*prev', '*next');  
    clear hprime acoeff* usf Verror nerror iext V n un
    if ii == numrep
        %        save data/configinfo stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength stimlocindex dapd desfile capprox stimstart
        save data/configinfo stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength stimlocindex dapd desfile capprox stimstart
        if exist('L')
            save data/configinfo L -append
        end
        xout = [Vnext; nnext];
        %altamp = diff(deltat*apds{end}); 
%         tempout = zeros(size([Vnext;nnext]));
%         if apds{numpart/2} == 2
%             if abs(diff(apds{numpart/2})*deltat) < 1
%         tempout = 1e-2*ones(size(tempout)); 
%             end
%         else
%         tempout = 1e-1*ones(size(tempout)); 
%         end
%         xout = [Vnext; nnext]+tempout;
    end
end % data write loop
toc
%k
if 1
    %        crange = [1 7 13 19 25 31]; %detlocindices
    %        crange = [1 4 7 10 13 16]; %detlocindices
    %        crange = [1 7 13 19 25 31 37]; %detlocindices
    %        crange = [37]; %detlocindices
    %        crange = [7 13 19 25 31 37]; %detlocindices
    %        crange = [31 61 91 121 151 181]; %detlocindices
            crange = [16 31 46 61 76 91]; %detlocindices
    %crange = [1:numpart]; %p2p
       % Define an offset colororder
             hfig = figure
        colororder=get(gca,'ColorOrder');
        close(hfig);
        colororderoffset = colororder([3:7 1:2],:);
        % Plot V and n
        h=figure;
        clear p1 p2;
        for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
            subplot(2,1,1)
            hold on;
            eval(['load data/' num2str(ii)])
            %    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
            set(gca,'ColorOrder',colororder);
            %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(detlocindices,1:writeintsteps));
            p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange,1:50:writeintsteps));
            %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(crange,1:writeintsteps));
            set(gca,'ColorOrder',colororderoffset);
            if dataflag ~= 'zoh'
                p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Vinter(crange,[((ii-1)*writeintsteps+1):ii*writeintsteps]),'--');
            else
                p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Verror(crange,1:50:writeintsteps)+V(crange,1:50:writeintsteps),'--');
            end
            clear hprime acoeff* usf Verror nerror iext V n un
        end
        %    yl=ylabel('membrane potential (mV)');
        yl=ylabel({'membrane', 'potential (mV)'});
        ti=title(['V and V^d vs. time']);
        legend([p1(end); p2(end)],'V','V^d')
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
        yl=ylabel({'refractory', 'variable'});
        ti=title(['n and n^d vs. time']);
        %    legend([p1(end); p2(end)],'n_i','n^d_i')
        legend([p1(end); p2(end)],'n','n^d')
        axis([0 finaltime 0 1.2])
        if poster
            set(gca,'FontSize',fs);
            set(p1,'LineWidth',lw1);
            set(p2,'LineWidth',lw2);
            set(xl,'FontSize',fs);
            set(yl,'FontSize',fs);
            set(ti,'FontSize',fs);
        end
    
    %Plot current and feedback terms
    h=figure;
    clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange],1:writeintsteps));
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    yl=ylabel('current (\muA/cm^2)');
    ti=title(['Stimulus and feedback terms vs. time']);
    %    axis([0 finaltime -20 20])
    axis([0 finaltime -1200 1200])
    %    axis([0 finaltime -6e-5 6e-5])
    %    axis([0 finaltime -2e-3 2e-3])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    set(gca,'Position',[0.13 0.583837 0.775 0.339794])
    if ~exist('psyms')
        hfig = figure
        colororder=get(gca,'ColorOrder');
        close(hfig);
        nrco = size(colororder,1); % number of rows (distinct colors)
        
        %    psyms = strvcat('g-o','--*','m-s','k-x','r-+','y-h','c-^','g-v','--x','m-+','k-h');
        %    psyms = strvcat('+','o','*','x','s','d','^','v','>','<','p','h');    %same order as in Matlab help file
        psyms = repmat(strvcat('o','*','s','+','x','d','^','v','>','<','p','h'),2,1);
    end
    
    clear p1;
    for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        %    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Verror(crange,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Verror(crange,1:50:writeintsteps));
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    %    yl=ylabel('membrane potential (mV)');
    yl=ylabel({'membrane', 'potential (mV)'});
    %    ti=title(['V and V^d vs. time']);
    ti=title(['Verror']);
    %    legend([p1(end); p2(end)],'V_i','V^d_i')
    axis([0 finaltime -100 100])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        %       set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    
    % Plot un and nerror
    h=figure;
    clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([crange],1:writeintsteps));
        %        plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,usf(1,1:writeintsteps),'g:');
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    yl=ylabel('perturbation to n_{ t}');
    ti=title(['Feedback term vs. time']);
    axis([0 finaltime -1.2 1.2])
    %    axis([0 finaltime -1200 1200])
    %    axis([0 finaltime -6e-5 6e-5])
    %    axis([0 finaltime -2e-3 2e-3])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    set(gca,'Position',[0.13 0.583837 0.775 0.339794])
    if ~exist('psyms')
        hfig = figure
        colororder=get(gca,'ColorOrder');
        close(hfig);
        nrco = size(colororder,1); % number of rows (distinct colors)
        
        %    psyms = strvcat('g-o','--*','m-s','k-x','r-+','y-h','c-^','g-v','--x','m-+','k-h');
        %    psyms = strvcat('+','o','*','x','s','d','^','v','>','<','p','h');    %same order as in Matlab help file
        psyms = repmat(strvcat('o','*','s','+','x','d','^','v','>','<','p','h'),2,1);
    end
    
    clear p1;
    for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        %    p1=plot(deltat*(1:1:numstep),V(detlocindices,1:1:numstep));
        set(gca,'ColorOrder',colororder);
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(detlocindices,1:writeintsteps));
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange,1:50:writeintsteps));
        %        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Verror(crange,1:writeintsteps));
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,nerror(crange,1:50:writeintsteps));
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    %    yl=ylabel('membrane potential (mV)');
    yl=ylabel({'refractory', 'variable error'});
    %    ti=title(['V and V^d vs. time']);
    ti=title(['nerror']);
    %    legend([p1(end); p2(end)],'V_i','V^d_i')
    axis([0 finaltime -1.2 1.2])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        %       set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    
 
    % Plot scaled values

    h0=figure;
    clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange],1:writeintsteps)/u1max);
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    yl=ylabel('scaled current');
    ti=title(['Scaled stimulus and feedback terms vs. time']);
%    axis([0 finaltime -1200 1200])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    set(gca,'Position',[0.13 0.583837 0.775 0.339794])
      clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([crange],1:writeintsteps)/u2max);
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    yl=ylabel('perturbation to n_{ t}');
    ti=title(['Scaled feedback term vs. time']);
%    axis([0 finaltime -1.2 1.2])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    
    h11=figure;
    clear p1 p2;
    for ii=1:numrep%(numrep-10):(numrep-1)%1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        set(gca,'ColorOrder',colororder);
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,V(crange,1:50:writeintsteps)/Vtmax);
        set(gca,'ColorOrder',colororderoffset);
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    yl=ylabel({'membrane', 'potential (scaled)'});
    ti=title(['V/Vtmax vs. time']);
%    axis([0 finaltime -100 40])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
         set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,2)
        hold on;
        eval(['load data/' num2str(ii)])
        set(gca,'ColorOrder',colororder);
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,n(crange,1:50:writeintsteps)/ntmax);
        set(gca,'ColorOrder',colororderoffset);
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    xl=xlabel('time (ms)');
    yl=ylabel({'refractory', 'variable'});
    ti=title(['n/ntmax vs. time']);
%    axis([0 finaltime 0 1.2])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(xl,'FontSize',fs);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end

    
    hh=figure;
    clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        %       p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps)*(-9/.0051),'g--');
        %       p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps)*(-9/.003),'g--');
        %            p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps),'--');
        %            p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un(1,1:writeintsteps),'--');
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange],1:writeintsteps));
        clear hprime acoeff* usf Verror nerror iext V n un
    end
    %        yl=ylabel('perturbation to n_{ t}');
    yl=ylabel('current (\muA/cm^2)');
    ti=title(['Stimulus and feedback terms vs. time']);
    %    axis([0 finaltime -20 20])
    axis([0 finaltime -1200 1200])
    %    axis([0 finaltime -6e-5 6e-5])
    %    axis([0 finaltime -2e-3 2e-3])
    if poster
        set(gca,'FontSize',fs);
        set(p1,'LineWidth',lw1);
        set(yl,'FontSize',fs);
        set(ti,'FontSize',fs);
    end
    set(gca,'Position',[0.13 0.583837 0.775 0.339794])
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
                axis([0 finaltime 0 300])%axis([0 10 .12 .14])%
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
    if ~isempty(apds{1})
    p(end+1)=plot(deltat*apendindices{i}, dapd*ones(size(apendindices{i})),'r--');
    end
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
    
    %al = (Vstar-(1.02)^M)
    %al = (Vstar-(.7)^M)
    al = (Vstar-(.5)^M)
    %al = (Vstar-(.6)^M)
    bet = al*ch3
    roots(al*[ch3 ch2 ch1 ch0] + [0 0 -1 Vb])
    r1=90
    %r2=75
    %r2=65
    r2=77
    r3=-35
    [bet*(r1+r2+r3)/al (1+bet*(r1*r2+r1*r3+r2*r3))/al (bet*r1*r2*r3 -Vb)/al]
    
    spaceconst = sqrt(epsil1*epsil2)
    if spaceconst < deltax
        disp('THE SPACE CONSTANT IS TOO SMALL COMPARED TO DELTAX')
    end
    
    figure
    Vtmp = -100:0.01:100;
    h1=ch0 + ch1*Vtmp + ch2*Vtmp.^2 + ch3*Vtmp.^3;
    f1=-Vtmp + Vb + (Vstar - (0/nB).^M).*h1;
    f2=-Vtmp + Vb + (Vstar - (0.5/nB).^M).*h1;
    f3=-Vtmp + Vb + (Vstar - (1/nB).^M).*h1;
    plot(Vtmp,f1,Vtmp,f2,'g:',Vtmp,f3,'r--')
    
    
    
end

%end

%---------------
function out1=karma_f(V,n,nB,M)

global Vstar Vb Vn Vh b ch0 ch1 ch2 ch3 capprox
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
global Vstar Vn Vh b ch0 ch1 ch2 ch3 capprox

%thet = V > Vn;
thet = 0.5*(1+tanh(capprox*(V-Vn)));

out2 = (1/b)*thet - n;

%function out=h2(V) % From Rappel's paper
%out = 72.83 - .83*V - .04976*V.^2 - (3.52e-4)*V.^3;