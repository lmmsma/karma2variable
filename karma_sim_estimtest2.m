% This is a multicell estimator code, based off karma_sim_estimtest_1cell.m
% I'm not using karma_sim_estimtest, since I believe it's farther out of
% date. Warning: This still uses the smooth-theta version of the Karma
% model. The justification for this was to make the simulation consistent with the
% single-cell analytical Jacobians. It's not clear whether this is
% necessary or useful for multi-cell simulations, but if Jacobian
% estimation is implemented for larger systems, it may make sense to retain
% it.
%
% Changed usf from a scalar to a vector of size numpart (to facilitate
% applying I-feedback to each cell)
% Changed > fbon_time to >= fbon_time
% Compared to _p2p version on 9/3. p2p differs from this one in the
% following ways: 
% reads in IC instead of having a fixed one
% reads in parameters from a file
% different final time
% default numpart is one cell
% dataflag not set to 'zoh'
% Vinter is shortened to Vinter(1:numpart,:) 
% Vd is not initialized to zero and its first value is not set equal to Vinter(:,1)
% xout is initialized
% neither computes acoeffs 
% Verrornext is defined in terms of Vinter instead of Vd
% it doesn't show a 1000s counter
% save and various other commands are different based on whether L, bcl, numpart have been loaded
% xout is assigned values
% plotting and endmatter differ
% * has detloc and stimoutloc commented out
% * Vd and nd resizing due to 'dtscale' is commented out
% * has eliminated unused early state feedback setup
% * has no sf-only condition for trueindex == 1
% * deleted comments on feedback gains 
% 
% "*" means these changes were propagated to this version. 
% The parameters (with the exception of epsil1 due to single vs. multicell) were verified to be the same as of 9/3. 
% The main difficulty in reconciling the two versions is the handling of
% Vinter and Vd, although this could be done. 

% Older comments:
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
function [] = karma_sim_statefbtest(bcl,varargin)%karma_sim(karma_test,pubyear,stimperiod)
% simulate 2-variable Karma model for 1-D geometry
% use parameter values and functions from either PRL '93 paper or Chaos '94
% paper
% varargin{1} = integral-proportional gain vector (Verror variables only)
% varargin{2} = number of cells in fiber
% if karma_test is nonzero, attempt to reproduce Fig. 3 in '93 or '94 paper
% use values from paper corresponding to 'pubyear': 93 or 94
% if pubyear == 'ext', use 1D external potential version of equations
% outlined in Niels' biodomain model handout
% 6/28/07: changed detection scheme to try to make it work in the case when
% sensor and actuator are in the same cell.
% 7/02/07: added capability for multiple detectors

%clear all;
%close all;

global Vstar Vb Vn Vh b ch0 ch1 ch2 ch3 capprox

% load data_80cell_b225_7500_ol/34.mat V n usf
% Vi = V;
% ni = n;
% usfi = usf;

mkdir('data');
dfl = input('Delete old files? y or n: ','s')
if dfl=='y'
    delete data/*.mat;
end

fblocs = [16 46 76]; % cell indices used for local feedback
%fblocs = [16 31 46 61 76 91]; % cell indices used for local feedback
karma_test = 0;
%pubyear = varargin{1};
pubyear = 'ext';
dataflag = 'zoh'; % Use zero-order hold if set to 'zoh', otherwise use linearly interpolated data
%writeint = 2000; % maximum array length in ms
%writeint = round(2000/bcl)*bcl; % maximum array length in ms
%writeint = round(700/bcl)*bcl; % maximum array length in ms
%writeint = round(1000/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(3250/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(325/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(500/bcl)*bcl; % maximum array length in ms, 40 cells
writeint = round(250/bcl)*bcl; % maximum array length in ms, 80 cells
%writeint = round(7500/bcl)*bcl; % maximum array length in ms
%writeint = round(13000/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(1000/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(1600/bcl)*bcl; % maximum array length in ms, 20 cells
%writeint = round(1200/bcl)*bcl; % maximum array length in ms, 20 cells
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
    %    deltat = 0.025; %ms
    %    deltat = 0.001; %ms
    %    deltat = 0.002; %ms
    %    deltat = 0.004; %ms
    deltat = 0.008; %ms
    dtold=deltat;
    %    deltat = 0.05; %ms
    dtscale = 1;%4; % integer scaling factor
    %    deltat = 0.025*dtscale; %ms
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
    %        fiberlength = 12*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 37*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 43*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 211*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    fiberlength = 106*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 19*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 71*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %        fiberlength = 10*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %    fiberlength = 20*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %    fiberlength = 40*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %    fiberlength = 80*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
    %    fiberlength = 160*deltax;%deltax*8;%deltax*2;%deltax*40;%deltax*3; %0.5%5%2.62; %0.3%deltax*3;% cm
end
%finaltime = 2250; % ms
%finaltime = 7500;
%finaltime = 3400;
%finaltime = 1200;
%finaltime = 1145;
%finaltime = 300;
%finaltime = 600;
%finaltime = 915;
%finaltime = 4200;
%finaltime = 4160;
%finaltime = 5000;
%finaltime = 3*bcl; 
finaltime = 2000;
%finaltime = 3000;
%finaltime = 1000;
%finaltime = 13000;
%finaltime = 40000;
%finaltime = 3250;
%finaltime = 1500;
%finaltime = 4500;
%finaltime = 2*7500;
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
    %    tauN = 250;
    %    tauN = 200;
    tauN = 170;
%    tauN = 150;
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
        %        M=10; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        M=10; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=37.5; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=.08*10/.3; % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=1 % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        M=30 % lowering to 10 softens repolarization edges but doesn't amplify alternans
        %        nB=0.707; %0.680;% increasing makes more "sensitive"
        nB=1; % nB=0.707;
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
    %    epsil1 = max([0 Dg*De/(Dg+De)]); %.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.0053; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 9.3*0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.009; %0.007;%.009; % increasing too much means no AP formation; some runtime error if pushed below .007
    %    epsil1 = 0.001/4; % cm^2/ms
    %    epsil1 = 0.001/2; % cm^2/ms
    %    epsil1 = 0.001; % cm^2/ms
    %    epsil1 = 0.003; % cm^2/ms
    %    epsil1 = 0.002; % cm^2/ms
    %    epsil1 = 0.004; % cm^2/ms
    epsil1 = 0.005; % cm^2/ms
    %    epsil1 = 0.006; % cm^2/ms
    %    epsil1 = 0.007; % cm^2/ms
    %    epsil1 = 0.01; % cm^2/ms
    %    epsil1 = 0.015; % cm^2/ms
    %    epsil1 = 0.008; % cm^2/ms
    %    epsil1 = 0.03; % cm^2/ms
    %    epsil1 = 0.02; % cm^2/ms
    %    epsil1 = 0.017; % cm^2/ms
    %    epsil1 = 0.1; % cm^2/ms
    %    epsil1 = 0.05; % cm^2/ms
    %    epsil1 = 0.04; % cm^2/ms
    %    epsil1 = 0.04; % cm^2/ms
    %    epsil1 = 0.033; % cm^2/ms
    %    epsil1 = 0.2; % cm^2/ms % with epsil2=1, the propagation time to the distal end is shorter than the measured time (CV too large)
    %    epsil1 = 0.0225; % cm^2/ms % suitable for epsil2=0.2 and Vstar=1.88
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
    %    epsil2 = 5;
    %    epsil2 = 1;
    %    epsil2 = 1;
    %    epsil2 = 0.5;
    epsil2 = 0.7;
    %    epsil2 = 0.8;
    %    epsil2 = 0.3;
    %    epsil2 = 0.1;
    %    epsil2 = 0.15;
    %    epsil2 = 0.07;
    %    epsil2 = 0.08;
    %    epsil2 = 0.05;
    %    epsil2 = 1;
    %    epsil2 = 0.2;
    % Fourier number
    Fom = epsil1*deltat/(deltax)^2
end

% Vstar = 1.5415; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
% Vn = 1; % lowering seems to make alternans go away faster; making too high causes breakup/block
% Vh = 3; % lowering seems to make alternans go away faster, effect of increasing is not clear
%Vstar = 1.88; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.2; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.3; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.4; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.5; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 1.6; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 2; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 2.2; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 2.8; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 3.0; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 3.5; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%%Vstar = 3.6; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 3.8; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
Vstar = 4.0; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 4.5; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 5; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 2.5; % increasing increases max AP height, lowering too much can kill AP formation; doesn't seem to affect alternans
%Vstar = 6; 
%Vstar = 18; 
%Vn = -60; % lowering seems to make alternans go away faster; making too high causes breakup/block
Vn = -65; % lowering seems to make alternans go away faster; making too high causes breakup/block
Vh = -21; % lowering seems to make alternans go away faster, effect of increasing is not clear
%Vb = -85; % mV (offset from Rappel)
Vb = -90; % mV (offset from Rappel)
%Vb = -95; % mV (offset from Rappel)
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
%   ch0 = 72.83;
%   ch1 = -0.83;
%   ch2 = -0.04976;
%ch3 = -3.52e-4; % default for CinC '10
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
ch0 = 107.8831; % default for CinC '10
ch1 = -0.1319; % default for CinC '10
ch2 = -0.0465; % default for CinC '10
ch3 = -3.52e-4; % default for CinC '10
%ch0 = (107.8831/2.5)/5;
%ch1 = -0.1319*8*2.3;
%ch2 = -0.0465*1.87;
%ch3 = (-3.52e-4/1.08)*2;
%ch0 = (107.8831/2.5)/4;
%ch1 = -0.1319*8*1.35;
%ch2 = -0.0465/1.01;
%ch3 = (-3.52e-4/1.08)/1.01;

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
%fbdet = 1% index of detector used for feedback measurments
stimoutdet = 4; % index of detector used to measure downstream current extraction point

stimloc = deltax;%0.1 % location of stimulating electrode, measured from upstream end (x=0) (cm)
%stimheight = (deltat/dtold)*1.6; %1.5/5;%
%stimheight = (deltat/dtold)*3.2; %1.5/5;%
%stimheight = (deltat/dtold)*3.2*2.8*2; %1.5/5;%
stimheight = (deltat/dtold)*3.2*2.8; %1.5/5;%
%stimheight = 1.6; %(deltat/dtold)*0.2;%mV %Rappel et al., 10-20mV? Using Rappel's h(V), the "minimum" stimulus for
% a clean-looking AP seems to be either 0.2 height for 1 ms or 0.1 for 2
% ms, when deltat=0.05/8. Stimheight should go up by a factor of deltat/dtold if
% deltat is increased. For the old approximate h(V), I expect it to be higher, but I haven't
% checked yet. The units of iext as defined in the code should be uF/cm^2
%stimduration = 1/250; %deltat% assumptions: 1) period = 1 corresponds to 250ms, 2) typical stimulus duration is 1-5ms
%stimduration = 7;%2.5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
%stimduration = 15;%2.5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
%stimduration = 5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
%stimduration = 3; %
%stimduration = 2.5; %In new units, a stimulus duration of 20 timesteps (old duration) would be equivalent to 1 ms when deltat=0.05
stimduration = 1;
% 1 cell: A stimduration of 1ms will work at 1.6 height for 1 cell, but need 2.5ms duration to
% get the depolarization edge "all the way up" without changes in slope
% 2 cells: The upstroke in cell 1 takes a dip around -30mV using the 1cell
% settings
% 3 cells: The upstroke in cell 1 takes a dip around -45mV using the 1cell
% settings. If dur is increased to 4ms, the dip moves to -15mV.
% 10 or 20 cells: some dur >5 but <10 ms is probably ideal
% 20 cells: 7ms seems ok (1.6 height)
%stimoutloc = 0.2%2*deltax; %2.5;% location where current is extracted
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
%stimoutlocindex = round(stimoutloc/deltax); % cell index of location where current is extracted
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
    %    fbon_time = round(2800/deltat); %round(25/deltat);% time index at which control is switched on
    %    fbon_time = round(1450/deltat); %round(25/deltat);% time index at which control is switched on
    %    fbon_time = round(475/deltat); %round(25/deltat);% time index at which control is switched on
%    fbon_time = round(450/deltat); %round(25/deltat);% time index at which control is switched on
    fbon_time = round(50/deltat); %round(25/deltat);% time index at which control is switched on
    %    fbon_time = round((1350-11.25)/deltat); %round(25/deltat);% time index at which control is switched on
    %    fbon_time = round(15*bcl/deltat); % test an integer multiple of bcl for use with lyap_proof_check
    %    fbon_time = round(19950/deltat); % test an integer multiple of bcl for use with lyap_proof_check
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
%load xf_106cell_b160_4160_init85_0;
%V(:,1) = Xf(1:numpart); % initial condition
%  n(:,1) = Xf((numpart+1):end); % initial condition
%  V(1,1) = 2.5; % initial condition
%  n(1,1) = 0.5; % initial condition
%V(:,1) = -85; % initial condition
%n(:,1) = 0; % initial condition
%V(:,1) = -90.5; % initial condition
%n(:,1) = 0.41; % initial condition
V(:,1) = -90; % initial condition
n(:,1) = 0.7; % initial condition
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
    %    desfile = 'desval_80cell_b225_d001_finaltime5_cutout'; % attempt to use only fiber data for more accurate depolarization slopes
    %*    desfile = 'desval_1cell_b225_d001_finaltime40';
    %    desfile = 'desval_1cell_b325_smooththeta_ol';
    %    eval(['load ' desfile ' dapd']);
    %    desfile = 'desval_1cell_b225_smooththeta_sf5_cltraj'; % set writeint to something small for this one
    %    eval(['load ' desfile ' Vd nd ']);
    %    desfile = 'desval_1cell_b400_smooththeta_ol';
    %    eval(['load ' desfile ' Vd nd dapd']);
    %    desfile = 'markdata_proc_50_70'; % use karma_sim_estimtest_1cell(1002)
    %    desfile = '121702control1c_proc_50_70';
    %    desfile = '121702control1c_proc_50_70_80cells_4_11_3_12';
    %    desfile = '121702control1c_proc_50_55_80cells_temp';
    %    desfile = '121702control1c_proc_50_55_37cells_deltax0p05';
    %    desfile = '121702control1c_proc_50_55_37cells_deltax0p05_deltat0p025';
    %    desfile = '121702control1c_proc_50_55_71cells_temp';
    %    desfile = '121702control1c_proc_50_52_37cells_deltax0p05_deltat0p001';
    %    desfile = '121702control1c_proc_50_52_37cells_deltax0p05_deltat0p002';
    %    desfile = '121702control1c_proc_50_52_37cells_deltax0p05_deltat0p004';
    %    desfile = '121702control1c_proc_50_55_37cells_deltax0p05_deltat0p008';
    %    desfile = '121702control1c_proc_50_55_19cells_deltax0p1_deltat0p008';
    if dataflag ~= 'zoh'
        desfile = '121702control1c_proc_50_52_43cells_deltax0p05_deltat0p008';
        eval(['load ' desfile ' Vinter']);
        Vd = Vinter;
    else
        %    desfile = '121702control1c_proc_50_55_43cells_deltax0p05_deltat1';
        %    desfile = '121702control1c_proc_50_55_211cells_deltax0p01_deltat1';
%        desfile = '121702control1c_proc_50_55_106cells_deltax0p02_deltat1';
        desfile = '121702control1a_proc_1895_1900_106cells_deltax0p02_deltat1';
%        desfile = '121702control1c_proc_5rep_75_80_106cells_deltax0p02_deltat1';
%        eval(['load ' desfile ' Vinter dtact']);

        load(desfile,'Vinter', 'dtact')
        %%%%%%%%%%%%%%
        if 0
crange = [16 31 46 61 76 91];
minV=min(Vinter(crange,:),[],2);
maxV=max(Vinter(crange,:),[],2);
yoffset = 90; % Changed due to Robert's advice
yscale = 170./abs(maxV-minV);
Vmeasscale = diag(yscale)*(Vinter(crange,:)-repmat(minV,1,size(Vinter,2)))+repmat(minV,1,size(Vinter,2));
Vinter(crange,:) = Vmeasscale; 
        end
        %%%%%%%%%%%%%%
    end
    %%    Vd = Vinter(:,1)';
    %%    Vd = Vinter(1:12,:);
    %%    Vd = Vinter(1:12,:);
    Vd = zeros(size(iextd,1),size(iextd,2)+1);
    Vd(:,1) = Vinter(:,1);
    nd = zeros(size(Vd));
    if ~exist('dapd')
        dapd = -1;
    end
    %    eval(['load ' desfile ' Vd nd']);
    %     if exist('dtscale');
    %         Vd = Vd(:,1:dtscale:end);
    %         nd = nd(:,1:dtscale:end);
    %     end
end

% % Expand Vd sizes if necessary
% if size(Vd,1) < numpart
%     %Vdt = zeros(numpart,numstep);
%     %ndt = zeros(numpart,numstep);
%     %Vdt(1:size(Vd,1),:) = Vd;
%     %ndt(1:size(nd,1),:) = nd;
%     Vdt = zeros(numpart,writeintsteps+1);
%     ndt = zeros(numpart,writeintsteps+1);
%     Vdt(1:size(Vd,1),:) = Vd(:,1:writeintsteps+1);
%     ndt(1:size(nd,1),:) = nd(:,1:writeintsteps+1);
%     Vd = Vdt;
%     nd = ndt;
%     clear Vdt ndt
% end
%
% % comment out condition if you want Vd reading not to wrap around
% if writeintsteps <= size(Vd,2)%numstep <= size(Vd,2)
%     %    Vd = Vd(1:end,1:numstep);
%     %    nd = nd(1:end,1:numstep);
%     if ~exist('icfile')
%         %        Vd = Vd(1:numpart,1:numstep);
%         %        nd = nd(1:numpart,1:numstep);
%         % This should be the default case; des traj are shortened to
%         % writeinsteps+1 length
%         Vd = Vd(1:numpart,1:writeintsteps+1);
%         nd = nd(1:numpart,1:writeintsteps+1);
%     else
%         % 10 is chosen arbitrarily as a stimulus index that occurs after the
%         % default fbon_time. If Vd is not sufficiently long, this won't work
%         Vd = Vd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
%         nd = nd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
%     end
% else
%     %    Vd = Vd(1:end,:);
%     %    nd = nd(1:end,:);
%     if ~exist('icfile')
%         Vd = Vd(1:numpart,:);
%         nd = nd(1:numpart,:);
%     else
%         % 10 is chosen arbitrarily as a stimulus index that occurs after the
%         % default fbon_time. If Vd is not sufficiently long, this won't work
%         Vd = Vd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
%         nd = nd(1:numpart,stimindices(10):numstep+stimindices(10)-1);
%     end
% end
% % figure
% % hold on;
% % plot(Vd')
% % plot(nd','g')
% % plot(iextd','c')
% % figure
% % hold
% % %plot((1:length(Vinter))*deltat,Vinter(:,1),'b')
% % plot((1:size(Vinter,2))*deltat,Vinter,'b')
% % for ii = 1:numrep
% % %    plot([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)]*deltat, Vinter([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)],1)','g--')
% %     plot([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)]*deltat, Vinter(:,[((ii-1)*writeintsteps+1):(ii*writeintsteps+1)])','g--')
% % end
% if size(Vd,2) < writeintsteps %numstep
%     % Interpolate -desired quantities if necesary
%     % correct shift induced by initial times
%     Vdnew = interp1(dtold*(1:size(Vd,2))-dtold+deltat,Vd',deltat*(1:numstep),'linear',0);
%     ndnew = interp1(dtold*(1:size(nd,2))-dtold+deltat,nd',deltat*(1:numstep),'linear',0);
%
%     figure
%     plot(dtold*(1:size(Vd,2))-dtold+deltat, Vd')
%     hold
%     plot(deltat*(1:numstep),Vdnew,'g.')
%     plot(deltat*(1:numstep),iextd,'c')
%
%     figure
%     plot(dtold*(1:size(nd,2))-dtold+deltat, nd')
%     hold
%     plot(deltat*(1:numstep),ndnew,'g.')
%     plot(deltat*(1:numstep),iextd,'c')
%
% %    Vd = Vdnew';
% %    nd = ndnew';
% %    Vd = Vdnew; % may need to change back for newer versions of matlab
% %    nd = ndnew;
%     Vd = [Vdnew zeros(numpart,1)]; % this was added to fix an out-of-bounds error that occurs when writeindex == finaltime
%     nd = [ndnew zeros(numpart,1)];
% end

% state feedback term
%usf = zeros(1,writeintsteps);
usf = zeros(numpart,writeintsteps);
if exist('icfile')
    iext(:,1) = iextinitnew(:,2);
%    usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
    usf(stimlocindex,1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
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
%usfprev = 0;
usfprev = zeros(numpart,1);
Vnext = zeros(numpart,1);
nnext = zeros(numpart,1);
Verrornext = zeros(numpart,1);
nerrornext = zeros(numpart,1);
%usfnext= 0;
usfnext= zeros(numpart,1);
Vdnnext = 0; % only needed for periodic BC
Vupnext = 0;

perttime = 19900; %time of perturbation for Lyapunov test
%pertval = -10; % reset value in mV
pertval = 0; % reset value in mV

% ff=figure;
% hold on;

if dataflag == 'zoh'
    vdctr = 1; % Counter for ZOH on data
end

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
        %        Vd = Vinter([((ii-1)*writeintsteps+1):(ii*writeintsteps+1)],1)';
        %        % use with older files that have Vinter not scaled to fiber
        %        length
        %        Vd = Vinter(1:12,((ii-1)*writeintsteps+1):(ii*writeintsteps+1));
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
%        usf = zeros(1,writeintsteps);
        usf = zeros(numpart,writeintsteps);
        if exist('icfile')
            iext(:,1) = iextinitnew(:,2);
%            usf(1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
            usf(stimlocindex,1) = iextd(stimlocindex,1)-iext(stimlocindex,1);
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

        if trueindex >= fbon_time %***k > fbon_time %& abs(V(stimlocindex,k)) < 4.5
            %***                sext=[Verror(:,k); Verror(:,k)-Verror(:,k-1)];
            sext=[Verror(:,k); Verror(:,k)-Verrorprev];
            KIred1 = zeros(1,2*numpart);
            %                     KIred1(3:4)=place(Ared,B1(1:2),0.
            %                     999*eig(Ared));
            fac =0.999;
            %                fac =0.998; % this value doesn't seem to improve the bad transients with -0.01, -0.005
            %                fac =0.997; % this value doesn't seem to improve the bad transients with -0.01
            %                usf(k) = usf(k-1)-KIred1*sext;
            % ***               sext2=[Verror(:,k); Verror(:,k)-fac*Verror(:,k-1)];
            sext2=[Verror(:,k); Verror(:,k)-fac*Verrorprev];
            %                    sext2=[nerror(:,k); nerror(:,k)-fac*nerrorprev]; % use this to apply n-based feedback to the V-dynamics
            % ***               usf(k) = 1*(fac*usf(k-1)-KIred1*sext2);
%            usf(k) = 1*(fac*usfprev-KIred1*sext2);
% Comment the following line back in to apply I-feedback to the cell at stimlocindex
%            usf(stimlocindex,k) = 1*(fac*usfprev(stimlocindex)-KIred1*sext2);
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
%         if isnan(usf(k))
%             usf(k) = 0;
%         end
        if isnan(usf(:,k))
            usf(:,k) = zeros(numpart,1);
        end

        %         if abs(usf(k)) > stimheight/deltat/cc_ext
        %             usf(k) = stimheight/deltat/cc_ext*sign(usf(k));
        %         end
        %
        if trueindex >= fbon_time

            %            usf(k) = 25*Verror(1,k);
            %            usf(k) = 35*Verror(1,k); % For epsil1=0.1,epsil2=1, this doesn't give a significant improvement over 25.
            %            usf(k) = 10*Verror(1,k);
            %            iext(7,k) = -15*Verror(7,k);
            %            iext(7,k) = -20*Verror(7,k);
            %            iext(13,k) = -50*Verror(13,k);
            %            iext(7,k) = -50*Verror(7,k);
            %            iext(19,k) = -50*Verror(19,k);
            %             iext(7,k) = -25*Verror(7,k);
            %            iext(13,k) = -20*Verror(13,k);
            %             iext(13,k) = -25*Verror(13,k);
            %             iext(19,k) = -25*Verror(19,k);
            %            iext(7,k) = -50*Verror(7,k);
            %            iext(7,k) = -25*Verror(7,k) -25*Verror(19,k);
            %            iext(7,k) = -75*Verror(7,k);
            %            iext(13,k) = -75*Verror(13,k);
            %            iext(13,k) = -50*Verror(13,k);
            %            iext(19,k) = -50*Verror(19,k);
            %            iext(19,k) = -20*Verror(19,k);
            %            iext(25,k) = -50*Verror(25,k);
            %            iext(25,k) = -25*Verror(25,k);
            %            iext(25,k) = -20*Verror(25,k);
            %            iext(31,k) = -50*Verror(31,k);
            %            iext(31,k) = -25*Verror(31,k);
            %            iext(31,k) = -20*Verror(31,k);
            %            iext(37,k) = -50*Verror(37,k);
            %            iext(37,k) = -20*Verror(37,k);
            %            iext(37,k) = -75*Verror(37,k);
            %            iext(7,k) = -100*Verror(7,k);
            %            iext(13,k) = -100*Verror(13,k);
            % %           iext(19,k) = -100*Verror(19,k);
            %            iext(25,k) = -100*Verror(25,k);
            %            iext(37,k) = -100*Verror(37,k);
            %            iext(7,k) = -150*Verror(7,k);
            %            iext(13,k) = -150*Verror(13,k);
            %            iext(25,k) = -150*Verror(25,k);
            %            iext(37,k) = -150*Verror(37,k);
            %             iext(31,k) = -50*Verror(31,k);
            %             iext(61,k) = -50*Verror(61,k);
            %             iext(91,k) = -50*Verror(91,k);
            %             iext(121,k) = -50*Verror(121,k);
            %             iext(151,k) = -50*Verror(151,k);
            %             iext(31,k) = -20*Verror(31,k);
            %             iext(61,k) = -20*Verror(61,k);
            %             iext(91,k) = -20*Verror(91,k);
            %             iext(121,k) = -20*Verror(121,k);
            %             iext(151,k) = -20*Verror(151,k);
            %            iext(31,k) = -7*Verror(31,k) -7*Verror(61,k) -7*Verror(91,k);
            %            iext(91,k) = -7*Verror(61,k) -7*Verror(91,k) -7*Verror(121,k);
            %             iext(16,k) = -50*Verror(16,k);
            %             iext(31,k) = -50*Verror(31,k);
            %             iext(46,k) = -50*Verror(46,k);
            %             iext(61,k) = -50*Verror(61,k);
            %             iext(76,k) = -50*Verror(76,k);
            %              iext(16,k) = -120*Verror(16,k);
            %              iext(46,k) = -120*Verror(46,k);
            %              iext(76,k) = -120*Verror(76,k);
%                         iext(16,k) = -50*Verror(16,k);
%                         iext(46,k) = -50*Verror(46,k);
%                         iext(76,k) = -50*Verror(76,k);
%                         iext(16,k) = -25*Verror(16,k);
%                         iext(46,k) = -25*Verror(46,k);
%                         iext(76,k) = -25*Verror(76,k);
%                           iext(16,k) = -13*Verror(16,k);
%                           iext(46,k) = -13*Verror(46,k);
%                           iext(76,k) = -13*Verror(76,k);
%                          iext(16,k) = -70*Verror(16,k);
%                          iext(46,k) = -70*Verror(46,k);
%                          iext(76,k) = -70*Verror(76,k);
%                          iext(16,k) = -100*Verror(16,k);
%                           iext(46,k) = -100*Verror(46,k);
%                           iext(76,k) = -100*Verror(76,k);
%                           iext(16,k) = -150*Verror(16,k);
%                           iext(46,k) = -150*Verror(46,k);
%                           iext(76,k) = -150*Verror(76,k);
%                           iext(16,k) = -200*Verror(16,k);
%                           iext(46,k) = -200*Verror(46,k);
%                           iext(76,k) = -200*Verror(76,k);
%                         usf(16,k) = 75*Verror(16,k);
%                         usf(16,k) = 1*Verror(16,k);
%                         usf(16,k) = 10*Verror(16,k);
%                         usf(16,k) = 25*Verror(16,k);
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-25];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[-0.1;-25];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[-0.01;-25];
%            usf(fblocs(1),k) = fac*usfprev(fblocs(1))-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(fblocs(2),k) = fac*usfprev(fblocs(2))-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(fblocs(3),k) = fac*usfprev(fblocs(3))-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(31,k) = fac*usfprev(31)-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(61,k) = fac*usfprev(61)-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(91,k) = fac*usfprev(91)-[0 0 0 -1 -1 -1]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(1,k) = fac*usfprev(1)-[0.01 -0.01 -0.01 -2 -2 -2]*[Verror(fblocs,k); Verror(fblocs,k)-fac*Verrorprev(fblocs)];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-5];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-10];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-1];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-20];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-.4];
%            usf(fblocs,k) = fac*usfprev(fblocs)-[Verror(fblocs,k) Verror(fblocs,k)-fac*Verrorprev(fblocs)]*[0;-.3];
%            usf((fblocs(1)+1):(fblocs(2)-1),k) = mean(usf(fblocs(1:2),k));
%            usf((fblocs(2)+1):(fblocs(3)-1),k) = mean(usf(fblocs(2:3),k));
%            usf(1:(fblocs(1)-1),k) = usf(fblocs(1),k);
%            usf((fblocs(1)+1):(fblocs(2)-1),k) = (1-(1:length((fblocs(1)+1):(fblocs(2)-1)))/length((fblocs(1)+1):(fblocs(2)-1)))*usf(fblocs(1),k)+((1:length((fblocs(1)+1):(fblocs(2)-1)))/length((fblocs(1)+1):(fblocs(2)-1)))*usf(fblocs(2),k);
%            usf((fblocs(2)+1):(fblocs(3)-1),k) = (1-(1:length((fblocs(2)+1):(fblocs(3)-1)))/length((fblocs(2)+1):(fblocs(3)-1)))*usf(fblocs(2),k)+((1:length((fblocs(2)+1):(fblocs(3)-1)))/length((fblocs(2)+1):(fblocs(3)-1)))*usf(fblocs(3),k);
%            usf((fblocs(3)+1):end,k) = usf(fblocs(3),k);
%usf(fblocs(1),k) = -.5*Verror(fblocs(1),k);             
%usf(fblocs(1),k) = -1*Verror(fblocs(1),k);             
%usf(fblocs(1),k) = -5*Verror(fblocs(1),k);             
%usf(fblocs,k) = -1*Verror(fblocs,k);             
%                         usf(16,k) = 2*Verror(16,k);

%                         usf(46,k) = 2*Verror(46,k);
%                         usf(76,k) = 2*Verror(76,k);
%            usf(16,k) = fac*usfprev(16)-[-0.001 -2]*[Verror(76,k); Verror(16,k)-fac*Verrorprev(16)];
            Vem = Verror(:,k);
            Vem(Vem<-5) = -40; 
            Vem(Vem>5) = 40; 
%                        usf(16,k) = 2*Vem(16);
%                        usf(46,k) = 2*Vem(46);
%                        usf(76,k) = 2*Vem(76);
            
            
        end
%            iext(fblocs,k) = -usf(fblocs,k); 
            iext(:,k) = -usf(:,k); 
%        iext([stimlocindex],k) = iextd(stimlocindex,k)-usf(k);
        iext([stimlocindex],k) = iextd(stimlocindex,k)-usf(stimlocindex,k);

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
    if trueindex >= fbon_time % This set should be consistent with b2 = -0.05
        %                un(:,k) = 0.005*Verror(:,k)/tauN; % proportional feedback term applied to n instead of V
        %                gains in foregoing line (using Verror for fb): negative
        %                gains, such as -0.005/tauN, don't appear to work.
        %                However, 0.005*deltat/tauN stabilizes, 0.003 stabilizes slowly, and 0.001 doesn't.
        %                un(:,k) = 0.05*Verror(:,k)/tauN; % very aggressive; no apparent oscillations
        %                un(:,k) = 0.01*Verror(:,k)/tauN; % several apparent oscillations
        %                un(:,k) = 0.03*Verror(:,k)/tauN; % almost 1 oscillation (debatable)
        %                   un(1:6,k) = repmat(0.05*Verror(1,k)/tauN,6,1); % putting the feedback on 3 cell ranges 1:6, 7:12, 13:18 isn't helpful; corrections don't propagate significantly, and the feedback causes secondary AP-like formations in affected cells. Suggest reducing gain if many cells are receiving feedback.
        %                   un(7:12,k) = repmat(0.05*Verror(7,k)/tauN,6,1); %
        %                   un(13:18,k) = repmat(0.05*Verror(13,k)/tauN,6,1); %
        %                  un(19,k) = repmat(0.05*Verror(19,k)/tauN,1,1); % almost 1 oscillation (debatable)
        %                  un(25,k) = repmat(0.05*Verror(25,k)/tauN,1,1); % almost 1 oscillation (debatable)
        %                  un(31,k) = repmat(0.05*Verror(31,k)/tauN,1,1); % almost 1 oscillation (debatable)
        %               un(:,k) = 0.04*Verror(:,k)/tauN; % almost 1 oscillation (debatable)

        %                un(:,k) = -0.75*nerror(:,k)/tauN; % proportional feedback term applied to n instead of V
        %                gains in foregoing line: -1/tauN and -0.75/tauN stabilize, -0.5/tauN stabilizes very
        %                slowly, -0.25/tauN does not stabilize
%                           un(16,k) = .05*Verror(16,k);
%                           un(46,k) = .05*Verror(46,k);
%                           un(76,k) = .05*Verror(76,k);
%                           un(16,k) = .0008*Verror(16,k);
%                           un(46,k) = .0008*Verror(46,k);
%                           un(76,k) = .0008*Verror(76,k);
%                             un(16,k) = .0001*Verror(16,k);
%                             un(46,k) = .0001*Verror(46,k);
%                             un(76,k) = .0001*Verror(76,k);
%                             un(fblocs,k) = .0001*Verror(fblocs,k);
%                             un(fblocs,k) = .0002*Verror(fblocs,k);
%                             un(fblocs,k) = .00001*Verror(fblocs,k);
%                             un(fblocs,k) = .00005*Verror(fblocs,k);
%                              un(fblocs,k) = .00002*Verror(fblocs,k);
%             un(1:(fblocs(1)-1),k) = .85*un(fblocs(1),k);
%             un((fblocs(1)+1):(fblocs(2)-1),k) = (1-(1:length((fblocs(1)+1):(fblocs(2)-1)))/length((fblocs(1)+1):(fblocs(2)-1)))*un(fblocs(1),k)+((1:length((fblocs(1)+1):(fblocs(2)-1)))/length((fblocs(1)+1):(fblocs(2)-1)))*un(fblocs(2),k);
%             un((fblocs(2)+1):(fblocs(3)-1),k) = 1.25*((1-(1:length((fblocs(2)+1):(fblocs(3)-1)))/length((fblocs(2)+1):(fblocs(3)-1)))*un(fblocs(2),k)+((1:length((fblocs(2)+1):(fblocs(3)-1)))/length((fblocs(2)+1):(fblocs(3)-1)))*un(fblocs(3),k));
%             un((fblocs(3)+1):end,k) = 1.5*un(fblocs(3),k);
%un(1:numpart,k) = 0.00005*Verror(16,k) + 0.00005*Verror(46,k); 
un(1:numpart,k) = 0.000005*Vem(16) + 0.000005*Vem(46)+ 0.00001*Vem(76); 

    end
    %        nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k); % This one is consistent with b2=1 assumption
    nnext = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN - deltat*un(:,k); % This should be consistent with b2 = -0.05
    %        n(:,k+1) = n(:,k) + deltat*karma_g(V(:,k),n(:,k))/tauN + un(:,k);
    % Perturb V for Lyapunov test
    %%       if trueindex==round(19900/deltat)
    %        if trueindex==round(perttime/deltat)
    %            Vnext = pertval;%-80;
    %        end

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
    if dataflag == 'zoh' & k > 1 & ~mod(k,1000*dtact/deltat)
        vdctr=vdctr+1;
    end
    if dataflag ~= 'zoh'
        Verrornext = Vd(:,k+1)-Vnext;
    else
        Verrornext = Vinter(:,vdctr)-Vnext;
    end
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
            %            falling edge, not when any detector detects a falling edge)
        end % if falling edge
    end % if deton
    if ~mod(trueindex,1000) % *** ~mod(k,1000)
        disp(trueindex) % *** disp(k)
    end
    Vprev = V(:,k);
    nprev = n(:,k);
    Verrorprev = Verror(:,k);
    nerrorprev = nerror(:,k);
%    usfprev = usf(k);
    usfprev = usf(:,k);
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
%    eval(['save data/' num2str(ii) ' un hprime acoeff* usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
%eval(['save data/' num2str(ii) ' un usf Verror nerror iext V n apds apendindices dctr rctr *prev *next']);  % s iint Vx_dn Vx_up Vdn Vup
    datafname = ['data/' num2str(ii)];
    save(datafname,'un', 'usf', 'Verror', 'nerror', 'iext', 'V', 'n', 'apds', 'apendindices', 'dctr', 'rctr', '*prev', '*next');  
clear hprime acoeff* usf Verror nerror iext V n un
if ii == numrep
    %        save data/configinfo B1 B2 stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength
    %        save data/configinfo pert* B1 B2 stimperiod nB M epsil1 epsil2 tauN Vstar Vn Vh Vb b ch0 ch1 ch2 ch3 cc_* finaltime deltat deltax numstep numpart stimduration stimheight writeint numrep detvec fbon_time KIred1 stimindices iextd writeintsteps fiberlength stimlocindex dapd desfile
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
if 1
    %        crange = [1 7 13 19 25 31]; %detlocindices
    %        crange = [1 4 7 10 13 16]; %detlocindices
    %        crange = [1 7 13 19 25 31 37]; %detlocindices
    %        crange = [37]; %detlocindices
    %        crange = [7 13 19 25 31 37]; %detlocindices
    %        crange = [31 61 91 121 151 181]; %detlocindices
    crange = [16 31 46 61 76 91]; %detlocindices
    h=figure;
    clear p p1;
    for ii=1:numrep % number of data writing cycles
        subplot(2,1,1)
        hold on;
        eval(['load data/' num2str(ii)])
        %       p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps)*(-9/.0051),'g--');
        %       p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps)*(-9/.003),'g--');
        %            p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([stimlocindex],1:writeintsteps),'--');
        %            p1(ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un(1,1:writeintsteps),'--');
        p1(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,iext([crange 1],1:writeintsteps));
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
else
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
    legend('i','i^d')
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

if 1    %    detvec = (71:80)';
    colorcounter = 1;
    markerindex = 1;
    hwhat=figure
    subplot(2,1,1)
%    subplot(2,1,2)
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
end

%        crange = [1 12 23 35 46 58 2]; %detlocindices
figure(hwhat)
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


% Define an offset colororder
colororderoffset = colororder([3:7 1:2],:);

%    numrep = 20;
h=figure;
%    crange = 1:10; %detlocindices
%    crange = [1 12]; %detlocindices
%    crange = [1 7 13 19 25 31]; %detlocindices
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
    %    p2=plot(deltat*(1:1:numstep),Vd(detlocindices,1:1:numstep),'--');
    %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Vd(detlocindices,1:writeintsteps),'--');
    %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vd(crange,1:50:writeintsteps),'--');
    %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vinter([((ii-1)*writeintsteps+1):50:ii*writeintsteps],crange),'--');
    %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):50:ii*writeintsteps])*deltat,Vinter(crange,[((ii-1)*writeintsteps+1):50:ii*writeintsteps]),'--');
    %        p2(:,ii)=plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,Verror(crange,1:writeintsteps)+V(crange,1:writeintsteps),'--');
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
%ti=title(['WARNING: You aren''t using the whole data trajectory, just a repeated portion']);
%    legend([p1(end); p2(end)],'V_i','V^d_i')
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
yl=ylabel({'ion channel', 'variable'});
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

eval(['load data/' num2str(max([numrep-1 1]))])
figure
hold
ii=numrep
plot(([((ii-1)*writeintsteps+1):ii*writeintsteps-1])*deltat,diff(V(37,1:writeintsteps))/deltat)
plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,V(37,1:writeintsteps),'g:')

h3=figure;
hold on;
%            eval(['load data_80cell_b225_15000_ol/configinfo'])
%            writeintsteps = round(writeint/deltat);
%            fiberlength = deltax*numpart;
reprange = (numrep-5):(numrep-1);%10:20%1:10%11:20%41:50%
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
caxis([-100 50])
colorbar
set(gca,'FontSize',fs);

h8=figure;
hold on;
for ii=1:numrep % number of data writing cycles
    eval(['load data/' num2str(ii)])
    %    plot((1:numstep)*deltat,un([stimlocindex],1:numstep))
    plot(([((ii-1)*writeintsteps+1):ii*writeintsteps])*deltat,un([crange],1:writeintsteps))
    xlabel('time')
    ylabel('npert')
    title(['Perturbations applied to n in cells ' num2str(crange)])
    clear hprime acoeff* usf Verror nerror iext V n un
end

%end


% mean(rankctr)
% max(rankctr)
% min(rankctr)
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

% pause
% fn= input('Enter a filename: ','s')
% eval(['save ' fn]);

% Record an OL trajectory for use as desired trajectory to try to fulfill assumptions of Lyapunov proof
% load data/1
% Vd = V(1,(stimindices(36)):stimindices(40)-1);
% nd = n(1,(stimindices(36)):stimindices(40)-1);
% iextd=iextd(:,1:length(Vd));
% figure
% plot(Vd')
% hold
% plot(nd')
% plot(iextd')
% save desval_1cell_b325_smooththeta_ol Vd nd iextd
%
% eval(['load data/' num2str(numrep)])
% Vd = V;
% nd = n;
% save desval_1cell_b400_smooththeta_ol Vd nd iextd
% figure
% plot(iextd')
% hold
% plot(Vd')
% plot(nd')
% Vdif = V(1,1:8000)-V(1,8001:16000);
% min(Vdif)
% max(Vdif)
%
% eval(['load data/' num2str(numrep-1)])
% Vd = V;
% nd = n;
% save desval_1cell_b400_smooththeta_ol_pert Vd nd iextd
% figure
% plot(iextd')
% hold
% plot(Vd')
% plot(nd')
% Vdif = V(1,1:8000)-V(1,8001:16000);
% min(Vdif)
% max(Vdif)
%
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