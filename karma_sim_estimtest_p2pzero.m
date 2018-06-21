function out = karma_sim_estimtest_p2pzero(x)
% Same as karma_sim_estimtest_p2p, but subtract off initial condition (since we want to
% find an x at which karma_sim_estimtest_p2p(x)-x = 0). 

xp1 = karma_sim_estimtest_p2p(x);
out=xp1-x;
