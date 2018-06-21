function out = karma_p2pzero_rev(x)
% Same as karma_p2p_rev, but subtract off initial condition (since we want to
% find an x at which karma_p2p_rev(x)-x = 0). 

xp1 = karma_p2p_rev2(x);
out=xp1-x;
