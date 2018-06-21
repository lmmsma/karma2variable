% Try to use Newton-Krylov solver to find stabilizing gains

clear variables
%initialvalue = -0.01; 
initialvalue = zeros(212,1);
initialvalue(16) = -0.01; 
%gaintol = 0.13; 
gaintol = 0.05; 

tstartL = tic
[Lnsoli, lit_hist, lierr, l_hist] = nsoli(initialvalue,'karmacleig_zero',[gaintol gaintol])
tendL=toc(tstartL)

%save L_106cell_b230_des0p85_pwrit2 Lnsoli lit_hist lierr l_hist initialvalue tendL gaintol 
%save L_106cell_b230_des0p85_tol0p05_pwrit2 Lnsoli lit_hist lierr l_hist initialvalue tendL gaintol 

clear variables
initialvalue = -10; 
gaintol = .5; 

tstartL = tic
[Lnsoli, lit_hist, lierr, l_hist] = nsoli(initialvalue,'karma_altamp',[gaintol gaintol])
tendL=toc(tstartL)


clear variables
initialvalue = 5; 
%opt_redtol=optimset('TolFun',1e-10); %default is 1e-6
tstartL = tic
%[Llsqn,resnorm,residual,exitflag,output] = lsqnonlin('karma_altamp',initialvalue,-50,-1,opt_redtol)
[Llsqn,resnorm,residual,exitflag,output] = lsqnonlin('karma_altamp',initialvalue,-50,5)
tendL=toc(tstartL)


clear variables
initialvalue = -0.01; 
tstartL = tic
[Llsqn,resnorm,residual,exitflag,output] = lsqnonlin('karmacleig_zero',initialvalue,-50,-0.001)
tendL=toc(tstartL)


clear variables
initialvalue = [-0.01]; 
tstartL = tic
[Lfmc,fval,exitflag,output] = fmincon('karmacleig_zero',initialvalue,[],[],[],[],[-10],[-0.001])
tendL=toc(tstartL)


clear variables
opt_intpt =optimset('Algorithm','interior-point');
initialvalue = [-0.01]; 
tstartL = tic
[Lfmc,fval,exitflag,output] = fmincon('karmacleig_zero',initialvalue,[],[],[],[],[-10],[-0.001],[],opt_intpt)
tendL=toc(tstartL)


clear variables
initialvalue = [-0.01;0]; 
tstartL = tic
[Lfmc,fval,exitflag,output] = fmincon('karmacleig_zero',initialvalue,[],[],[],[],[-5;-0.001],[-0.001;0.001])
tendL=toc(tstartL)

