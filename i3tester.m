cforce=80000;%N
mspan=90;%m
ustrength=17000000;%Pa
beamratio=4; %Maximum height:width ratio of solid beam.
bc=100; % $/rad
vc=4301; % $/m^3

[fval,x,exitflag]=iteration3(cforce,mspan,ustrength,beamratio,bc,vc);
fval
x
exitflag