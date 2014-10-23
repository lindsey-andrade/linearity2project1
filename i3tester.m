cforce=80000;%N
mspan=90;%m
ustrength=670000000;%Pa
beamratio=4; %Maximum height:width ratio of solid beam.
bc=10; % $/rad
vc=4301; % $/m^3

%1 car = 1500kg and 5m long
%Span vs cost
%Max load vs cost


%[fval,x,exitflag]=iteration3(cforce,mspan,ustrength,beamratio,bc,vc);
%fval
%x
%exitflag

SPAN=linspace(50,100,30);
LOAD=linspace(10000,100000,30);
COST=zeros(length(SPAN),length(LOAD));
EXIT=zeros(length(SPAN),length(LOAD));
for i=1:length(SPAN)
    for j=1:length(LOAD)
        [fval,x,exitflag]=iteration3(LOAD(j),SPAN(i),ustrength,beamratio,bc,vc);
        if (exitflag>0)
            COST(i,j)=fval;
        else
            COST(i,j)=0;
        end
        disp(i)
        disp(j)
    end
end
figure(1);
surf(SPAN,LOAD,COST);
