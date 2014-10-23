function [fval,x,exitflag]=iteration3(cforce,mspan,ustrength,beamratio,bc,vc)

P = cforce; %Linear Force applied @ the crown
min_span=mspan; %minimum distance the bridge must span.
ultimate_strength=ustrength; %maximum stress the material can take.

    function [vol] = volume(X)
        theta = X(1); 
        a = X(2); 
        b = X(3); 
        c = X(4);
        %Some simplifications are used here.
        vol=(theta/2)*(a+c)*(b*(c-a));
    end

    function [t_c] = cost(X)
        theta = X(1); 
        a = X(2); 
        b = X(3); 
        c = X(4);
        
        b_c=bc; %The cost per bend degree
        v_c=vc; % $/m^3 from http://www.roymech.co.uk/Useful_Tables/Matter/Costs.html
        
        vol=(theta/2)*(a+c)*(b*(c-a));
        
        %The statement on the next line says...
        %The beam costs v_c per unit of volume
        %And additional costs are incurred for pre-bending.
        t_c=b_c*theta+(v_c*vol); 
    end

    function [c, ceq] = stress(X)
        theta = X(1); 
        a = X(2); 
        b = X(3); 
        c = X(4); 
        
        e_x = sin(theta)*((a + c)/2); %X coordinate of the bridge endpoint
        e_y = cos(theta)*((a+c)/2); %Y coordinate of the bridge endpoint
        
        A = (c-a)*b; %Cross section Area
        
        N1 = P/2; 
        N2 = P/2; 
        N = N1 + N2; 
        
        R = (c + a)/2;
        
        M1 = N1*e_x; 
        M2 = N2*e_x; 
        M = M1 + M2; 
        
        Am = b*log(c/a);
        
        tStress = N/A + M*(A - a*Am)/(A*a*(R*Am - A)); %will be +
        
        cStress = N/A + M*(A - c*Am)/(A*c*(R*Am - A)); %will be -
        
        % For optimization
        
        span = 2*e_x;
        
        c = tStress-ultimate_strength;
        
        ceq = min_span-span;
    end

%Start with some reasonable init values.
i_theta=pi/8;
i_a=(min_span/sin(i_theta))/2;
i_b=1;
i_c=((min_span/sin(i_theta))/2)+1;

op = optimset('fmincon');
op = optimset(op,'MaxFunEvals',50000,'MaxIter',10000);
[x,fval,exitflag]=fmincon(@cost,[i_theta,i_a,i_b,i_c],[0 1 0 -1;0 -1 -beamratio 1;],[0;0],[],[],[0;0;0;0],[pi/2;1e+10;1e+10;1e+10],@stress,op);
end