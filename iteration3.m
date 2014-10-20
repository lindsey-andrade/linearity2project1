function iteration3()

P = 1000; 

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
        
        b_c=1*b*((a+c)/2)^3; %The cost per bend degree
        v_c=4301; % $/m^3 from http://www.roymech.co.uk/Useful_Tables/Matter/Costs.html
        
        vol=(theta/2)*(a+c)*(b*(c-a));
        
        %The statement on the next line says...
        %The beam costs v_c per unit of volume
        %And additional costs are incurred for pre-bending.
        t_c=b_c*theta*+(v_c*vol); 
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
        
        tStress = N/A + M*(A - a*Am)/(A*a*(R*Am - A));
        
        cStress = N/A + M*(A - c*Am)/(A*c*(R*Am - A)); 
        
        totalStress = tStress + cStress;  
        
        % For optimization
        
        span = 2*e_x;
        
        c = totalStress-670000;
        
        ceq = 10-span;
    end

% data=zeros(1,90);
% test_angle=0;
% for j=1:90
%     data(j)=stress([test_angle/57.3, 30, 50, 80]);
%     test_angle=test_angle+1;
% end
% plot(data)

%test_data=[0 4 3 10]
%s=stress(test_data)
%v=volume(test_data)

[x,fval]=fmincon(@cost,[pi/6;1;1;2],[0 1 0 -1],[0],[],[],[0;0;0;0],[pi/2;1e+10;1e+10;1e+10],@stress);
disp(x)
disp(fval)
end