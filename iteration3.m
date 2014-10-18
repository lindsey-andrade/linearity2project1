function iteration3()

P = 1000; 


    function [c, ceq] = stress(X)
        theta = stress(1); 
        a = stress(2); 
        b = stress(3); 
        c = stress(4); 
        
        z = sin(theta)*(a + (c-a)/2);
        
        A = (c-a)*b; 
        
        N1 = P/2; 
        N2 = P/2; 
        N = N1 + N2; 
        
        R = (c + a)/2;
        
        M1 = N1*z; 
        M2 = N2*z; 
        M = M1 + M2; 
        
        Am = b*log(c/a);
        
        tStress = N/A + M(A - a*Am)/(A*a(R*Am - A));
        
        cStress = N/A + M(A - c*Am)/(A*c(R*Am - A)); 
        
        totalStress = tStress + cStress;  
        
        % For optimization
        
        span = 2*z; 
        
        c = 670 - totalStress; 
        
        ceq = 100-span;
    end

end