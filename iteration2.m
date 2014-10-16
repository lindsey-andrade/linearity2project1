function res = iteration2()

    %The main script for optimization project 1.
    %Model equation sourced from the article @ http://www2.hmc.edu/~dym/19-JSE2011Arches.pdf.
    
    q = 100000; 
    span = 50;


    function [f]=cost_func(vals)
        b=vals(1);
        h=vals(2);
        alpha=vals(3);
        
        R=abs(span/(2*sin(alpha)));
        arc_l=((2*alpha)/(2*pi))*pi*2*R; %Arc length
        %Calculating the cost of the beam based on a 1"x1"x36" beam costs $60
        base_beam_volume = 0.9144; % 36" = 0.9144m
        cost_per_unit_vol = 60/base_beam_volume; % $/m^3

        volume_of_beam = b * h * arc_l;

        f = volume_of_beam * cost_per_unit_vol;
    end

    function [c,ceq]=compute_stress(vals)
        %In the paper these are a sum, but I split them up because they are easier
        %to read that way :)
        b=vals(1);
        h=vals(2);
        alpha=vals(3);
        
        u_strength=670;
        I=b*h^3/12;%Second Area moment of the beam.
        A=b*h; %Cross sectional area
        R=abs(span/(2*sin(alpha)));
        ibar=I/(A*R^2); % Some weird way of consolidating geometric factors??
        lam=(alpha^4)/(12*ibar);
        %lam=(2.*R./h).*(1-cos(alpha)); %Another geometric adjustment.
        z= R-R*cos(alpha); %The segment height.
        arc_l=((2*alpha)/(2*pi))*pi*2*R; %Arc length
        str_cr_comp=((q.*R./A)./(1+(8/5).*lam.^2)).*((8/5)*lam.^2);
        str_cr_bend=((q.*R./A)./(1+(8/5).*lam.^2)).*((2.*z./h)*3.*lam);

        comp = str_cr_comp/1000000;
        bend = str_cr_bend/1000000;
        c = -comp-bend+u_strength;
        ceq=[];
    end
    x0=[1, 1, 10];
    [x, fval] = fmincon(@cost_func, x0, [], [], [], [], [0.1,0.1,0], [], @compute_stress)
end
