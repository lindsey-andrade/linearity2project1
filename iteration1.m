function res = iteration1(b, h)

    %The main script for optimization project 1.
    %Model equation sourced from the article @ http://www2.hmc.edu/~dym/19-JSE2011Arches.pdf.


    %Constraint constants:
    %For consistency, let's say this is in MGS. Anything in here is set and
    %forget.
    q=100; %The force applied radially at the crown.
    span=10; %The width of the bridge
    %Decision Variables:
    % b=.1; %The base of the rectangular cross section
    % h=.1; %The height of the rectangular cross section
    alpha=linspace(0.01,pi/2,1000)'; %Half the arc angle. Radians plz.

    %Expression Definitions:
    I=b.*h.^3./12;%Second Area moment of the beam.
    A=b.*h; %Cross sectional area
    R=span./(2*sin(alpha));
    ibar=I./(A.*R.^2); % Some weird way of consolidating geometric factors??
    lam=(alpha.^4)./(12.*ibar);
    %lam=(2.*R./h).*(1-cos(alpha)); %Another geometric adjustment.
    z= R-R.*cos(alpha); %The segment height.
    arc_l=((2*alpha)/(2*pi)).*pi*2.*R;

    %In the paper these are a sum, but I split them up because they are easier
    %to read that way :)
    str_cr_comp=((q.*R./A)./(1+(8/5).*lam.^2)).*((8/5)*lam.^2);
    str_cr_bend=((q.*R./A)./(1+(8/5).*lam.^2)).*((2.*z./h)*3.*lam);

    comp = str_cr_comp/1000000;
    bend = str_cr_bend/1000000;
    sum = comp + bend;

    res = [comp, bend, sum];

    figure(1);
    clf();
    hold on;
    plot(alpha, str_cr_comp/1000000,'b');
    plot(alpha, str_cr_bend/1000000,'r');
    plot(alpha, (str_cr_bend+str_cr_comp)/1000000,'g');
    xlabel('Alpha (rad)');
    ylabel('Stress (MPa)');
    legend('Compressive Stress','Bending Stress','Total Stress');
    title('Stress as a function of Alpha');

    figure(2);
    plot(alpha,arc_l);
    title('Arc Length as a function of Alpha');
    xlabel('alpha(rad)');
    ylabel('length (m)');

    %The yield strength of steel of 434MPa! That's A LOT.
end
