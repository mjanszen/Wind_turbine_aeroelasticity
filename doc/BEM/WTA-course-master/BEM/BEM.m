
function [Rx,FN,FT,P]=BEM(v0,omega,pitch)
%------------------------------------------------
% Blade Element Momentum
%------------------------------------------------
% this code is for the course: wind turbine aeroelasticity
% Delft university of technology

%-------------------STATEMENT-------------------%
% a: axial induction
% a_prime: tangential induction
% Phi: inflow angle ¦Õ
% Alpha: local attack angle ¦Á
% Theta: twist angle
% pitch: blade pitch angle
% Sigma: solidity
% Cl: lift coefficient
% Cd: drag coefficient
% Cn: components of n (along wind direction) in Cartesian coordinate system
% Ct: components of t (perpendicular to wind direction) in Cartesian coordinate system
% v0: inflow wind speed
% omega: blade rotation angular speed
% r: radius of blade
%-----------------------------------------------%

%---------------START SIMULATION----------------%

%------------------------------------------------
% fixed parameters
%------------------------------------------------
B=3;           %number of blades
R=63;          %rotor radius
hubrad=1.5;    %hub radius
rou=1.225;     %density of air
EPS=0.00001;    %iterative precision tolerance

%------------------------------------------------
% Initialization & Iteration
%------------------------------------------------
%initialization: initial value of inductions
a=0;a_prime=0;

%import Blade section file
BS=importdata('Blade\Blade section\Blade section.dat').data;

%import Aero data files
Readfiles = dir(fullfile('Blade\Aero data\','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade\Aero data\',Readfiles(i).name));
end

NBS=length(BS);    %Number of blade sections
% define vectors for blade section locations and loads in two directions
Rx=zeros(NBS,1);FN=zeros(NBS,1);FT=zeros(NBS,1);

% LOOP: from the root section to the tip section
for i=1:NBS
    ADofBS=BS(i,2); % read airfoil number for each section
    r=BS(i,3);      % read radius
    Rx(i)=r;        % record radius
    dr=BS(i,4);     % read segment length
    Theta=BS(i,5);  % read twist angle
    chord=BS(i,6);   % chord length
    alpha=AD{ADofBS}(:,1);% coefficients table: AOA
    Cl=AD{ADofBS}(:,2); % coefficients table: lift coe
    Cd=AD{ADofBS}(:,3); % coefficients table: drag coe
    Sigma=chord*B/(2*pi*r); % solidity
    ax=a;                     %change value
    ax_prime=a_prime;         %change value
    a=ax-10*EPS;              %generate error, active iteration
    a_prime=ax_prime-10*EPS;  %generate error, active iteration
    
    numite=0; % iteration counter
    %iteration, stop when error is smaller than EPS
    while abs(ax-a)>=EPS || abs(ax_prime-a_prime)>=EPS
        numite=numite+1;
        
        % record results of last step
        a=ax;
        a_prime=ax_prime;
        
        % inflow angle
        Phi=atan((1-a)*v0/((1+a_prime)*r*omega));
        Phi=rad2deg(Phi);
        
        %AOA
        Alpha=Phi-Theta-pitch;
        
        % find Cl and Cd
        Cla=interp1(alpha,Cl,Alpha);
        Cda=interp1(alpha,Cd,Alpha);
        
        %projection in and out of plane
        Cn=Cla*cosd(Phi)+Cda*sind(Phi);
        Ct=Cla*sind(Phi)-Cda*cosd(Phi);
        
        %Prandtl Loss
        f_tiploss = B/2*(R-r)/(r*sind(Phi));
        F_tiploss = (2/pi)*acos(exp(-f_tiploss));
        f_hubloss = B/2*(r-hubrad)/(r*sind(Phi));
        F_hubloss = (2/pi)*acos(exp(-f_hubloss));
        F = F_tiploss*F_hubloss;
        
        %Glauert Correction
        ac=0.2;
        if ax>ac
            K=4*F*sind(Phi)^2/(Sigma*Cn);
            ax=0.5*(2+K*(1-2*ac)-sqrt((K*(1-2*ac)+2)^2+4*(K*ac^2-1)));
        else
            ax=1/(4*F*(sind(Phi))^2/(Sigma*Cn)+1);
        end
        ax_prime=1/(4*F*sind(Phi)*cosd(Phi)/(Sigma*Ct)-1);
        
        % in case of iterative convergence failure
        if numite>=100
            ax=0.3;
            ax_prime=0.1;
        end
    end
    
%------------------------------------------------
% Result
%------------------------------------------------
    % update value
    a=ax;
    a_prime=ax_prime;
    
    % force in two directions
    FN(i)=0.5*rou*((r*omega*(1+a_prime))^2+(v0*(1-a))^2)*chord*Cn*dr;
    FT(i)=0.5*rou*((r*omega*(1+a_prime))^2+(v0*(1-a))^2)*chord*Ct*dr;
    % bending moment
    Mx(i)=FT(i)*r;
end

M=sum(Mx); % rotor torque from one blade
P=M*omega*3*0.944;  % Power
end