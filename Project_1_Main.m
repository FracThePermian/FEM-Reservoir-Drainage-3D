close all
clear all
clc
[para] = reservoir; %class to help clean up program
 x = linspace(para.dx,para.L-para.dx/2,para.NX); y = linspace(para.dy,para.W-para.dy/2,para.NY); %for graphing
P = 1000*ones(para.N,1);    BC = zeros(para.N,1);      P_B = zeros(para.N,1);
well_location = [10000,10000]; q_well = [-1000]; %preppring for later use and when there is multiple wells.

para_wells = 5050; %well_location(1)/para.dx+well_location(2)/para.dy*para.NX; %which is 5050
 
[ T, Q,B] =  TBQ_box_f(BC,P_B,para_wells,q_well); %a box that spits out matrices T, B, Q
 
dt = 1; t = 0; t_end = 51; n = 1;
P_days = zeros(t_end/dt,para.N);
while t < t_end
    P_2 = P; P = (T+B/dt)\(B*P_2/dt+Q);
    P_days(n,:) = P;
    t = t + dt; n = n + 1;
end

[I,J] = meshgrid(x,y); %Makes the meshgrid
P_plot50 = reshape(P_days(50,:),para.NX,para.NY);
surf(I,J,P_plot50) %This is the 3-D plot of well pressure wrt res. dimensions
colorbar
xlabel('Reservoir Length (ft)','FontSize',14)
ylabel('Reservoir Width (ft)','FontSize',14)
zlabel('Reservoir Pressure(psi)','FontSize',14)
title('3-D Pressure Visualization','FontSize',20)

figure %This prints our the contour map of pressures over the well 
contourf(I,J,P_plot50)
colorbar
xlabel('Reservoir Length(ft)','FontSize',14)
ylabel('Reservoir Width(ft)','FontSize',14)
title('Reservoir Contour','FontSize',20)

figure %This plots the final analytical and numerical reservoir pressure wrt length 
P_d1 = P_days(1,:);  P_d10 = P_days(10,:);  P_d50 = P_days(50,:);
P1=P_d1(5001:5100);  P10=P_d10(5001:5100);  P50=P_d50(5001:5100);
plot(x,P1,'r+',x,P_analytical(x-10000,1),'r-',x,P10,'b+',x,P_analytical(x-10000,10),'b-',x,P50,'k+',x,P_analytical(x-10000,50),'k-');
xlabel('Reservoir Length(ft)','FontSize',14)
ylabel('Reservoir Pressure(psi)','FontSize',14)
title('Reservoir Pressure vs Reservoir Length (Numerical vs. Analytical)','FontSize',20)
legend('t = 1 day','Analytical 1 day','t = 10 days','Analytical 10 days','t = 50days','Analytical 50 days')

