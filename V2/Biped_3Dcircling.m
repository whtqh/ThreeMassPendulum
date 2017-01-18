clear;
clc;
close all;

%constant
g = 9.8;     %gravity
time = 0.01;  %the simulation time

a = 10; %positive weight of the evaluation function 4.58
b = 1;

%input
zc = 0.8;       %height
Tc = sqrt(zc/g);

Tsup = 0.8;     %support period
N = round(Tsup/time);
C = cosh(Tsup/Tc);
S = sinh(Tsup/Tc);

%walking parameters
sx = [0  ,0.25,0.25,0.25,0]; %step length
sy = [0.2,0.2 ,0.2 ,0.2 ,0.2];  %step width
%so = [0  ,20  ,40  ,60  ,60];  %step arc
so = [0, 20 , 20, 20 ,0];
% sx = [0,  0.3,0.3,0.3,0]; %step length
% sy = [0.2,0.2,0.2,0.2,0.2];  %step width
% so = [0,0,0,0,0];
so = so./180* pi ;
%%

len = length(sx);
px = zeros(1,len+1);   %foot placement
py = zeros(1,len+1);
px_star = zeros(1,len+1);
py_star = zeros(1,len+1);

px(1) = 0;  %initial foot placement
py(1) = 0;

x_bar = zeros(1,len);
y_bar = zeros(1,len);
vx_bar = zeros(1,len);
vy_bar = zeros(1,len);


%step 1
xi = 0;    %set initial position of CoM(x,y)
yi = 0.1;
xi_d = 0;   %initial velocity
yi_d = 0.24;


px_star(1) = px(1);
py_star(1) = py(1);
T = 0;

route_x = zeros(len-1,N);
route_xd = zeros(len-1,N);
route_y = zeros(len-1,N);
route_yd = zeros(len-1,N);

for i = 1:len-1  %caution: the starting index of the list is 1 rather than 0
    %step 4
    T = T + Tsup;
    
    %step 5  modified by 4.60
    px(i+1) = px(i) + cos(so(i)) * sx(i) - sin(so(i)) * (-(-1)^(i) * sy(i));    %caution: the starting index of the list is 1 rather than 0
    py(i+1) = py(i) + sin(so(i)) * sx(i) + cos(so(i)) * (-(-1)^(i) * sy(i));
    
    %step 6 modified by 4.61
    x_bar(i+1) = cos(so(i+1)) * sx(i+1)/2 - sin(so(i+1)) * ((-1)^(i) * sy(i+1)/2);     %the next walk primitive
    y_bar(i+1) = sin(so(i+1)) * sx(i+1)/2 + cos(so(i+1)) * ((-1)^(i) * sy(i+1)/2);
%   vx_bar(i+1) = cos(so(i+1)) * (C+1) / (Tc*S) * x_bar(i+1) - sin(so(i+1)) * (C-1) / (Tc*S) * y_bar(i+1);
%   vy_bar(i+1) = sin(so(i+1)) * (C+1) / (Tc*S) * x_bar(i+1) + cos(so(i+1)) * (C-1) / (Tc*S) * y_bar(i+1);
    vx_bar(i+1) = (C+1) / (Tc*S) * x_bar(i+1);
    vy_bar(i+1) = (C-1) / (Tc*S) * y_bar(i+1);
     

    %step 7
    xd = px(i+1) + x_bar(i+1);  %target state x
    xd_d = vx_bar(i+1);
    yd = py(i+1) + y_bar(i+1);  %target state y
    yd_d = vy_bar(i+1);
    
    %step 8
    D = a*(C-1)^2 + b*(S/Tc)^2;     %equation 4.59
    px_star(i+1) = -a*(C-1)/D * (xd - C*xi - Tc*S*xi_d) ...
        - b*S/(Tc*D) * (xd_d - S/Tc*xi - C*xi_d);
    py_star(i+1) = -a*(C-1)/D * (yd - C*yi - Tc*S*yi_d) ...
        - b*S/(Tc*D) * (yd_d - S/Tc*yi - C*yi_d);
        
    %step 3  integrate equation
    delt_x_r = cos(-so(i+1)) * (xi - px_star(i+1)) - sin(-so(i+1)) * (yi - py_star(i+1));   %relative
    delt_y_r = sin(-so(i+1)) * (xi - px_star(i+1)) + cos(-so(i+1)) * (yi - py_star(i+1));
    xi_d_r = cos(-so(i+1)) * xi_d - sin(-so(i+1)) * yi_d;
    yi_d_r = sin(-so(i+1)) * xi_d + cos(-so(i+1)) * yi_d;
    
    delt_xt = delt_x_r * cosh(Tsup/Tc) + Tc * xi_d_r * sinh(Tsup/Tc);
    delt_xt_d = delt_x_r/Tc * sinh(Tsup/Tc) + xi_d_r * cosh(Tsup/Tc);   
    delt_yt = delt_y_r * cosh(Tsup/Tc) + Tc * yi_d_r * sinh(Tsup/Tc);
    delt_yt_d = delt_y_r/Tc * sinh(Tsup/Tc) + yi_d_r * cosh(Tsup/Tc);
    
    xt = cos(so(i+1)) * delt_xt - sin(so(i+1)) * delt_yt + px_star(i+1); %final state of (i+1)th step
    xt_d = cos(so(i+1)) * delt_xt_d - sin(so(i+1)) * delt_yt_d; 
    yt = sin(so(i+1)) * delt_xt + cos(so(i+1)) * delt_yt + py_star(i+1);  %final state of (i+1)th step
    yt_d = sin(so(i+1)) * delt_xt_d + cos(so(i+1)) * delt_yt_d; 
    
    for j =1:N  %caculate the route alongside
        route_x(i,j) = (xi - px_star(i+1)) * cosh(j*time/Tc) + Tc * xi_d * sinh(j*time/Tc) + px_star(i+1); %final state of (i+1)th step
        route_xd(i,j) = (xi - px_star(i+1))/Tc * sinh(j*time/Tc) + xi_d * cosh(j*time/Tc);
        route_y(i,j) = (yi - py_star(i+1)) * cosh(j*time/Tc) + Tc * yi_d * sinh(j*time/Tc) + py_star(i+1); %final state of (i+1)th step
        route_yd(i,j) = (yi - py_star(i+1))/Tc * sinh(j*time/Tc) + yi_d * cosh(j*time/Tc);

    end

    %update
    xi = xt;
    xi_d = xt_d;
    yi = yt;
    yi_d = yt_d;
end

routex = [route_x(1,:),route_x(2,:),route_x(3,:),route_x(4,:)];
routey = [route_y(1,:),route_y(2,:),route_y(3,:),route_y(4,:)];
figure(1);plot(routex,routey);
hold on;
plot(px_star,py_star,'o');
plot(px,py,'bx');