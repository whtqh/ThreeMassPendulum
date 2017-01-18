clear;
clc;
close all;

%constant
g = 9.8;     %gravity
time = 0.005;  %the simulation time

a = 1; %positive weight of the evaluation function 4.58
b = 10;

%input
zc = 0.6;       %height
Tc = sqrt(zc/g);

Tsup = 0.9;     %support period
N = round(Tsup/time);
C = cosh(Tsup/Tc);
S = sinh(Tsup/Tc);
count = 41;     %Step Num

%walking parameters
% sx = [0,  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0]; %step length
% sy = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];  %step width
sx = zeros(1,count);
sx(:) = 0.1;
sx(1) = 0;sx(end) = 0;
sy = zeros(1,count);
sy(:) = 0.2;
% sx = [0,  0.2,0.2,0.2,0]; %step length
% sy = [0.2,0.3,0.1,0.3,0.2];  %step width

%%
len = length(sx);
px = zeros(1,len+1);   %foot placement
py = zeros(1,len+1);
px_star = zeros(1,len+1);   %foot placement
py_star = zeros(1,len+1);
px(1) = 0;  %initial foot placement
py(1) = 0;

x_bar = zeros(1,len);
y_bar = zeros(1,len);
vx_bar = zeros(1,len);
vy_bar = zeros(1,len);

xd = zeros(1,len);
xd_d = zeros(1,len);
yd = zeros(1,len);
yd_d = zeros(1,len);
%step 1
xi = 0;    %set initial position of CoM(x,y)
yi = 0.1;
xi_d = 0;   %initial velocity
yi_d = 0.4;


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
    
    %step 5
    px(i+1) = px(i) + sx(i);    %caution: the starting index of the list is 1 rather than 0
    py(i+1) = py(i) -(-1)^(i) * sy(i);
    
    %step 6
    x_bar(i+1) = sx(i+1)/2;     %the next walk primitive
    y_bar(i+1) = (-1)^(i) * sy(i+1)/2;
    vx_bar(i+1) = (C+1) / (Tc*S) * x_bar(i+1);
    vy_bar(i+1) = (C-1) / (Tc*S) * y_bar(i+1);
     
    %step 7
    xd(i+1) = px(i+1) + x_bar(i+1);  %target state x
    xd_d(i+1) = vx_bar(i+1);
    yd(i+1) = py(i+1) + y_bar(i+1);  %target state y
    yd_d(i+1) = vy_bar(i+1);
    
    %step 8
    D = a*(C-1)^2 + b*(S/Tc)^2;     %equation 4.59
    px_star(i+1) = -a*(C-1)/D * (xd(i+1) - C*xi - Tc*S*xi_d) ...
        - b*S/(Tc*D) * (xd_d(i+1) - S/Tc*xi - C*xi_d);
    py_star(i+1) = -a*(C-1)/D * (yd(i+1) - C*yi - Tc*S*yi_d) ...
        - b*S/(Tc*D) * (yd_d(i+1) - S/Tc*yi - C*yi_d);
        
    %step 3  integrate equation
    xt = (xi - px_star(i+1)) * cosh(Tsup/Tc) + Tc * xi_d * sinh(Tsup/Tc) + px_star(i+1); %final state of (i+1)th step
    xt_d = (xi - px_star(i+1))/Tc * sinh(Tsup/Tc) + xi_d * cosh(Tsup/Tc);
    yt = (yi - py_star(i+1)) * cosh(Tsup/Tc) + Tc * yi_d * sinh(Tsup/Tc) + py_star(i+1); %final state of (i+1)th step
    yt_d = (yi - py_star(i+1))/Tc * sinh(Tsup/Tc) + yi_d * cosh(Tsup/Tc);
    
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

route_x = route_x';
route_y = route_y';
routex = route_x(:)';
routey = route_y(:)';
% routex = [route_x(1,:),route_x(2,:),route_x(3,:),route_x(4,:),route_x(5,:),route_x(6,:),route_x(7,:),route_x(8,:)];
% routey = [route_y(1,:),route_y(2,:),route_y(3,:),route_y(4,:),route_y(5,:),route_y(6,:),route_y(7,:),route_y(8,:)];
figure(1);plot(routex,routey);
hold on;
plot(px_star,py_star,'o');
plot(px,py,'bx');
plot(xd,yd,'ro')

