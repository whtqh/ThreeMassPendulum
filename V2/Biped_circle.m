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
Zc = zeros(1,4*N);
Z0 = zeros(1,4*N);
for i = 1:4*N
Zc(i) = zc;
end
%walking parameters
sx = [0,  0.25,0.25,0.25,0  ]; %step length
%sx = [0,  0.3,0.3,0.3,0  ]; %step length
sy = [0.2,0.2, 0.2, 0.2, 0.2];  %step width
ro = [20, 20, 40, 60, 60];
%ro = [0, 0, 0, 0, 0];
ro = ro./180.*pi;

%%

len = length(sx);
px = zeros(1,len+1);   %foot placement
py = zeros(1,len+1);
px_star = zeros(1,len+1);   %foot placement
py_star = zeros(1,len+1);
px(1) = 0;  %initial foot placement
py(1) = 0;
sx_r = zeros(1,len);
sy_r = zeros(1,len);
    
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
yi = 0.08;
xi_d = 0;   %initial velocity
yi_d = 0.5;
xd(1) = xi;
yd(1) = yi;

px_star(1) = px(1);
py_star(1) = py(1);
T = 0;

route_x = zeros(len-1,N);
route_xd = zeros(len-1,N);
route_y = zeros(len-1,N);
route_yd = zeros(len-1,N);
foot_x = zeros(1,N*4);
foot_y = zeros(1,N*4);
for i =1:len
    sx_r(i) = sx(i) * cos(ro(i)) - (-(-1)^(i) * sy(i)) * sin(ro(i));
    sy_r(i) = sx(i) * sin(ro(i)) + (-(-1)^(i) * sy(i)) * cos(ro(i));    
end
for i = 1:len-1  %caution: the starting index of the list is 1 rather than 0
    %step 4
    T = T + Tsup;
    
    %step 5
    px(i+1) = px(i) + sx_r(i);    %caution: the starting index of the list is 1 rather than 0
    py(i+1) = py(i) + sy_r(i);
    
    %step 6

    x_bar(i) = sx_r(i+1)/2;     %the next walk primitive
    y_bar(i) = sy_r(i+1)/2;
% ������
%     vx_bar(i) = (C+1) / (Tc*S) * x_bar(i);
%     vy_bar(i) = (C-1) / (Tc*S) * y_bar(i);
% The final bug.
    vx_bar_r =  (C+1)/(Tc*S)*sx(i+1)/2;
    vy_bar_r =  (C-1)/(Tc*S)*(-(-1)^(i+1) * sy(i+1))/2;
    vx_bar(i) = vx_bar_r * cos(ro(i+1)) - vy_bar_r*sin(ro(i+1));
    vy_bar(i) = vx_bar_r * sin(ro(i+1)) + vy_bar_r*cos(ro(i+1));
    
    %step 7
    xd(i+1)   = px(i+1) + x_bar(i);  %target state x
    xd_d(i+1) = vx_bar(i);
    yd(i+1)   = py(i+1) + y_bar(i);  %target state y
    yd_d(i+1) = vy_bar(i);
    
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
        foot_x(N*(i-1)+j) = px_star(i+1);
        foot_y(N*(i-1)+j) = py_star(i+1);
    end

    
    
    %update
    xi = xt;
    xi_d = xt_d;
    yi = yt;
    yi_d = yt_d;
end
%plot 2D Figure
routex = [route_x(1,:),route_x(2,:),route_x(3,:),route_x(4,:)];
routey = [route_y(1,:),route_y(2,:),route_y(3,:),route_y(4,:)];
figure(1);plot(routex,routey);
hold on;
plot(px_star,py_star,'o');
plot(px,py,'bx');
plot(xd,yd,'r+');
%Plot 3D Figure
figure(2);
h = plot3(routex,routey,Zc,'b','MarkerSize',20);
hold on;
h = plot3(routex,routey,Z0,'r','MarkerSize',20);
h2 = line([routex(1),foot_x(1)],[routey(1),foot_y(1)],[Zc(1),Z0(1)]);
for i = 1:N/4
line([routex(i*4),foot_x(i*4)],[routey(i*4),foot_y(i*4)],[Zc(i*4),Z0(i*4)],'Color',[.8 .1 .1]);
end
for i = N/4+1:N/2
line([routex(i*4),foot_x(i*4)],[routey(i*4),foot_y(i*4)],[Zc(i*4),Z0(i*4)],'Color',[.1 .1 .8]);
end
for i = N/2+1:3*N/4
line([routex(i*4),foot_x(i*4)],[routey(i*4),foot_y(i*4)],[Zc(i*4),Z0(i*4)],'Color',[.9 .1 .0]);
end
for i = N/4*3+1:N
line([routex(i*4),foot_x(i*4)],[routey(i*4),foot_y(i*4)],[Zc(i*4),Z0(i*4)],'Color',[.1 .1 .9]);
end
