clear;
clc;
close all;

%constant
g = 9.8;     %gravity
time = 0.003;  %the simulation time

%input
z = 0.8;       %height
Tc = sqrt(z/g);

x_0 = 1;    %initial position
xd0 = 9;    %initial velocity
xf0 = x_0;  
xf1 = 1.5;
xf2 = x_0;

s1 = 3;     %the length of first step 
s2 = 2;     %the length of second step

E_0 = 1/2 * xd0^2 -g/(2*z) * x_0^2;  % initial orbital energy
E_1 = 1/2 * xd0^2 - g/(2*z) * (s1 - x_0)^2;
xd1 = sqrt(2*E_1 + g/z * xf1^2);
E_2 = 1/2 * xd1^2 - g/(2*z) * (s2-xf1)^2;
xd2 = sqrt(2*E_2 + g/z * xf2^2);

%eq 4.11
torq1 = Tc * log((xf0 - s1 - Tc * xd0)/(xf1 - Tc * xd1));    %¦Ó1 is the period that the CoM takes a trip
torq2 = Tc * log((xf1 - s2 - Tc * xd1)/(xf2 - Tc * xd2));    %¦Ó2 is the period that the CoM takes a trip

len1 = round(torq1/time);     %caculate the number of the key points
len2 = round(torq2/time);
route1 = ones(1,len1);        %save the position along the route1 at each key points
velocity1 = ones(1,len1);     %save the velocity along the route1 at each key points
route2 = ones(1,len2);
velocity2 = ones(1,len2);

for i = 1:1:len1    %caculate position and velocity
   route1(i) = (xf0 - s1) * cosh(i*time/Tc) + Tc * xd0 * sinh(i*time/Tc);   
   route1(i) = route1(i) + s1;
   velocity1(i) = (xf0-s1)/Tc* sinh(i*time/Tc) + xd0 * cosh(i*time/Tc);
end
for i = 1:1:len2
   route2(i) = (xf1 - s2) * cosh(i*time/Tc) + Tc * xd1 * sinh(i*time/Tc);
   route2(i) = route2(i) + s2 + s1;
   velocity2(i) = (xf1-s2)/Tc* sinh(i*time/Tc) + xd1 * cosh(i*time/Tc);
end

%%draw the model
h = plot(route1(1),z,'bo','MarkerSize',20);
axis([0,10,0,1]);    

h2 = line([route1(1),s1],[z,0]);     
ax = gca;
set(gca,'Children',[h,h2]);

for i = 2:len1    %the first step
    set(h,'xdata',route1(i));
    set(h2,'xdata',[route1(i),s1]);
    drawnow;
end
for i = 1:len2    %the second step
    set(h,'xdata',route2(i));
    set(h2,'xdata',[route2(i),s1+s2]);
    set(h2,'Color',[.1,.1,0.1]);
    drawnow;
end

figure(2);     %plot velocity mat
velocity = [velocity1,velocity2];
plot(velocity);