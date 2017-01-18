clear;
clc;
close all;
h = 0.6;
sup_x      = 0.3;   sup_z      = 0.08; 
ddot_sup_x = 0;   ddot_sup_z = 0;


swg_x      = 0;   swg_z      = 0.08;   
ddot_swg_x = 0;   ddot_swg_z = 0;

total_ZMP_x = -0.1;%-0.08;
pend_x = 0.0       ;

simu_time = 0.75;
step_time = 0.01;
round = round(simu_time / step_time);
total_ZMP_x_t = [linspace(-0.1,-0.1,round/3),linspace(-0.1,0.2,round/3),linspace(0.2,0.45,round/3)];

pend_x_t         = zeros(round,1);
ddot_pend_x_t    = zeros(round,1);
dot_pend_x_t     = zeros(round,1);
dot_pend_x_t(1) = 0;
pend_ZMP_x_t     = zeros(round,1);
for i = 1:round-1
    [ddot_pend_x_t(i+1),pend_ZMP_x_t(i+1)] = ThreeMassCal(sup_x,ddot_sup_x,swg_x,ddot_swg_x,sup_z,ddot_sup_z,swg_z,ddot_swg_z,...
    total_ZMP_x_t(i),pend_x_t(i));
    dot_pend_x_t(i+1) = dot_pend_x_t(i) + ddot_pend_x_t(i+1)*step_time;
    pend_x_t(i+1) = pend_x_t(i) +dot_pend_x_t(i+1)*step_time;
end
% figure(1);
% hold on;
% plot(ddot_pend_x_t,'k');hold off;
% figure(2);hold on;
% plot(pend_x_t,'r');
% plot(pend_ZMP_x_t,'b');
%%draw the model
figure;hold on;
h_pend_ZMP = plot(pend_x_t(1),h,'bo','MarkerSize',20);
h_total_ZMP = plot(total_ZMP_x_t(1),0,'ro','MarkerSize',5,'MarkerFaceColor','r');
h_foot_right = rectangle('Position',[-0.1 0 0.25 0.02]);
h_foot_right = rectangle('Position',[0.2 0 0.25 0.02]);

h_pend = line([pend_x_t(1),pend_ZMP_x_t(1)],[h,0]);   
axis([-1,2,0,0.7]);   
% ax = gca;
% set(gca,'Children',[h_pend_ZMP,h_total_ZMP,h_pend]);

for i = 2:round    %the first step
    set(h_pend_ZMP,'xdata',pend_x_t(i));
    set(h_total_ZMP,'xdata',total_ZMP_x_t(i));
    set(h_pend,'xdata',[pend_x_t(i),pend_ZMP_x_t(i)]);
    set(h_pend,'Color',[.1,.1,0.1]);
    drawnow;
end



function [ ddot_pend_x, pend_ZMP_x] = ThreeMassCal(sup_x,ddot_sup_x,swg_x,ddot_swg_x,sup_z,ddot_sup_z,swg_z,ddot_swg_z,total_ZMP_x,pend_x)
    g = 10;
    h = 0.6;
    m_sup = 3;
    m_swg = 3;
    m_pend = 24;
    m_feet = m_sup + m_swg;
    m_total = m_feet + m_pend;
    M_feet = m_sup*(sup_x -total_ZMP_x)*(g + ddot_sup_z) ...
                - m_sup * ddot_sup_x * (sup_z - 0) ...
                + m_swg * (swg_x - total_ZMP_x)*(g + ddot_swg_z)  ...
                - m_swg * ddot_swg_x * (swg_z - 0);
            

    feet_ZMP_x   = M_feet / m_feet / g + total_ZMP_x;
    
    pend_ZMP_x = m_total / m_pend * total_ZMP_x - m_feet / m_pend *feet_ZMP_x;

     ddot_pend_x = g/h*(pend_x - pend_ZMP_x);
end