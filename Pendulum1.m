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

simu_time = 0.7;
step_time = 0.01;
simu_counts = round(simu_time / step_time);
total_ZMP_x_t = [linspace(0.15,0.2,round(simu_counts/7)),linspace(0.2,0.2,round(5*simu_counts/7)),linspace(0.2,0.45,round(1*simu_counts/7))];

pend_x_t         = zeros(simu_counts,1);
ddot_pend_x_t    = zeros(simu_counts,1);
dot_pend_x_t     = zeros(simu_counts,1);
pend_x_t(1) = 0.2;
dot_pend_x_t(1) = 0.1;
pend_ZMP_x_t     = zeros(simu_counts,1);
[ swg_x_t,ddot_swg_x_t,swg_z_t,ddot_swg_z_t ] = FootTrack(simu_time*6/7,step_time,0.6,0.05,10);
for i = 1:simu_counts-1
    if(i  > round(simu_counts/7))
        [ddot_pend_x_t(i+1),pend_ZMP_x_t(i+1)] = ThreeMassCal(sup_x,ddot_sup_x,swg_x_t(i - round(simu_counts/7)),ddot_swg_x_t(i - round(simu_counts/7)), ...
            sup_z,ddot_sup_z,swg_z_t(i - round(simu_counts/7)),ddot_swg_z_t(i - round(simu_counts/7)),...
            total_ZMP_x_t(i),pend_x_t(i));
    else
        [ddot_pend_x_t(i+1),pend_ZMP_x_t(i+1)] = ThreeMassCal(sup_x,ddot_sup_x,swg_x,ddot_swg_x,sup_z,ddot_sup_z,swg_z,ddot_swg_z,...
            total_ZMP_x_t(i),pend_x_t(i));
    end
    
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
h_foot_left = rectangle('Position',[-0.1 0 0.25 0.02],'Curvature',[0.5,1],'EdgeColor','b');
h_foot_right = rectangle('Position',[0.2 0 0.25 0.02],'Curvature',[0.5,1],'EdgeColor','r');
h_foot_mass_left = plot(0,0.08,'bo','MarkerSize',5,'LineWidth',3);
h_foot_mass_right = plot(0.3,0.08,'ro','MarkerSize',5,'LineWidth',3);

h_pend = line([pend_x_t(1),pend_ZMP_x_t(1)],[h,0]);   
axis([-0.3,0.95,0,0.7]);   
ax = gca;
set(gca,'Children',[h_pend_ZMP,h_total_ZMP,h_foot_left,h_foot_right,h_foot_mass_left,h_foot_mass_right,h_pend]);
vid = VideoWriter('myPeaks2.avi');
vid.Quality = 100;
vid.FrameRate = 15;
open(vid);

for i = 2:simu_counts    %the first step
    set(h_pend_ZMP,'xdata',pend_x_t(i));
    set(h_total_ZMP,'xdata',total_ZMP_x_t(i));
    set(h_pend,'xdata',[pend_x_t(i),pend_ZMP_x_t(i)]);
    set(h_pend,'Color',[.1,.1,0.1]);
    if(i  > round(simu_counts/7))
        set(h_foot_left,'Position',[swg_x_t(i - round(simu_counts/7))-0.1,swg_z_t(i - round(simu_counts/7)),0.25,0.02]);
        set(h_foot_mass_left,'xdata',swg_x_t(i - round(simu_counts/7)),'ydata',swg_z_t(i - round(simu_counts/7)) + 0.08);
    end
    
    drawnow;
    writeVideo(vid, getframe(gcf));
end
for i = 1: 30
     drawnow;
    writeVideo(vid, getframe(gcf));
end
close(vid);

winopen('myPeaks2.avi')
