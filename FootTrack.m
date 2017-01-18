function [ swg_x_t,ddot_swg_x_t,swg_z_t,ddot_swg_z_t ] = FootTrack(T,step_time,L_Step,H_Step,ddot_x_1)
    t1 = T/2 - sqrt(ddot_x_1^2*T^2 - 4*ddot_x_1*L_Step)/2/ddot_x_1;
    ddot_z_1 = 4*H_Step/(T*t1);
    ddot_z_2 = - ddot_z_1*t1/(T/2-t1);
    
    simu_counts = round(T / step_time);
    
    swg_x_t = zeros(simu_counts,1);
    swg_z_t = zeros(simu_counts,1);
    ddot_swg_x_t = zeros(simu_counts,1);
    ddot_swg_z_t = zeros(simu_counts,1);
    for i = 1:simu_counts
        t_now = i * step_time;
        if t_now < t1
            swg_x_t(i) = 1/2 * ddot_x_1 * t_now^2;
            swg_z_t(i) = 1/2 * ddot_z_1 * t_now^2;
            ddot_swg_x_t(i) = ddot_x_1;
            ddot_swg_z_t(i) = ddot_z_1;
        elseif t_now < T - t1
            swg_x_t(i) = 1/2 * ddot_x_1 * (2*t1*t_now - t1^2);
            swg_z_t(i) = 1/2 * ddot_z_1 * t_now^2 ...
                            - 1/2 * (ddot_z_1 - ddot_z_2) * (t_now - t1)^2;
            ddot_swg_x_t(i) = 0;
            ddot_swg_z_t(i) = ddot_z_2;
        else
            swg_x_t(i) = L_Step - 1/2 * ddot_x_1 * (T - t_now)^2;
            swg_z_t(i) = 1/2 * ddot_z_1 * (T - t_now)^2;
            ddot_swg_x_t(i) = - ddot_x_1;
            ddot_swg_z_t(i) =   ddot_z_1;
        end
    end

end