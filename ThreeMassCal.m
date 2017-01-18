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