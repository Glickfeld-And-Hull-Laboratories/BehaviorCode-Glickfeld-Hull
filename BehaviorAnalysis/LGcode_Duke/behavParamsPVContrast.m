function pv = behavParamsPVContrast
    pv.task_mat = strvcat('Con', 'Ori');  
    pv.pos_mat = [35];
    pv.mouse_mat = [502 503];
    pv.chr2_mat = [1 1];
    pv.arch_mat = [0 0];
    pv.power_mat = [0.25];
    pv.basecon_mat = [0 0.5 1];
    pv.mark_str = strvcat('o', 's');
    pv.col_str = strvcat('b', 'r');
end