function pv = behavParamsPV
    pv.task_mat = [1 2 3 4];  
    pv.pos_mat = [35 15 0 -35];
    pv.mouse_mat = [101 103 104 108 110];
    pv.chr2_mat = [1 0 1 0 1];
    pv.arch_mat = [0 1 0 1 0];
    pv.power_mat = [10.24 5.12 4 2.56 2 1.28 1 0.64 0.5 0.32 0.25 0.16 0.15 0.125 0.1 0.08 0.0625 0.04 0.03 0.02 0.01 0.015 0.0075 0.005 0.002];
    pv.basecon_mat = [1 0.75 .5 .25 .125 0];
    pv.mark_str = strvcat('o', 's', '+', 'x', '*');
    pv.col_str = strvcat('b', 'c', 'g', 'k', 'r', 'm');
end