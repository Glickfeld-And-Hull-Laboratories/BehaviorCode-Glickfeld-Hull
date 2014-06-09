function plotAllAttnData(itask)    

    mouse = createExptStructAttn(itask);
    rc = behavConstsAttn;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAttn;

    a = [2 1; 1 2];
    b = [1 1; 2 2];
    block_mat = cat(3, a, b);
    
    figure;
    for imouse = 1:length(pv.mouse_mat)
        subplot(length(pv.mouse_mat),2,1+((imouse-1)*2))
        for ipos = 1:length(pv.pos_mat)
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
                if mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).well_fit == 1
                    plot(ipos, mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).thresh, 'ok')
                    hold on
                end
            end
        end
        xlim([0 3])
        ylim([0 .3])
        title(['Threshold contrast at each position for ' num2str(pv.mouse_mat(imouse))])
        ylabel('Threshold contrast')
        xlabel('Position')
    end
    
    col_mat = strvcat('k', 'r');
    for imouse = 1:length(pv.mouse_mat)
        subplot(length(pv.mouse_mat),2,2+((imouse-1)*2))
        for ipos = 1:2
            for iblock = 1:2
                ub = mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).ci95(1) - mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).pctCorrect;
                lb = mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).pctCorrect - mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).ci95(2);                
                plot(ipos, mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).pctCorrect, ['o' col_mat(block_mat(iblock,2,ipos))])
                hold on
                errorbarxy(ipos, mouse(imouse).task(itask).pos(block_mat(iblock,1,ipos)).block(block_mat(iblock,2,ipos)).pctCorrect, [], lb,[],ub,col_mat(block_mat(iblock,2,ipos)), col_mat(block_mat(iblock,2,ipos)));
                hold on
            end
            text(ipos,0.1, ['P =' num2str(chop(mouse(imouse).task(itask).pos(ipos).p,2))]);
        end
        xlim([0 3])
        ylim([0 1]) 
    end

    title(['Attn task ' num2str(itask) ' summary for ' num2str(pv.mouse_mat(imouse))])
    ylabel('Fraction Correct')
    xlabel('Position')
    pn = fullfile(rc.fitOutputSummary, ['Attn-Task' num2str(itask) '-summary-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

            
        
    col_mat = strvcat('k', 'r');
    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            figure;
            iexp = 1;
            for ind = mouse(imouse).task(itask).pos(ipos).ind
                for iblock = 1:2
                    for icon = 1:length(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities)
                        subplot(ceil(sqrt(length(mouse(imouse).task(itask).pos(ipos).ind))),ceil(sqrt(length(mouse(imouse).task(itask).pos(ipos).ind))), iexp) 
                        plot(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon), mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).pctCorrect,['o' col_mat(iblock)]);
                        hold on
                        ub = mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).ci95(2) - mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).pctCorrect;
                        lb = mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).pctCorrect - mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).ci95(1); 
                        errorbar(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon), mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).pctCorrect, lb,ub, col_mat(iblock))
                        hold on
                    end
                end
                xlim([0 1])
                ylim([0 1])
                title([' P = ' num2str(chop(mouse(imouse).task(itask).pos(ipos).expt(iexp).p(1),2))]);
                iexp = iexp +1;
            end
            pct = num2str(100-(1200/80));
            suptitle(['Position ' pct '% at ' num2str(pv.pos_mat(ipos)) 'deg for ' num2str(pv.mouse_mat(imouse))])
            y = pv.pos_mat(ipos);
            if y >0
                pos = num2str(y);
            else
                pos = ['neg' num2str(y)];
            end
            pn = fullfile(rc.fitOutputSummary, ['Attn-Task' num2str(itask) '-Pos-' pos 'allexpt-summary-i' num2str(pv.mouse_mat(imouse)) '-' date '.pdf']);
            exportfig_print(gcf, pn, 'FileFormat', 'pdf');
        end
    end
    
    figure;
    for imouse = 1:length(pv.mouse_mat)
        subplot(1,length(pv.mouse_mat),imouse)
        for iblock = 1:2
            plot(1:2, [mouse(imouse).task(itask).block(iblock).half(1).pctCorrect mouse(imouse).task(itask).block(iblock).half(2).pctCorrect],['-' col_mat(iblock)]);
            hold on
        end
        xlim([0 3])
        ylim([0 1])
        title(['Mouse ' num2str(pv.mouse_mat(imouse))])
        xlabel('Half')
        ylabel('% correct')
    end
    pn = fullfile(rc.fitOutputSummary, ['Attn-Task' num2str(itask) '-byhalf-allexpt-summary-' date '.pdf']);
            exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end
    