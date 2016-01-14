rc = behavConstsAttn;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
pv = behavParamsAttn;

xs = 0:.01:1;
for imouse = 1:length(pv.mouse_mat)
    pos_ind = find(abs(diff(mouse(imouse).posN)));
    modelFun = mouse(imouse).modelFun;
    for itask = 1:2
        task_ind = intersect(pos_ind, find(mouse(imouse).taskN==itask));
        figure;
        for iexp = task_ind(2:end)
            if mouse(imouse).expt(iexp).well_fit
                if mouse(imouse).expt(iexp-1).well_fit
                    break
                end
            end
        end
        for ipos = 1:2
            subplot(1,2,ipos)
            plot(mouse(imouse).expt(iexp-2+ipos).intensityX, mouse(imouse).expt(iexp-2+ipos).fractCorrY, 'ok')
            hold on;
            pctCorr = mouse(imouse).expt(iexp-2+ipos).block(2).binofit(1).pctCorrect;
            err = [pctCorr-mouse(imouse).expt(iexp-2+ipos).block(2).binofit(1).ci95(1) mouse(imouse).expt(iexp-2+ipos).block(2).binofit(1).ci95(2)-pctCorr];
            errorbar(mouse(imouse).expt(iexp-2+ipos).block(2).intensities(1), pctCorr, err(1), err(2), 'or');
            hold on;
            p =  mouse(imouse).expt(iexp-2+ipos).coefEsts;
            plot(xs, modelFun(p,xs),'-k')
            hold on
            title(['Expt ' num2str(iexp-2+ipos) ' Block 1 Position: ' num2str(mouse(imouse).expt(iexp-2+ipos).pos)])
            set(gca, 'xscale', 'log');
        end
        suptitle(['Mouse ' num2str(pv.mouse_mat(:,imouse)) ' Task ' num2str(itask)])
    end
end

xs = 0:.01:1;
col_mat = {[0 0 0], [0.5 0.5 0.5]; [0 0 1], [0 1 1]}; 
for imouse = 1:length(pv.mouse_mat)
    modelFun = mouse(imouse).modelFun;
    for itask = 1:2
        task_ind = find(mouse(imouse).taskN==itask);
        figure;
        start = 1;
        for iexp = task_ind
            if iexp+1 <= task_ind(end)
                if mouse(imouse).expt(iexp).pos ~= mouse(imouse).expt(iexp+1).pos
                    if mouse(imouse).expt(iexp).well_fit
                        if mouse(imouse).expt(iexp+1).well_fit
                            subplot(5,5,start)
                            for ipos = 1:2
                                pos = mouse(imouse).expt(iexp-1+ipos).pos;
                                x(ipos) = plot(mouse(imouse).expt(iexp-1+ipos).intensityX, mouse(imouse).expt(iexp-1+ipos).fractCorrY, 'o');
                                set(x(ipos),'Color', col_mat{pos,1})
                                hold on;
                                pctCorr = mouse(imouse).expt(iexp-1+ipos).block(2).binofit(1).pctCorrect;
                                err = [pctCorr-mouse(imouse).expt(iexp-1+ipos).block(2).binofit(1).ci95(1) mouse(imouse).expt(iexp-1+ipos).block(2).binofit(1).ci95(2)-pctCorr];
                                y(ipos) = errorbar(mouse(imouse).expt(iexp-1+ipos).block(2).intensities(1), pctCorr, err(1), err(2), 'o');
                                set(y(ipos), 'Color', col_mat{pos,2})
                                hold on;
                                p =  mouse(imouse).expt(iexp-1+ipos).coefEsts;
                                z(ipos) = plot(xs, modelFun(p,xs),['-' col_mat{pos,1}]);
                                set(z(ipos), 'Color', col_mat{pos,1})
                                hold on
                                set(gca, 'xscale', 'log');
                            end
                            title(['Expt ' num2str(iexp) '-' num2str(iexp+1)])
                            start = start+1;
                        end
                    end
                end
            end
        end
        suptitle(['Mouse ' num2str(pv.mouse_mat(:,imouse)) ' Task ' num2str(itask)])
    end
end