function [easy_hr_expt hard_hr_expt expt_names] = plotAVHR_expt(mouse,imouse,use_ind,mice_hr_fig,colors_mice,mouse_name);
    nExp = length(use_ind);
    start = 1;
    for iexp = 1:nExp
        input = mouse(imouse).input(use_ind(iexp));
        
        missedIx = strcmp(input.trialOutcomeCell, 'ignore');
        successIx = strcmp(input.trialOutcomeCell, 'success');
        earliesIx = strcmp(input.trialOutcomeCell, 'failure');
        FAIx = strcmp(input.catchTrialOutcomeCell, 'FA');
        CRIx = strcmp(input.catchTrialOutcomeCell, 'CR');
        earliesIx(FAIx) = 0;
        earliesIx(CRIx) = 0;

        catchDirectionDeg = chop(double(celleqel2mat_padded(input.tCatchGratingDirectionDeg)),2);
        catchOris = unique(catchDirectionDeg(~isnan(catchDirectionDeg)));
        targetDirectionDeg = chop(double(celleqel2mat_padded(input.tGratingDirectionDeg)),2);
        
        if ~isempty(catchOris)
            easy_deg = catchOris(end);
            hard_deg = catchOris(1);

            easy_inv_ind = catchDirectionDeg == easy_deg;
            hard_inv_ind = catchDirectionDeg == hard_deg;
            easy_val_ind = targetDirectionDeg == easy_deg;
            hard_val_ind = targetDirectionDeg == hard_deg;

            easy_inv_hr = sum(easy_inv_ind(FAIx))./(sum(easy_inv_ind(FAIx))+sum(easy_inv_ind(CRIx)));
            hard_inv_hr = sum(hard_inv_ind(FAIx))./(sum(hard_inv_ind(FAIx))+sum(hard_inv_ind(CRIx)));
            easy_val_hr = sum(easy_val_ind(successIx))./(sum(easy_val_ind(successIx))+sum(easy_val_ind(missedIx)));
            hard_val_hr = sum(hard_val_ind(successIx))./(sum(hard_val_ind(successIx))+sum(hard_val_ind(missedIx)));
            
                figure(mice_hr_fig);
                subplot(1,2,1)
                hold on
                h = plot([1 2],[easy_val_hr easy_inv_hr],'o');
                h.MarkerEdgeColor = [1 1 1];
                h.MarkerFaceColor = colors_mice(imouse,:);
                h.LineStyle = '-';
                h.Color = [0.5 0.5 0.5];

                figure(mice_hr_fig);
                subplot(1,2,2)
                hold on
                h = plot([1 2],[hard_val_hr hard_inv_hr],'o');
                h.MarkerEdgeColor = [1 1 1];
                h.MarkerFaceColor = colors_mice(imouse,:);
                h.LineStyle = '-';
                h.Color = [0.5 0.5 0.5];

            easy_hr_expt(:,start) = cat(1,easy_val_hr,easy_inv_hr);
            hard_hr_expt(:,start) = cat(1,hard_val_hr,hard_inv_hr);
            expt_names(start) = mouse_name;
            start = start+1;
        end


    end
end