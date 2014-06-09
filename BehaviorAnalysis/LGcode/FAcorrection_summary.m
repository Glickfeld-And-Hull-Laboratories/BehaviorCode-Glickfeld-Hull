rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
pv = behavParamsPV;

for itask = 1:4;
pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '_mouse.mat']);
load(pn);

pct_all = [];
if itask == 1
    for imouse = find(pv.chr2_mat == 1)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                    outS = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).outS;
                                    rem = sum(round(outS.estFalseCorrPerImage(:)));
                                    tot = size(outS.whichTrials,2);
                                    pct = rem/tot;
                                    pct_all = [pct_all pct];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
elseif itask == 2
    for imouse = find(pv.chr2_mat == 1)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                icon = find(pv.basecon_mat == 1);
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1;
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                    outS = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).outS;
                                    rem = sum(round(outS.estFalseCorrPerImage(:)));
                                    tot = size(outS.whichTrials,2);
                                    pct = rem/tot;
                                    pct_all = [pct_all pct];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    elseif itask > 2
    for imouse = 1
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                icon = find(pv.basecon_mat == 0.5);
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1;
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                    outS = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).outS;
                                    rem = sum(round(outS.estFalseCorrPerImage(:)));
                                    tot = size(outS.whichTrials,2);
                                    pct = rem/tot;
                                    pct_all = [pct_all pct];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

FAcorr(itask).pct_all = pct_all;
FAcorr(itask).pct_avg = mean(pct_all,2);
FAcorr(itask).pct_med = median(pct_all,2);
FAcorr(itask).pct_range = [min(pct_all,[],2) max(pct_all,[],2)];
end

pn = fullfile(rc.fitOutputSummary, ['Chr2-FAcorr.mat']);
        save(pn, 'FAcorr');
        
FAcorr_all = [FAcorr(1).pct_all FAcorr(2).pct_all FAcorr(3).pct_all FAcorr(4).pct_all]
FAcorr_range = [min(FAcorr_all,[],2) max(FAcorr_all,[],2)]
FAcorr_avg = mean(FAcorr_all,2)
FAcorr_med = median(FAcorr_all,2)