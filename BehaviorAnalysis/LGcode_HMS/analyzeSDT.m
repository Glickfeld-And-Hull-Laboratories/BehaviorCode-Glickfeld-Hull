rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
pv = behavParamsPV;

itask = 1;
pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '_mouse.mat']);
load(pn);

for imouse = find(pv.chr2_mat == 1)
	for ipos = [find(pv.pos_mat == 35) find(pv.pos_mat == -35)]
    	for ipow = 1:length(pv.power_mat)
            diffc_all = [];
            diffdp_all = [];
            pctFA_all = [];
         	for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                fName = fullfile(rc.pathStr,['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                                load(fName);
                                %backwards compat:
                                if isfield(input,'tTotalReqHoldTimeMs')
                                    reqHold = cell2mat(input.tTotalReqHoldTimeMs);
                                elseif isfield(input,'reqHoldTimeMs')
                                    reqHold = cell2mat(input.reqHoldTimeMs);
                                end
                                outS = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).outS;
                                [max_hold ind] = max(cell2mat(input.holdTimesMs),[],2);
                                max_req = reqHold(ind)-550;
                                crop_ind = find(reqHold<max_req);
                                % crop_ind = find(reqHold<(0.7*max_hold));
                                cons = intersect(outS.intensitiesC{1}, outS.intensitiesC{2});
                                bad_con = [];
                                for icon = 1:length(cons)
                                    for iblock = 1:2
                                        block_ind = intersect(crop_ind, find(outS.block2V==(iblock-1)));
                                        con_ind = intersect(block_ind,find(outS.intensityV == cons(icon)));
                                        if length(con_ind) < 10
                                            bad_con = [bad_con icon];
                                        end
                                    end
                                end
                                if length(bad_con)>1
                                    cons(bad_con) = [];
                                end
                                for iblock = 1:2
                                    hits = zeros(1,length(cons));
                                    misses = zeros(1,length(cons));
                                    CRs = zeros(1,length(cons));
                                    FAs = zeros(1,length(cons));
                                    holds = zeros(1,length(cons));
                                    block_ind = intersect(crop_ind, find(outS.block2V==(iblock-1)));
                                    hold_ind = intersect(block_ind, find(cell2mat(input.holdTimesMs)>400));
                                    pctFA = sum(outS.earlyIx(hold_ind))./length(hold_ind);
                                    for icon = 1:length(cons)
                                        con_ind = intersect(block_ind,find(outS.intensityV == cons(icon)));
                                        hits(1,icon) = sum(outS.successIx(con_ind));
                                        misses(1,icon) = sum(outS.missedIx(con_ind));
                                        tempCR = 0;
                                        tempFA = 0;
                                        reqHolds = [];
                                        for ind = con_ind
                                            if outS.earlyIx(ind)==0
                                                reqHolds = [reqHolds reqHold(ind)];
                                                long_hold_inds = find(reqHold(block_ind)>(reqHold(ind)+550));
                                                least_hold_inds = find(cell2mat(input.holdTimesMs)>(reqHold(ind)));
                                                for icatch = intersect(long_hold_inds, least_hold_inds)
                                                    holdT = input.holdTimesMs{icatch};
                                                    if holdT>(reqHold(ind)+550)
                                                        tempCR = tempCR+1;
                                                    else
                                                        tempFA = tempFA+1;
                                                    end
                                                end
                                            end
                                        end
                                        FAs(1,icon) = tempFA;
                                        CRs(1,icon) = tempCR;
                                        holds(1,icon) = mean(reqHolds,2);
                                    end
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).FAs = FAs;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).CRs = CRs;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).Misses = misses;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).Hits = hits;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).cons = cons;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).holds = holds;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).pctFA = pctFA;
                                    %code from NS to calculate dp and c
                                    hitRate = hits./(hits+misses);
                                    faRate = FAs./(FAs+CRs);
                                    [i j] = find(hitRate==1);
                                    if length(i)>0
                                        for ii = 1:length(i)
                                            hitRate(i(ii),j(ii)) = 1-0.5./(hits(i(ii),j(ii))+misses(i(ii),j(ii)));
                                        end
                                    end
                                    [i j] = find(hitRate==0);
                                    if length(i)>0
                                        for ii = 1:length(i)
                                            hitRate(i(ii),j(ii)) = 0.5./(hits(i(ii),j(ii))+misses(i(ii),j(ii)));
                                        end
                                    end
                                    [i j] = find(faRate==1);
                                    if length(i)>0
                                        for ii = 1:length(i)
                                            faRate(i(ii),j(ii)) = 0.5./(FAs(i(ii),j(ii))+CRs(i(ii),j(ii)));
                                        end
                                    end
                                    [i j] = find(faRate==0);
                                    if length(i)>0
                                        for ii = 1:length(i)
                                            faRate(i(ii),j(ii)) = 1-0.5./(FAs(i(ii),j(ii))+CRs(i(ii),j(ii)));
                                        end
                                    end
                                    dp = norminv(hitRate)-norminv(faRate);
                                    c = -norminv(faRate);
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).hitRate = hitRate;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).faRate = faRate;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).c = c;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(iblock).dp = dp;
                                end
                                meandiffc = nanmean((mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(1).c-mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(2).c),2);
                                meandiffdp = nanmean((mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(1).dp-mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(2).dp),2);
                                fractFA = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(2).pctFA./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT(1).pctFA;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT_meandiffc = meandiffc;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT_meandiffdp = meandiffdp;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).SDT_fractFA = fractFA;
                                diffc_all = [diffc_all meandiffc];
                                diffdp_all = [diffdp_all meandiffdp];
                                pctFA_all = [pctFA_all fractFA];
                            end
                        end
                    end
                end
            end
            if length(diffc_all)>0
                mouse(imouse).task(itask).pos(ipos).pow(ipow).c_all = diffc_all;
                mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_all = diffdp_all;
                mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_all = pctFA_all;
                mouse(imouse).task(itask).pos(ipos).pow(ipow).c_avg = mean(diffc_all,2);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_avg = mean(diffdp_all,2);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_avg = mean(pctFA_all,2);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).c_sem = std(diffc_all,[],2)./sqrt(length(diffc_all));
                mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_sem = std(diffdp_all,[],2)./sqrt(length(diffdp_all));
                mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_sem = std(pctFA_all,[],2)./sqrt(length(pctFA_all));
            end
        end
    end
end
for imouse = find(pv.chr2_mat == 1)
	for ipos = [find(pv.pos_mat == 35) find(pv.pos_mat == -35)]
    	for ipow = 1:length(pv.power_mat)
         	for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_all)>0
                                    col_mat = strvcat('b','k','g','k','r');
                                    if itask == 1;
                                        if ipos == 1;
                                            subplot(2,3,1)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('Dprime Contra')
                                            subplot(2,3,2)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).c_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).c_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('C Contra')
                                            subplot(2,3,3)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('FA Contra')
                                        elseif ipos == 4;
                                            subplot(2,3,4)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).dp_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('Dprime Ipsi')
                                            subplot(2,3,5)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).c_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).c_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('C Ipsi')
                                            subplot(2,3,6)
                                            errorbar(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_avg, mouse(imouse).task(itask).pos(ipos).pow(ipow).pctFA_sem, col_mat(imouse));
                                            hold on
                                            xlim([-3 1])
                                            title('FA Ipsi')
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


