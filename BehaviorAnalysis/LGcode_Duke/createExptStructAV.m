function mouse = createExptStructAV;
rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
av = behavParamsAV;

for imouse = 1:size(av,2);
    mouse_name = av(imouse).mouse;
    disp(mouse_name)
    ind = find(xd.Subject == mouse_name);
    date_mat = [];
    early_mat = [];
    HR_ori_mat = [];
    HR_amp_mat = [];
    catch_mat = zeros(2, length(ind));
    for iexp = 1:length(ind)
        disp(iexp)
        idate = xd.DateStr{ind(iexp)};
        n  = dir(fullfile(rc.pathStr, ['data-i' num2str(mouse_name) '-' idate '-*']));
        if ~isnan(str2num(xd.ChooseMatFile{ind(iexp)}))
            n = n(str2num(xd.ChooseMatFile{ind(iexp)}));
        end
        for ifile = 1:size(n,1)
            if ifile == 1
                input_temp = mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], []);
            else
                input_temp = [input_temp mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], [])];
%             try
%                 input_temp = [input_temp mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], [])];
%             catch
%                 input2 = mwLoadData(fn_mworks, [], []);
%                 inpNames1 = fieldnames(input);
%                 inpNames2 = fieldnames(input2);
%                 inpLong = gt(length(inpNames1),length(inpNames2));
%                 if inpLong == 1
%                     inpPlusInd = ismember(inpNames1,inpNames2);
%                     inpPlus = inpNames1(~inpPlusInd);
%                     for i = 1:length(inpPlus)
%                         input2.(genvarname(inpPlus(i))) = cell(1,length(inpNames1));
%                     end
%                 else
%                     inpPlusInd = ismember(inpNames2,inpNames1);
%                     inpPlus = inpNames2(~inpPlusInd);
%                     for i = 1:length(inpPlus)
%                         input_temp.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames2));
%                     end
%                 end
%                 input_temp = [input_temp input2];
%             end
            end
        end
        if size(n,1)>1
            input_temp = concatenateDataBlocks(input_temp);
        end
        range = str2num(xd.TrialRangeToUse{ind(iexp)});
        if ~isnan(range(1))
            input_temp = trialChopper(input_temp, [range(1) range(end)]);
        end
       date_mat = [date_mat; idate];
        
        failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
        missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
        successIx = strcmp(input_temp.trialOutcomeCell, 'success');
        pctEarly = sum(failureIx,2)./length(failureIx);
        early_mat = [early_mat; pctEarly];
        
        gratingDirectionDeg = celleqel2mat_padded(input_temp.tGratingDirectionDeg);
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);

        oris = unique(gratingDirectionDeg);
        amps = unique(soundAmplitude);
        maxOriTrials = find(gratingDirectionDeg == max(oris,[],2));
        maxAmpTrials = find(soundAmplitude == max(amps,[],2));
        pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
        pctCorr_maxAmp = sum(successIx(maxAmpTrials),2)./(sum(successIx(maxAmpTrials),2)+sum(missedIx(maxAmpTrials),2));
        HR_ori_mat = [HR_ori_mat; pctCorr_maxOri];
        HR_amp_mat = [HR_amp_mat; pctCorr_maxAmp];
        
        if ~isfield(input_temp, 'catchTrialOutcomeCell')
            for trN = 1:length(input_temp.trialOutcomeCell)
                if input_temp.tShortCatchTrial{trN}
                    if input_temp.tFalseAlarm{trN}
                        input_temp.catchTrialOutcomeCell{trN} = 'FA';
                    end
                    if isfield(input_temp, 'cCatchOn')
                        if isempty(input_temp.cCatchOn{trN})
                            input_temp.cCatchOn{trN} = NaN;
                            input_temp.catchTrialOutcomeCell{trN} = 'failure';
                        end
                        if (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})>input_temp.nFramesReact
                            input_temp.catchTrialOutcomeCell{trN} = 'CR';
                        end
                        if (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})<input_temp.nFramesTooFast
                            input_temp.catchTrialOutcomeCell{trN} = 'failure';
                        end
                    else
                        if isempty(input_temp.tCatchTimeMs{trN})
                            input_temp.cCatchOn{trN} = NaN;
                            input_temp.catchTrialOutcomeCell{trN} = 'failure';
                        end
                        if (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})>input_temp.reactTimeMs
                            input_temp.catchTrialOutcomeCell{trN} = 'CR';
                        end
                        if (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})<input_temp.tooFastTimeMs
                            input_temp.catchTrialOutcomeCell{trN} = 'failure';
                        end
                    end
                else
                    input_temp.catchTrialOutcomeCell{trN} = 'NaN';
                end
            end
        end
        if ~isfield(input_temp, 'tSoundCatchAmplitude')
            input_temp.tSoundCatchAmplitude = cell(size(input_temp.trialOutcomeCell));
        end
        mouse(imouse).input(iexp).trialOutcomeCell = input_temp.trialOutcomeCell;
        mouse(imouse).input(iexp).tGratingDirectionDeg = input_temp.tGratingDirectionDeg;
        mouse(imouse).input(iexp).tSoundTargetAmplitude = input_temp.tSoundTargetAmplitude;
        mouse(imouse).input(iexp).tCatchGratingDirectionDeg = input_temp.tCatchGratingDirectionDeg;
        mouse(imouse).input(iexp).tSoundCatchAmplitude = input_temp.tSoundCatchAmplitude;
        mouse(imouse).input(iexp).catchTrialOutcomeCell = input_temp.catchTrialOutcomeCell;
        mouse(imouse).input(iexp).date = idate;
        mouse(imouse).input(iexp).pctEarly = pctEarly;
        mouse(imouse).input(iexp).pctCorr_maxOri = pctCorr_maxOri;
        mouse(imouse).input(iexp).pctCorr_maxAmp = pctCorr_maxAmp;
        mouse(imouse).input(iexp).reactTimeMs = input_temp.reactTimesMs;
        mouse(imouse).input(iexp).leverUpTimeMs = input_temp.tLeverReleaseTimeMs;
        mouse(imouse).input(iexp).leverDownTimeMs = input_temp.tLeverPressTimeMs;
        mouse(imouse).input(iexp).stimOnTimeMs = input_temp.stimOnTimeMs;
        mouse(imouse).input(iexp).stimOffTimeMs = input_temp.stimOffTimeMs;
        mouse(imouse).input(iexp).catchCyclesOn = input_temp.catchCyclesOn;
        mouse(imouse).input(iexp).tCyclesOn = input_temp.tCyclesOn;
      
        if length(unique(cell2mat_padded(input_temp.tCatchGratingDirectionDeg)))>1
            catch_mat(1,iexp) = 1;
        end
        if length(unique(cell2mat_padded(input_temp.tSoundCatchAmplitude)))>1
            catch_mat(2,iexp) = 1;
        end
        mouse(imouse).input(iexp).catch = catch_mat(:,iexp);
    end
    mouse(imouse).catch_mat = catch_mat;
    mouse(imouse).early_mat = early_mat;
    mouse(imouse).HR_ori_mat = HR_ori_mat;
    mouse(imouse).HR_amp_mat = HR_amp_mat;
    mouse(imouse).date_mat = date_mat;
    
end
save(fullfile(rc.fitOutputSummary, [date '_i613_i614_i616_CatchSummary.mat']), 'mouse');