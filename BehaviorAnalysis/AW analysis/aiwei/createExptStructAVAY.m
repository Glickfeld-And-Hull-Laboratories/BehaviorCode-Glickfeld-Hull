function mouse = createExptStructAVAY; %create structure array
rc = behavConstsAV; %the short file, doesn't change, just has AW13 and AW14
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
av = behavParamsAV;%compiling all the data

for imouse = 1:size(av,2);
    mouse_name = av(imouse).mouse;
    ind = find(xd.Subject == mouse_name);
    date_mat = [];
    early_mat = [];
    HR_ori_mat = [];
    HR_amp_mat = [];
    catch_mat = zeros(2, length(ind));%putting together the data
    for iexp = 1:length(ind)
        idate = xd.DateStr{ind(iexp)};
        n  = dir(fullfile(rc.pathStr, ['data-i' num2str(mouse_name) '-' idate '-*']));%directory of file
        if ~isnan(str2num(xd.ChooseMatFile{ind(iexp)}))
            n = n(str2num(xd.ChooseMatFile{ind(iexp)}));
        end
        for ifile = 1:size(n,1)
            if ifile == 1
                input_temp = mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], []);
            else
                input_temp = [input_temp mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], [])];
            end
        end
        if size(n,1)>1
            input_temp = concatenateDataBlocks(input_temp);%concatenating
        end
        range = str2num(xd.TrialRangeToUse{ind(iexp)});%figuring out the range of trials inputed 
        if ~isnan(range(1))
            input_temp = trialChopper(input_temp, [range(1) range(end)]);
        end
       date_mat = [date_mat; idate];
        
        failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
        missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
        successIx = strcmp(input_temp.trialOutcomeCell, 'success');
        pctEarly = sum(failureIx,2)./length(failureIx); %percentage earlys
        early_mat = [early_mat; pctEarly];
        
        gratingDirectionDeg = cell2mat_padded(input_temp.tGratingDirectionDeg);
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);
        
        reactTimes = celleqel2mat_padded(input.reactTimesMs);
        
        oris = unique(gratingDirectionDeg);
        amps = unique(soundAmplitude);
        maxOriTrials = find(gratingDirectionDeg == max(oris,[],1));
        maxAmpTrials = find(soundAmplitude == max(amps,[],2));
        pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
        pctCorr_maxAmp = sum(successIx(maxAmpTrials),2)./(sum(successIx(maxAmpTrials),2)+sum(missedIx(maxAmpTrials),2));
       
        HR_ori_mat = [HR_ori_mat; pctCorr_maxOri];%hit rates
        HR_amp_mat = [HR_amp_mat; pctCorr_maxAmp];
        
        if ~isfield(input_temp, 'catchTrialOutcomeCell')%true if field is in a structure array 
            for trN = 1:length(input_temp.trialOutcomeCell)%identify what trial is what
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
        mouse(imouse).input(iexp).reactTimesMs = input_temp.reactTimesMs; 
        mouse(imouse).input(iexp).date = idate;
        mouse(imouse).input(iexp).pctEarly = pctEarly;
        mouse(imouse).input(iexp).pctCorr_maxOri = pctCorr_maxOri;
        mouse(imouse).input(iexp).pctCorr_maxAmp = pctCorr_maxAmp;
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
save(fullfile(rc.fitOutputSummary, [date '_i613_i614_CatchSummary.mat']), 'mouse');