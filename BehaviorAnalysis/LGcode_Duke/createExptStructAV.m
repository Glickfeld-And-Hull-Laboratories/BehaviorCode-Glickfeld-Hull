function mouse = createExptStructAV;
rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
av = behavParamsAV;
mice = unique(xd.Subject);
nMice = length(mice);

for imouse = 1:nMice;
    mouse_name = mice(imouse);
    ind = find(xd.Subject == mouse_name);
    date_mat = [];
    early_mat = [];
    HR_ori_mat = [];
    HR_amp_mat = [];
    catch_mat = zeros(2, length(ind));
   
    for iexp = 1:length(ind)
        iexp
        idate = xd.DateStr{ind(iexp)};
        n  = dir(fullfile(rc.pathStr, ['data-i' num2str(mouse_name) '-' idate '-*']));
%         if ~isnan(str2num(xd.ChooseMatFile{ind(iexp)}))
%             n = n(str2num(xd.ChooseMatFile{ind(iexp)}));
%         end
        for ifile = 1:size(n,1)
            if ifile == 1
                input_temp = mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], []);
            else
                input_temp = [input_temp mwLoadData(fullfile(rc.pathStr, n(ifile).name), [], [])];
            end
        end
        if size(n,1)>1
            input_temp = concatenateDataBlocks(input_temp);
        end
        range = str2num(xd.TrialRangeToUse{ind(iexp)});
        if ~isnan(range(1))
            if range(end) > input_temp.trialSinceReset
                input_temp = trialChopper(input_temp, [range(1) input_temp.trialSinceReset]);
            else
                input_temp = trialChopper(input_temp, [range(1) range(end)]);
            end
        end
        date_mat = [date_mat; idate];
        
        failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
        missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
        successIx = strcmp(input_temp.trialOutcomeCell, 'success');
        pctEarly = sum(failureIx,2)./length(failureIx);
        early_mat = [early_mat; pctEarly];
        
        gratingDirectionDeg = celleqel2mat_padded(input_temp.tGratingDirectionDeg);
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);
        
%         catchgratingdirectiondeg = celleqel2mat_padded(input_temp.tCatchGratingDirectionDeg);
       
        

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
        mouse(imouse).input(iexp).tCyclesOn = input_temp.tCyclesOn
        mouse(imouse).input(iexp).tCatchTimeMs = input_temp.tCatchTimeMs;
        mouse(imouse).input(iexp).uniqueCatchDeg = unique(celleqel2mat_padded(input_temp.tCatchGratingDirectionDeg));
        mouse(imouse).input(iexp).uniqueDeg= unique(mouse(imouse).input(iexp).uniqueCatchDeg(isnan(mouse(imouse).input(iexp).uniqueCatchDeg) ==0));
        mouse(imouse).input(iexp).FAIx =strcmp(input_temp.catchTrialOutcomeCell, 'FA');
        mouse(imouse).input(iexp).CRIx =strcmp(input_temp.catchTrialOutcomeCell, 'CR');
        mouse(imouse).input(iexp).Ignorex =strcmp(input_temp.trialOutcomeCell, 'ignore');
        mouse(imouse).input(iexp).SuccessIx =strcmp(input_temp.trialOutcomeCell, 'success');
        mouse(imouse).input(iexp).hitrates = zeros(1,length(mouse(imouse).input(iexp).uniqueDeg));
        
        
        mouse(imouse).input(iexp).countDegs = length(mouse(imouse).input(iexp).uniqueDeg);
        
        for level = 1: length(mouse(imouse).input(iexp).uniqueDeg)
            
            deg = mouse(imouse).input(iexp).uniqueDeg(level);
            mouse(imouse).input(iexp).find = find(celleqel2mat_padded(mouse(imouse).input(iexp).tCatchGratingDirectionDeg)==deg);
            mouse(imouse).input(iexp).FAs =  mouse(imouse).input(iexp).FAIx( mouse(imouse).input(iexp).find)
            mouse(imouse).input(iexp).levelFA = sum( mouse(imouse).input(iexp).FAs)/ length(mouse(imouse).input(iexp).find)
            mouse(imouse).input(iexp).hitrates(level) = mouse(imouse).input(iexp).levelFA;
            
        end
        
          

        
%         mouse(imouse).input(iexp).levelone= find(celleqel2mat_padded(input_temp.tCatchGratingDirectionDeg)==22.5000);
%         mouse(imouse).input(iexp).leveltwo= find(celleqel2mat_padded(input_temp.tCatchGratingDirectionDeg)==90);
%         
% 
%         
%         mouse(imouse).input(iexp).oneFA =  mouse(imouse).input(iexp).FAIx( mouse(imouse).input(iexp).levelone);
%         mouse(imouse).input(iexp).twoFA= mouse(imouse).input(iexp).FAIx(mouse(imouse).input(iexp).leveltwo);
%         mouse(imouse).input(iexp).leveloneFA = sum( mouse(imouse).input(iexp).oneFA)/ length(mouse(imouse).input(iexp).levelone);
%         mouse(imouse).input(iexp).leveltwoFA = sum( mouse(imouse).input(iexp).twoFA)/ length(mouse(imouse).input(iexp).leveltwo);
      
        totalCatchLevels = double(input_temp.catchTrPerB2Level1+input_temp.catchTrPerB2Level2 + input_temp.catchTrPerB2Level3 +input_temp.catchTrPerB2Level4+input_temp.catchTrPerB2Level5+input_temp.catchTrPerB2Level6+input_temp.catchTrPerB2Level7+input_temp.catchTrPerB2Level8)./80;
        mouse(imouse).input(iexp).catchvalue = totalCatchLevels;
        
        mouse(imouse).input(iexp).measuredcatchvalue = (length(input_temp.catchTrialOutcomeCell) - sum(strcmp(input_temp.catchTrialOutcomeCell, 'NaN')))/ length(input_temp.catchTrialOutcomeCell) ;
        
        
        mouse(imouse).input(iexp).measured = (mouse(imouse).input(iexp).FAIx+  mouse(imouse).input(iexp).CRIx)/ (mouse(imouse).input(iexp).FAIx+  mouse(imouse).input(iexp).CRIx+ mouse(imouse).input(iexp).Ignorex+ mouse(imouse).input(iexp).SuccessIx);
        
        mouse(imouse).input(iexp).measuredDegs = ones(1,mouse(imouse).input(iexp).countDegs) *mouse(imouse).input(iexp).measured;
        
        
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
save(fullfile(rc.fitOutputSummary, [date '_CatchSummary.mat']), 'mouse');
