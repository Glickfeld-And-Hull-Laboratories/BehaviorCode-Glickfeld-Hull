function [FAIx CRIx] = findCatchTrialOutcomeFSAV(catchGratingDirection,isFalseAlarm)

catchTrials_ind = find(cell2mat_padded(catchGratingDirection) > 0 );

catchFA = cell2mat(isFalseAlarm(catchTrials_ind));

FAIx = catchTrials_ind(find(catchFA > 0));
CRIx = catchTrials_ind(find(catchFA == 0));

end