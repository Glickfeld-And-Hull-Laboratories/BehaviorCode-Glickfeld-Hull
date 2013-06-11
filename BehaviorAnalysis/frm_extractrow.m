function xNd = frm_extractrow(xd, desNs)
%FRM_EXTRACTROW Extract rows from a data frame and return a new frame
%
%   xNd = frm_extractrow(xd, desNs)
%
%   
%
% histed 120531

assert(isvector(desNs), 'desNs must be a vector');


fNames = fieldnames(xd);


% this algorithm is simple; if it becomes speed limiting we can optimize it
% (e.g. rewrite in Java)

nDes = length(desNs);
for iD = 1:nDes
    desN = desNs(iD);
    
    for iF = 1:length(fNames)
        tFN = fNames{iF};
        if any(ismember({'colNames', 'nCols', 'nRows' }, tFN))
            continue
        else
            tF = xd.(tFN);
            if iscell(tF)
                if nDes == 1 % special case cell singleton
                    xNd(1).(tFN) = tF{desN};
                else
                    xNd(1).(tFN)(iD) = tF(desN);
                end
                    
            elseif isnumeric(tF)
                xNd(1).(tFN)(iD) = tF(desN);
            else
                error('unknown field type');
            end
        end
    end
end


