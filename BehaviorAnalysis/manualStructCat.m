function temp = manualStructCat(temp, input);

% adds field from one struct (input) as new dimension in another struct (temp). 

    nrun = size(temp,2) + 1;

    fNamesA = fieldnames(input);
    nFA = length(fNamesA);

    for iF = 1:nFA
        a = eval(cell2mat(['input.' fNamesA(iF)]));
        eval(cell2mat(['temp(nrun).' fNamesA(iF) ' = a']));
    end
end