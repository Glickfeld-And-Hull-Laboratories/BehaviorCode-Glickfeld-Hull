function [ ds ] = dataComparator(a, b)
%DATACOMPARATOR: reads all the fields in a, compares those values to b,
%looks for differences and quantifies those differences. A will always be
%compared against B in a A-B fashion.
%%
ds = struct;
ds.a = a;
ds.b = b;
fieldNames = fieldnames(a);
nFNs = length(fieldNames);
nTrs = length(a.trialOutcomeCell);

for i=1:nFNs,
    c = strcat('input.', fieldNames(i));
    field = eval(c{1});
    if isfield(b, fieldNames{i})
        disp('Cool beans')
    else
        disp('Sour grapes')
    end
end
end
