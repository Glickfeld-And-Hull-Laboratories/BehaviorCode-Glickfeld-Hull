function cat_struct = CatStructFields(dim,varargin)
fields = cellfun(@fieldnames,varargin,'un',0);
fields = fields{1};
% fields{1} = ('trialOutcomeCell');
% fields{2} = ('tLeftTrial');

T = [varargin{:}];
cat_struct = struct();

for k = 1:numel(fields)
  aField     = char(fields(k));
  cat_struct.(aField) = cat(dim, T.(aField));
end