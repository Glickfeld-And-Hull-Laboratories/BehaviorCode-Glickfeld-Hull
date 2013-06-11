function cellFieldValues = extract_fields(inStruct, desiredFields)
%EXTRACT_FIELDS (ps-utils): retun vals of specified struct fields in cell array
%  CELLFIELDVALUES = EXTRACT_FIELDS(STRUCTURE, DESIREDFIELDS)
%  cellFieldValues is a cell vector of same length as desiredFields.  The value
%  of each field in field names is put into one cell of this cell vector.
%
%  MH - http://github.com/histed/tools-mh

allNames = fieldnames(inStruct);
allVals = struct2cell(inStruct);

desNameIx = ismember(allNames, desiredFields);
cellFieldValues = allVals(desNameIx);
