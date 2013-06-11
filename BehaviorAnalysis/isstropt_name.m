function bool = isstropt_name (optin, names)
%ISSTROPT_NAME (ps-utils): does option name exist in stropt?
%   BOOL = ISSTROPT_NAME (OPTIN, NAMES)
%   OPTIN: stropt (a cell array of name/value pairs, see STROPT_SET)
%   NAMES: cellstr of names to check for existence in OPTIN
%
%   See also STROPT_SET, STROPT_GET, STROPT2STRUCT, STROPT_DEFAULTS
%
%  MH - http://github.com/histed/tools-mh

% chkstropt(optin); done in stropt_names
namelist = stropt_names(optin);
bool = ismember(names, namelist);

