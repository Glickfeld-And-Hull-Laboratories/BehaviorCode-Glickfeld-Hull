function set_struct_in_map(inMap, key, varargin);
%SET_STRUCT_IN_MAP (mh): set fields of structures that are the values of a Containers.Map
%
%  outMap = set_struct_in_map(inMap, key, 'structFieldName1', value1, ...);
%    key is the map key for the struct
%    
%  Note!  No value is returned because as a Java object, MATLAB's Containers.Map is 
%  altered in place (passed by reference, not by value like all other matlab objects).  
%   
% histed 130124

if nargout > 0, error('Containers.Map is altered *in-place*!  (It''s a Java object)'); end

temp0 = inMap(key);
assert(mod(length(varargin)/2,1)==0, 'must provide pairs');
nToSet = length(varargin)/2;
for iN = 1:nToSet;
    fieldN = (nToSet-1)*2+1;
        tFieldName = varargin{fieldN};
    tVal = varargin{fieldN+1};
    temp0.(tFieldName) = tVal;
end
inMap(key) = temp0;


