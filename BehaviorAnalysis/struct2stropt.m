function stropts = struct2stropt (struct)
%STRUCT2STROPT (ps-utils): convert opt struct to get/set type stropts
%  STROPTS = STRUCT2STROPT(STRUCT)
%
%  MH - http://github.com/histed/tools-mh

s=fieldnames(struct);
v=struct2cell(struct);
stropts=cell(1,length(s).*2);
stropts(1:2:end)=s;
stropts(2:2:end)=v;
