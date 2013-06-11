function[area]=roc_forloop(pspikes,nspikes)
%pspikes and nspikes are row vectors of spike counts for preferred and
%null choices or stimuli
%
% from MRC 090519
%  MH - http://github.com/histed/tools-mh
    
cmax=150;
rip=zeros(cmax+1,1);
rin=zeros(cmax+1,1);
total_p=length(pspikes);
total_n=length(nspikes);
if total_p>3 && total_n>3
    for criteria=0:cmax
        rip(criteria+1)=(sum(pspikes>=criteria))/total_p;
        rin(criteria+1)=(sum(nspikes>=criteria))/total_n;
    end
    rip=flipud(rip);
    rin=flipud(rin);
    temp=diff(rin);
    temp2=(rip(2:end)+rip(1:end-1))/2;
    area=sum(temp2.*temp);
else area=NaN;
end

%% debug
%fprintf(1, 'data: sum %4d %4d, length %4d %4d\n', ...
%        sum(pspikes), sum(nspikes), length(pspikes), length(nspikes));
