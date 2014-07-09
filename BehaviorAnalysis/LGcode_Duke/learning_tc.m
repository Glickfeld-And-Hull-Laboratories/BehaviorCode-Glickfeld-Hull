datapath = 'Z:\home\andrew\Behavior\Data';
cd(datapath)
mouse_list = [508; 509; 511; 512];
s = struct;
for imouse = 2:length(mouse_list)
list = dir(['data-i' num2str(mouse_list(imouse,:)) '-*']);
siz = size(list);
ndays = 0;
pct_corr = [];
pct_miss = [];
flash = [];
reqhold = [];
min_ori = [];
for ifile = 1:siz(1)
    load(fullfile(datapath,list(ifile).name));
    if isfield(input, 'reactTimeMs')
        ndays = ndays+1;
        A=cell2mat(input.reactTimesMs);
        ind_corr = intersect(find(A>100),find(A<550));
        ind_miss = find(A>550);
        pct_corr = [pct_corr length(ind_corr)./length(A)];
        pct_miss = [pct_miss length(ind_miss)./length(A)];
        if isfield(input, 'stimOnTimeMs')
            flash = [flash 1];
        else
            flash = [flash 0];
        end
        if isfield(input, 'randReqHoldMaxMs')
            reqhold = [reqhold (input.randReqHoldMaxMs + input.fixedReqHoldTimeMs)];
        else
            reqhold = [reqhold ((input.stimOnTimeMs+input.stimOffTimeMs)*input.maxCyclesOn)];
        end
        min_ori = [min_ori min(cell2mat(input.tGratingDirectionDeg),[],2)];
    end
end

figure;
subplot(3,1,1)
plot(pct_corr, 'k');
hold on
plot(pct_miss, 'r');
hold on
plot(1+(flash.*0.1), 'b')
subplot(3,1,2)
plot(reqhold);
subplot(3,1,3)
plot(min_ori);
suptitle(num2str(mouse_list(imouse,:)))

s(imouse).pct_corr = pct_corr;
s(imouse).pct_miss = pct_miss;
s(imouse).flash = flash;
s(imouse).reqhold = reqhold;
s(imouse).min_ori = min_ori;
end