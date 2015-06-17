function mouse = createExptStruct2AFC
%takes all data from xls, sorts according to task type and gets thresholds
%etc
rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
twoafc = behavParams2AFC;
for imouse = 1:length(twoafc.mouse_mat)
    mouse(imouse).ind = find(xd.Subject ==  twoafc.mouse_mat(imouse));
    mouse(imouse).name = twoafc.mouse_mat(imouse);
    for ipow = 1:length(twoafc.power_mat)
        pow = twoafc.power_mat(ipow);
        mouse(imouse).pow(ipow).ind = intersect(find(xd.Power==pow), mouse(imouse).ind);
        mouse(imouse).pow(ipow).name = pow;
    end
end

for imouse = 1:length(twoafc.mouse_mat)
    mouseStr = num2str(twoafc.mouse_mat(imouse));
    for ipow = 1:length(twoafc.power_mat)
        mouse(imouse).pow(ipow).expt = [];
        nexp = length(mouse(imouse).pow(ipow).ind);
        thresh = zeros(nexp,2);
        for iexp = 1:nexp
            dateStr = cell2mat(xd.DateStr(mouse(imouse).pow(ipow).ind(iexp)));
            disp([dateStr ' i' mouseStr])
            fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(twoafc.mouse_mat(imouse)) '-' dateStr '-' cell2mat(xd.DataBlock(mouse(imouse).pow(ipow).ind(iexp))) '.mat']);
            load(fn);
            mouse(imouse).pow(ipow).expt(iexp).thresh = [fitS(1).thresh fitS(2).thresh];
            mouse(imouse).pow(ipow).expt(iexp).slope = [fitS(1).slope fitS(2).slope];
            mouse(imouse).pow(ipow).expt(iexp).ci95 = [fitS(1).bootStats.ci95; fitS(2).bootStats.ci95];
            subplot(2,2,ipow)
            plot(1:2,mouse(imouse).pow(ipow).expt(iexp).thresh,'-k')
            hold on
            title([num2str(mouse(imouse).pow(ipow).name) ' mW'])
        end
        ylim([0 1.5])
    end
    suptitle(['i' mouseStr])
    fname = ['ThresholdSummary_i' num2str(twoafc.mouse_mat(imouse)) '_2AFCtask.pdf'];
    fn_out = fullfile(rc.fitOutputSummary, fname);
    exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
end

x = date;
fn_out = fullfile(rc.fitOutputSummary, ['2AFC_summarystruct.mat']);
save(fn_out, 'mouse');