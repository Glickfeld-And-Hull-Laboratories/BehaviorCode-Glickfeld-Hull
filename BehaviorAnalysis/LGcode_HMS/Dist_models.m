mu_sig = [0.1:.05:.45];
std_sig = [0.1*ones(1,length(mu_sig))];
std_noise = [.1];
mu_noise = [0];

offset_shift = [0 .03];
r_noise = zeros(length(mu_noise),length(offset_shift),1000000);    
for ishift = 1:length(offset_shift)
    for idist = 1:length(std_noise)
        r_noise(idist,ishift,:) = normrnd(mu_noise(idist)-offset_shift(ishift),std_noise(idist), [1 1000000]);
    end
end
r_sig = zeros(length(mu_sig),length(offset_shift),1000000);
for ishift = 1:length(offset_shift)
    for idist = 1:length(mu_sig)
        r_sig(idist,ishift,:) = normrnd(mu_sig(idist)-offset_shift(ishift),std_sig(idist), [1 1000000]);
    end
end

col_mat = strvcat('r', 'm');
siz = ceil(sqrt(length(mu_sig)));
for ishift = 1:(length(offset_shift))
    figure;
    for idist_sig = 1:(length(mu_sig))
        subplot(siz,siz,idist_sig)
        hist(squeeze(r_noise(:,ishift,:)))
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',col_mat(ishift,:))
        hold on
        hist(squeeze(r_sig(idist_sig,ishift,:)))
        xlim([-0.5 1])
        title(['Mu = ' num2str(mu_sig(idist_sig))]);
    end
    if ishift == 1
        suptitle(['Red: Offset shift = -' num2str(offset_shift(ishift))])
    elseif ishift == 2
        suptitle(['Magenta: Offset shift = -' num2str(offset_shift(ishift))])
    end
    pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['histograms_sig_offset_shift_pt0' num2str(offset_shift(ishift)*100) '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end

crit = [0.1:0.05:0.35];
siz = ceil(sqrt(length(crit)));
hit = zeros(length(mu_sig), length(crit), length(offset_shift));
FAs = zeros(length(crit),length(offset_shift));
thresh = zeros(length(crit),length(offset_shift));
d = zeros(length(mu_sig),length(crit),length(offset_shift));
abs_crit = zeros(length(mu_sig),length(crit),length(offset_shift));
figure;
for ishift = 1:length(offset_shift)
    for icrit = 1:length(crit)   
        FAs(icrit, ishift) = length(find(r_noise(:,ishift,:)>crit(icrit)))./size(r_noise,3);
        for idist_sig = 1:length(mu_sig)
            hit(idist_sig,icrit,ishift) = length(find(r_sig(idist_sig,ishift,:)>crit(icrit)))./size(r_sig,3);
            [d(idist_sig,icrit,ishift) abs_crit(idist_sig,icrit,ishift)] = dprime(hit(idist_sig,icrit,ishift), FAs(icrit,ishift));
        end
        subplot(subplot(siz,siz,icrit))
        plot(mu_sig', squeeze(hit(:,icrit,ishift)), ['o' col_mat(ishift,:)])
        hold on
        ylim([0 1])
        xlim([0 .5])
        [y x] = sort(abs(hit(:,icrit,ishift)-.5));
        m = y(1)+y(2)./(mu_sig(x(1))-mu_sig(x(2)));
        b = y(1)-(m.*mu_sig(x(1)));
        thresh(icrit,ishift) = -(b/m);
        title(['criteria = ' num2str(crit(icrit))]);
    end
end
suptitle(['Red: Offset shift = -' num2str(offset_shift(1)) '  Magenta: Offset shift = -' num2str(offset_shift(2))])
pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['threshold-by-crit-sig-offset-shift.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

figure;
for ishift = 1:length(offset_shift)
    subplot(2,2,1)
    scatter(crit, FAs(:,ishift), col_mat(ishift,:))
    hold on
    xlabel('Criterion (% contrast)')
    ylabel('False alarm rate')
    xlim([0 0.4])
    subplot(2,2,2)
    scatter(crit, thresh(:,ishift),col_mat(ishift,:))
    hold on
    xlabel('Criterion (% contrast)')
    ylabel('Threshold')
    xlim([0 0.4])
    ylim([0 0.5])
    subplot(2,2,3)
    for icrit = 1:length(crit)
        diffs = abs(bsxfun(@minus, mu_sig, thresh(icrit,ishift)));
        [val con] = min(diffs);
        scatter(crit(:,icrit), d(con,icrit,ishift),col_mat(ishift,:))
        hold on
    end
    xlabel('Criterion (% contrast)')
    ylabel('dprime at threshold')
    xlim([0 0.4])
    ylim([0 4])
    subplot(2,2,4)
    con = 2;
    for icrit = 1:length(crit)
        scatter(crit(:,icrit), d(con,icrit,ishift),col_mat(ishift,:))
        hold on
    end
    xlabel('Criterion (% contrast)')
    ylabel('dprime at 15% contrast')
    xlim([0 0.4])
    ylim([0 4])
end
suptitle(['Red: Offset shift = -' num2str(offset_shift(1)) '  Magenta: Offset shift = -' num2str(offset_shift(2))])
pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['criteria-noise-plots-sig-offset-shift.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

gain_shift = [1 .7];
r_noise = zeros(length(mu_noise),length(gain_shift),1000000);
for ishift = 1:length(gain_shift)
    for idist = 1:length(std_noise)
        r_noise(idist,ishift,:) = normrnd(mu_noise(idist)*gain_shift(ishift),std_noise(idist), [1 1000000]);
    end
end
r_sig = zeros(length(mu_sig),length(gain_shift),1000000);
for ishift = 1:length(gain_shift)
    for idist = 1:length(mu_sig)
        r_sig(idist,ishift,:) = normrnd(mu_sig(idist)*gain_shift(ishift),std_sig(idist), [1 1000000]);
    end
end

col_mat = strvcat('r', 'm');
siz = ceil(sqrt(length(mu_sig)));
for ishift = 1:(length(gain_shift))
    figure;
    for idist_sig = 1:(length(mu_sig))
        subplot(siz,siz,idist_sig)
        hist(squeeze(r_noise(:,ishift,:)))
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',col_mat(ishift,:))
        hold on
        hist(squeeze(r_sig(idist_sig,ishift,:)))
        xlim([-0.5 1])
        title(['Mu = ' num2str(mu_sig(idist_sig))]);
    end
    if ishift == 1
        suptitle(['Red: Gain shift = x' num2str(gain_shift(ishift))])
    elseif ishift == 2
        suptitle(['Magenta: Gain shift = x' num2str(gain_shift(ishift))])
    end
    pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['histograms_sig_gain_shift_pt' num2str(gain_shift(ishift)) '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end

crit = [0.1:0.05:0.35];
siz = ceil(sqrt(length(crit)));
hit = zeros(length(mu_sig), length(crit), length(gain_shift));
FAs = zeros(length(crit),length(gain_shift));
thresh = zeros(length(crit),length(gain_shift));
d = zeros(length(mu_sig),length(crit),length(gain_shift));
abs_crit = zeros(length(mu_sig),length(crit),length(gain_shift));
figure;
for ishift = 1:length(gain_shift)
    for icrit = 1:length(crit)   
        FAs(icrit, ishift) = length(find(r_noise(:,ishift,:)>crit(icrit)))./size(r_noise,3);
        for idist_sig = 1:length(mu_sig)
            hit(idist_sig,icrit,ishift) = length(find(r_sig(idist_sig,ishift,:)>crit(icrit)))./size(r_sig,3);
            [d(idist_sig,icrit,ishift) abs_crit(idist_sig,icrit,ishift)] = dprime(hit(idist_sig,icrit,ishift), FAs(icrit,ishift));
        end
        subplot(subplot(siz,siz,icrit))
        plot(mu_sig', squeeze(hit(:,icrit,ishift)), ['o' col_mat(ishift,:)])
        hold on
        ylim([0 1])
        xlim([0 .5])
        [y x] = sort(abs(hit(:,icrit,ishift)-.5));
        m = y(1)+y(2)./(mu_sig(x(1))-mu_sig(x(2)));
        b = y(1)-(m.*mu_sig(x(1)));
        thresh(icrit,ishift) = -(b/m);
        title(['criteria = ' num2str(crit(icrit))]);
    end
end
suptitle(['Red: Gain shift = x' num2str(gain_shift(1)) '  Magenta: Gain shift = x' num2str(gain_shift(2))])
pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['threshold-by-crit-sig-gain-shift.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

figure;
for ishift = 1:length(gain_shift)
    subplot(2,2,1)
    scatter(crit, FAs(:,ishift), col_mat(ishift,:))
    hold on
    xlabel('Criterion (% contrast)')
    ylabel('False alarm rate')
    xlim([0 0.4])
    subplot(2,2,2)
    scatter(crit, thresh(:,ishift),col_mat(ishift,:))
    hold on
    xlabel('Criterion (% contrast)')
    ylabel('Threshold')
    xlim([0 0.4])
    ylim([0 0.5])
    subplot(2,2,3)
    for icrit = 1:length(crit)
        diffs = abs(bsxfun(@minus, mu_sig, thresh(icrit,ishift)));
        [val con] = min(diffs);
        scatter(crit(:,icrit), d(con,icrit,ishift),col_mat(ishift,:))
        hold on
    end
    xlabel('Criterion (% contrast)')
    ylabel('dprime at threshold')
    xlim([0 0.4])
    ylim([0 4])
    subplot(2,2,4)
    con = 2;
    for icrit = 1:length(crit)
        scatter(crit(:,icrit), d(con,icrit,ishift),col_mat(ishift,:))
        hold on
    end
    xlabel('Criterion (% contrast)')
    ylabel('dprime at 15% contrast')
    xlim([0 0.4])
    ylim([0 4])
end
suptitle(['Red: Gain shift = x' num2str(gain_shift(1)) '  Magenta: Gain shift = x' num2str(gain_shift(2))])
pn = fullfile('~/Documents/Figures/Signal_noise_simulations', ['criteria-noise-plots-sig-gain-shift.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
    
        