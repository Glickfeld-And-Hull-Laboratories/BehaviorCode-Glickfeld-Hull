%% code for making behavior datasets with different strategies
%% params
t_spo = 2;
t_max = 1;
d_max = 100;
d_spo = 0.37; %i459's SPO
n = 1:4;

%calculate all targets and differences
t =  fliplr(t_max./(2.^((n-1)./(t_spo))));
d =  fliplr(d_max./(2.^((n-1)./(d_spo)))+1);

%random trial generator
s = rng('default');

ntrials = 2516;
x_t = randi(n(end),[1 ntrials]); %which target [1:4]
x_d = randi(n(end),[1 ntrials]); %which diff [1:4]
x_s = randi([0 1],[1 ntrials]); %target location: 0-left 1-right

x_tcon = t(x_t); %actual target con
x_dcon = x_tcon./d(x_d); %actual distractor con
r_con = x_tcon;
r_con(find(x_s==0)) = x_dcon(find(x_s==0)); %actual right con
l_con = x_tcon;
l_con(find(x_s==1)) = x_dcon(find(x_s==1)); %actual left con

%% contrast discrim
%probability of correct choice by contrast diff
% d_p = [0.6 0.7 0.85 1]; % based on arbitrary input
d_p = [0.4777 0.8024 0.9733 0.9711]; % based on i459's actual prob right
y_s = zeros(1,ntrials);
for i = 1:ntrials
    if rand(1)<=d_p(x_d(i))
        y_s(i) = x_s(i);
    else
        y_s(i) = abs(x_s(i)-1);
    end
end

%% right side con
%proportional to right side contrast
y_s = zeros(1,ntrials);
for i = 1:ntrials
    if rand(1)<r_con(i) 
        y_s(i) = 1;
    else
        y_s(i) = 0;
    end
end

%% proportional to right side contrast based on fract targets
y_fit = fit(r_con',x_s','poly2');
x = 0:0.001:1;
p = y_fit.p1.*(x.^2)+ y_fit.p2.*x+y_fit.p3;
figure; scatter(r_con, x_s,'o');
hold on, plot(x,p)

y_s = zeros(1,ntrials);
for i = 1:ntrials
    p = y_fit.p1.*(r_con(i).^2)+ y_fit.p2.*r_con(i)+y_fit.p3;
    if rand(1)<p 
        y_s(i) = 1;
    else
        y_s(i) = 0;
    end
end


%% plot
x_d_sign = x_d;
x_d_sign(find(x_s==0)) = x_d(find(x_s==0)).*-1;
d_sign = unique(x_d_sign);
n_d = length(d_sign);
y_pctright = zeros(1,n_d);
y_CIright = zeros(2,n_d);
for i = 1:n_d
    ind = find(x_d_sign == d_sign(i));
    [y_pctright(:,i), y_CIright(:,i)] = binofit(sum(y_s(ind)),length(ind));
end
figure; errorbar(d_sign, y_pctright, y_pctright-y_CIright(1,:), y_CIright(2,:)-y_pctright);




