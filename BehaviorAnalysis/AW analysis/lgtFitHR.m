function [HRfit,xFit,stats] = lgtFitHR(x_all,tr_ind,hits,misses)
    ind = tr_ind & (hits | misses);
    xFit = x_all(ind);
    y = hits(ind);
    
    [b,~,stats] = glmfit(xFit,y,'binomial');
    HRfit = glmval(b,xFit,'logit');
end