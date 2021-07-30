function [ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(ang1,sp1,ang2,sp2);

x0 = 0;
y0 = 0;
[x1 y1] = pol2cart(ang1,sp1);
[x2 y2] = pol2cart(ang2,sp2);
ang3 = (ang1+ang2)./2;
sp3 = (sp1+sp2)/2;
[x3 y3] = pol2cart(ang3,sp3);

figure;
plot([x0 x1],[y0 y1])
hold on
plot([x0 x2],[y0 y2])

x3 = (x1+x2)/2;
y3 = (y1+y2)/2;
plot([x0 x3],[y0 y3])
[av_ang,av_sp] = cart2pol(x3,y3);


if rem(ang1,pi/2) ~= 0
    x1p = -40:40;
    s1 = y1/x1;
    s1p = -1/s1;
    b1 = y1-(s1p.*x1);
    y1p = (s1p.*x1p) + b1;
elseif rem(ang1,pi) == 0
    x1p = sp1.*ones(size([-40:40]));
    y1p = [-40:40];
    s1p = y1p(end)/x1p(end);
    b1 = y1p(end)-(s1p.*x1p(end));
elseif rem(ang1,pi/2) == 0
    y1p = sp1.*ones(size([-40:40]));
    x1p = [-40:40];
    s1p = y1p(end)/x1p(end);
    b1 = y1p(end)-(s1p.*x1p(end));
end
plot(x1p,y1p,'--k')


if rem(ang2,pi/2) ~= 0
    x2p = -40:40;
    s2 = y2/x2;
    s2p = -1/s2;
    b2 = y2-(s2p.*x2);
    y2p = (s2p.*x2p) + b2;
elseif rem(ang2,pi) == 0
    x2p = sp2.*ones(size([-40:40]));
    y2p = [-40:40];
    s2p = y2p(end)/x2p(end);
    b2 = y2p(end)-(s2p.*x2p(end));
elseif rem(ang2,pi/2) == 0
    y2p = sp2.*ones(size([-40:40]));
    x2p = [-40:40];
    s2p = y2p(end)/x2p(end);
    b2 = y2p(end)-(s2p.*x2p(end));
end
plot(x2p,y2p,'--k')

if rem(ang1,pi/2) ~= 0 & rem(ang2,pi/2) ~= 0
    xint = (b2-b1)./(s1p-s2p);
    yint = (s1p.*xint) + b1;
elseif rem(ang1,pi) == 0 || rem(ang2,pi)  == 0
    if rem(ang1,pi) == 0
        xint = sp1;
        yint = (s2p.*xint) + b2;
    else
        xint = sp2;
        yint = (s1p.*xint) + b1;
    end
elseif rem(ang1,pi/2) == 0 || rem(ang2,pi/2)  == 0
    if rem(ang1,pi/2) == 0
        yint = sp1;
        xint = (yint-b2)./s2p;
    else
        yint = sp2;
        xint = (yint-b1)./s1p;
    end
end
plot([x0 xint],[y0 yint])
[ioc_ang,ioc_sp] = cart2pol(xint,yint);
xlim([-40 40])
ylim([-40 40])
hline(0)
vline(0)

end


