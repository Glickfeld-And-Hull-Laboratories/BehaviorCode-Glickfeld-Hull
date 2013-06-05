from pylab import *

xvals = [5, 10,20,30,50]
yvals = [300,425,530,700,950]

figure
plot(xvals, yvals)

xlabel('time open (ms)')
ylabel('uL per 500 corrects')


#canvas = gcf().canvas
#canvas.start_event_loop(timeout=10)

