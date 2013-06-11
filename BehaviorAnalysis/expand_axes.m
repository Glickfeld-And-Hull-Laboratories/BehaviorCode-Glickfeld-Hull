function expand_axes(axH);
%
%   Make axes as big as possible within the OuterPosition.  Note that you can 
%
% histed 111016

if nargin < 1
    axH = gca;
end


op = get(axH, 'OuterPosition');
p = get(axH, 'Position');
newP = op;
oSz = op(3:4);

%xMargins = [0.2 0.2];
%yMargins = [0.2 0.2]; 
xMargins = [0.06 0.01];
yMargins = [0.085 0.05]; 

xFSize = 1 - sum(xMargins);
yFSize = 1 - sum(yMargins);


newP(1:2) = op(1:2) + oSz.*[xMargins(1) yMargins(1)];
newP(3:4) = oSz.*[xFSize yFSize];
%ti = get(axH, 'TightInset')
%ewP = op + [ti(1:2), -(ti(1:2)+ti(3:4))]

set(gca, 'Position', newP);
%set(gca, 'ActivePositionProperty', 'outerposition');

fprintf(1, 'new, old: [%4.2f %4.2f, %4.2f %4.2f, %4.2f %4.2f, %4.2f %4.2f]\n', ...
        chop([newP(1) op(1) newP(2) op(2), ...
              newP(3) op(4) newP(3) op(4)],2));
