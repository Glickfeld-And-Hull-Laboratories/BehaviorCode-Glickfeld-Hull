function textH = suptitle2(str) 
%SUPTITLE2 (ps-utils): title over all subplots; improved version (2)
%   SUPTITLE2('text') adds text to the top of the figure
%   above all subplots (a "super title"). 
% 
%   This is much simpler than the old SUPTITLE because it uses the new
%   OuterPosition property.  
%   It must force a redraw for this property to be computed correctly.
%   (as of Matlab 2011a)
%
%   You should call this after all subplots are drawn, for speed reasons, but
%   it should resize all your subplots correctly even if used between plots.
%
%   MH - http://github.com/histed/tools-mh

suptitleNormHeight = 0.04;

% save old axes
savedEntryFigH = gcf;
savedEntryAxesH = gca;

% look for old
supH = findobj(gcf, 'Type', 'axes', '-and', 'Tag', 'suptitle2');

if ~isempty(supH)
    if length(supH) > 1
        error('Bug in suptitle2: more than one invis axes?');
    end
    
    % restore old positions
    ud = get(supH, 'UserData');
    validHList = ud.oldSubplotHandles(ishandle(ud.oldSubplotHandles));
    if ~isempty(validHList)
        set(validHList, {'OuterPosition'}, ud.oldOuterPositions);
    end
    

    % remove so we can make a new one below
    delete(supH);
end


% find existing subplots
allAxH = findobj(gcf, 'Type', 'axes');
otherH = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
subH = setdiff(allAxH, otherH);


if any(strcmp('Colorbar', get(allAxH, 'Tag')))
    % this function needs to be updated to handle colorbars.  The problem is
    % that matlab's automatic colorbar positioning code is very complicated, 
    % resizes the axis, and it's hard to guess what it's going to do.  It
    % might be possible to call colorbar in some way to get it to update the
    % axis position-- but I couldn't figure out an easy way to do this.
    error('Draw colorbars after calling suptitle2');
end



%% goal here: scale y height and y position to shrink all axes.  To scale the
% whole set of plots correctly we must use the y end position not the start
% (i.e. pos(:,2)+pos(:,4) not pos(:,2)

% bug in matlab 2011a: does not compute outerposition correctly until after a
% drawnow, when(I think) activepositionproperty is 'position'.  It's possible
% to try computing from Position and TightInset but the most robust solution
% is to draw.
drawnow;

opC = get(subH, 'OuterPosition');
if ~iscell(opC), opC = {opC}; end  % one subplot
opM = cat(1, opC{:});
yEndOld = opM(:,2)+opM(:,4);
yHtNew = opM(:,4)*(1-suptitleNormHeight);
yEndNew = yEndOld*(1-suptitleNormHeight);
opNew = [opM(:,1) (yEndNew-yHtNew) opM(:,3) yHtNew];
newC = mat2cell(opNew, ones([1,length(subH)]), 4);
set(subH, {'OuterPosition'}, newC);


invisH = axes('Units', 'normalized', ...
              'Position', [0 1-suptitleNormHeight 1 suptitleNormHeight], ...
              'Visible', 'off', ...
              'Tag', 'suptitle2');
textH = text(0.5, 0.35, str, ...
             'Units', 'normalized', ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Tag', 'suptitle2');

% save old positions in axes for later update
ud.oldOuterPositions = opC;
ud.oldSubplotHandles = subH;
set(invisH, 'UserData', ud);

% restore old axis, without using a drawnow
set(0, 'CurrentFigure', savedEntryFigH);
set(gcf, 'CurrentAxes', savedEntryAxesH);



