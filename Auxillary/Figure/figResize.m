function [ axi ] = figResize( width, height, pad, subpad, fontsize, mm )
%RESIZEFIGURETIGHT - resize figure on screen


% Padding incase we need extra space around the figure
if nargin < 3
    padding = [20,40,20,20];
elseif numel(pad) == 1
    padding = [pad,pad,pad,pad];
else
    padding = pad;
end

if nargin < 4
    subpad = 10;
end

if nargin < 5
    fontsize = get(gca,'FontSize');
end

if nargin < 6
    mm = 0;
end

% Find any axes
subplots = getappdata(gcf,'SubplotGrid');

% Font size
% Not subplots
if numel(subplots) == 0
    set(gca,'FontSize',fontsize)
    set(findall(gcf,'type','text'),'FontSize',fontsize)
% Subplots
else
    [r,c] = size(subplots); 
    for ri = 1:r
        for ci = 1:c
            axInfo = get_axis_info(subplots(ri,ci),mm);
            set(axInfo.ax,'fontsize',fontsize);
            if numel(axInfo.leg) > 0 
                set(axInfo.ax,'fontsize',fontsize);
            end
        end
    end
end


% % Get the current position and scale of the figure in pixels
set(gca,'LooseInset',[0 0 0 0]);
if ~mm
    set(gcf,'Units','Pixels');
else
    set(gcf,'Units','centimeters');
end
pFig = get(gcf,'Position');


pFig(1) = pFig(1) - (width - pFig(3))/2;     % Left offset
if pFig(1) < 0
    pFig(1) = 0;
end
pFig(2) = pFig(2) - (height - pFig(4));       % Right offset
if pFig(2) < 0
    pFig(2) = 0;
end
pFig(3) = width;                       % Width
pFig(4) = height;                       % Height


newPos = pFig;
newPos(1) = padding(1);
newPos(2) = padding(2);
newPos(3) = newPos(3) - padding(1) - padding(3);
newPos(4) = newPos(4) - padding(2) - padding(4);

% Not subplots
if numel(subplots) == 0
    axs = findall(gcf,'type','axes');
    ax = axs(1);
    axi = get_axis_info(ax,mm);
    move_axis(axi, newPos);
    
% Subplots
else
    [r,c] = size(subplots); 
    
    for ri = 1:r
        for ci = 1:c
            axInfo = get_axis_info(subplots(ri,ci),mm);
            set(axInfo.ax,'fontsize',fontsize);
            if numel(axInfo.leg) > 0 
                set(axInfo.ax,'fontsize',fontsize);
            end
        end
    end
    
    height = (newPos(4) - (subpad*(r-1))) / r;
    width = (newPos(3) - (subpad*(c-1))) / c;
    
    p = newPos;
    for ri = 1:r
        for ci = 1:c
            p(1) = newPos(1) + (width + subpad)*(ci-1);
            p(2) = newPos(2) + (height + subpad)*(ri-1);
            p(3) = width;
            p(4) = height;
            axInfo = get_axis_info(subplots(ri,ci),mm);
            move_axis(axInfo,p);
        end
    end
end
setappdata(gcf, 'SubplotGrid',subplots);
set(gcf,'Position',pFig);
axi = 0;
            

end


function [ axInfo ] = get_axis_info(ax,mm)

% Get the legend
leg = getappdata(ax,'LegendPeerHandle');
legPos = [];
legSide = 0;

% Get the axes position
if mm == 0
    set(ax,'Units','Pixels');
else
    set(ax,'Units','centimeters');
end
tInset = get(ax,'TightInset');
axisPos = get(ax,'Position');

% Get the tight bounding box around the axes
boundPos = axisPos;
boundPos(1) = boundPos(1) - tInset(1);              % Left edge
boundPos(2) = boundPos(2) - tInset(2);              % Bottom edge
boundPos(3) = boundPos(3) + tInset(1) + tInset(3);  % Width
boundPos(4) = boundPos(4) + tInset(2) + tInset(4);  % Height

if numel(leg) > 0
    leg = leg(1);
    
    % Get the legend position
    if mm == 0
        set(leg,'Units','Pixels');
    else
        set(leg,'Units','centimeters');
    end
    legPos = get(leg,'Position');
    
    % Check if legend outside and add inset    
    %   Left edge     Legend left edge 
    t = boundPos(1)     - legPos(1);
    if t > 0
        boundPos(1) = legPos(1);    % Left edge now legend left edge
        tInset(1) = tInset(1) + t;  % Inset increased
        legSide = 1;
    end
    %   Legend right edge          Right edge
    t = (legPos(1) + legPos(3)) - (boundPos(1) + boundPos(3));
    if t > 0
        boundPos(3) = boundPos(3) + t;  % Width now goes to legend right edge
        tInset(3) = tInset(3) + t;      % Inset increased
        legSide = 3;
    end
    %   Bottom edge     Legend bottom edge 
    t = boundPos(2)     - legPos(2);
    if t > 0
        boundPos(2) = legPos(2);    % Left edge now legend bottom edge
        tInset(2) = tInset(2) + t;  % Inset increased
        legSide = 2;
    end
    %   Legend right edge          Right edge
    t = (legPos(2) + legPos(4)) - (boundPos(2) + boundPos(4));
    if t > 0
        boundPos(4) = boundPos(4) + t;  % Width now goes to legend right edge
        tInset(4) = tInset(4) + t;      % Inset increased
        legSide = 4;
    end
end

axInfo.ax = ax;
axInfo.axisPos = axisPos;
axInfo.boundPos = boundPos;
axInfo.tightInset = tInset;
axInfo.leg = leg;
axInfo.legPos = legPos;
axInfo.legSide = legSide;
end


function [ axInfo ] = move_axis( axInfo, newPos )

    rShift = newPos(3) - axInfo.boundPos(3);
    uShift = newPos(4) - axInfo.boundPos(4);

    % Initial shifts to keep legend centred properly
    %axInfo
    if numel(axInfo.leg) > 0 ...
       && axInfo.legSide == 3 ...
       && axInfo.legPos(2) > axInfo.axisPos(2) + 10 ...
       && (axInfo.legPos(2) + axInfo.legPos(4))  < (axInfo.axisPos(2) + axInfo.axisPos(2) - 10)
        axInfo.legPos(2) = axInfo.legPos(2) + 0.5*uShift;
    end
    if numel(axInfo.leg) > 0 ...
       && axInfo.legSide == 4 ...
       && axInfo.legPos(1) > axInfo.axisPos(1) + 10 ...
       && (axInfo.legPos(1) + axInfo.legPos(3))  < (axInfo.axisPos(1) + axInfo.axisPos(3) - 10)
        axInfo.legPos(1) = axInfo.legPos(1) + 0.5*rShift;
    end
    
    % Shift left edge
    lShift = newPos(1) - axInfo.boundPos(1);
    axInfo.boundPos(1) = newPos(1);
    axInfo.axisPos(1) = newPos(1) + axInfo.tightInset(1);
    if numel(axInfo.leg) > 0
        axInfo.legPos(1) = axInfo.legPos(1) + lShift;
    end
    
    % Shift bottom edge
    bShift = newPos(2) - axInfo.boundPos(2);
    axInfo.boundPos(2) = newPos(2);
    axInfo.axisPos(2) = newPos(2) + axInfo.tightInset(2);
    if numel(axInfo.leg) > 0
        axInfo.legPos(2) = axInfo.legPos(2) + bShift;
    end
    
    % Shift right edge
    axInfo.boundPos(3) = newPos(3);
    axInfo.axisPos(3) = newPos(3) - axInfo.tightInset(1) - axInfo.tightInset(3);
    if numel(axInfo.leg) > 0 %&& axInfo.legSide == 3
        axInfo.legPos(1) = axInfo.legPos(1) + rShift;
    end
    
    % Shift upper edge
    axInfo.boundPos(4) = newPos(4);
    axInfo.axisPos(4) = newPos(4) - axInfo.tightInset(2) - axInfo.tightInset(4);
    if numel(axInfo.leg) > 0 %&& axInfo.legSide == 4
        axInfo.legPos(2) = axInfo.legPos(2) + uShift;
    end
    
    
    
    % Move
    set(axInfo.ax,'Position',axInfo.axisPos);
    if numel(axInfo.leg) > 0
        set(axInfo.leg,'Position',axInfo.legPos);
    end
    %axInfo

end



