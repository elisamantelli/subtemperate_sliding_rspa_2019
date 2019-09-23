function hout = plotpanel(ii,parameters,labels,lims,fontsize)
%hout = plotpanel(ii,parameters,labels,lims,fontsize)
%Set up subplots with minimal whitespace. Input variables are
%   parameters: structure, containing the following optional fields
%       n_rows:     number of rows in grid of panels (default 1)
%       n_cols:  number of columns in grid of panels (default 3)
%       col_tickoffset: scaled (to total width) space between columns of
%                   panels (default 0.025)
%       row_tickoffset: scaled (to total height) space between rows of
%                   panels (default 0.025)
%       col_offset: scaled (to total width) space to the left of the columns to
%                   accommodate labels (default 0.1)
%       row_offset: scaled (to total height) space below the tows to
%                   accommodate labels (default 0.1)
%       box:        boolean, 'true' forces 'box on' (default false)
%       equal:      boolean, 'true' sets 'axis equal' (default false)
%       hold:       boolean, 'true' sets 'hold on' (default false)
%       log:        creats logscale on both axes
%       semilogx:   creates logscale on x axis
%       semilogy:   creates logscale on y axis
%   labels: structure containing the following optional fields
%       xlab:       x-axis label, displayed only if panel specified through
%                   ii is at the bottom of a column (default '')
%       ylab:       y-axis label, displayed only if panel specified through
%                   ii is at the left-hand end of a row (default '')
%       panel_lab:  panel label, placed top left hand corner (default '')
%   lims: structure containing the following optional fields
%       x:          one-by-two-vector defining x limits
%       y:          one-by-two-vector defining y limits
%   fontsize:       Font size to be used (in pts)
%   ii:   identifies which panel to generate, counting across from top left
%The output variable is
%   hout:           axes handle for axes created

if nargin < 2, parameters.dummy = []; end
if isfield(parameters,'n_rows')
    n_rows = parameters.n_rows;
else
    n_rows = 1;
end
if isfield(parameters,'n_cols')
    n_cols = parameters.n_cols;
else
    n_cols = 3;
end
if isfield(parameters,'col_offset')
    col_offset = parameters.col_offset;
else
    col_offset = .1;
end
if isfield(parameters,'row_offset')
    row_offset = parameters.row_offset;
else
    row_offset = .1;
end
if isfield(parameters,'col_tickoffset')
    col_tickoffset = parameters.col_tickoffset;
else
    col_tickoffset = .025;
end
if isfield(parameters,'row_tickoffset')
    row_tickoffset = parameters.row_tickoffset;
else
    row_tickoffset = .025;
end
if isfield(parameters,'box')
    boxflag = parameters.box;
else
    boxflag = false;
end
if isfield(parameters,'equal')
    equalflag = parameters.equal;
else
    equalflag = false;
end
if isfield(parameters,'hold')
    holdflag = parameters.hold;
else
    holdflag = false;
end
if isfield(parameters,'log')
    logflag = parameters.log;
else
    logflag = false;
end
if isfield(parameters,'semilogx')
    semilogxflag = parameters.semilogx;
else
    semilogxflag = false;
end
if isfield(parameters,'semilogy')
    semilogyflag = parameters.semilogy;
else
    semilogyflag = false;
end

boxwidth = (1-col_offset)/n_cols;
boxheight = (1-row_offset)/n_rows;

plotwidth = boxwidth-col_tickoffset;
plotheight = boxheight-row_tickoffset;

col = 1+mod(ii-1,n_cols);
row = n_rows-floor((ii-1)/n_cols);

if row<0 && col<0, error('plot outside figure'), end

hout = subplot('position',[col_offset+(col-1)*boxwidth row_offset+(row-1)*boxheight plotwidth plotheight]);

%set plot type
if logflag
    loglog(NaN,NaN)
elseif semilogxflag
    semilogx(NaN,NaN)
elseif semilogyflag
    semilogy(NaN,NaN)
end

%Font size
if nargin >= 5;
    set(gca,'FontSize',fontsize)
end

%finish labelling etc
if boxflag, box on, end
if equalflag, axis equal, end
if holdflag, hold on, end

if nargin < 3
    xlab  = ''; ylab = ''; panel_lab = '';
else
    if isfield(labels,'xlab')
        xlab = labels.xlab;
    else
        xlab = '';
    end
    if isfield(labels,'ylab')
        ylab = labels.ylab;
    else
        ylab = '';
    end
    if isfield(labels,'panel_lab')
        panel_lab = labels.panel_lab;
    else
        panel_lab = '';
    end
end
if nargin >= 4
    if ~isfield(lims,'x')
        lims.x(1) = -1; lims.x(2) = 1;
    end
    if ~isfield(lims,'y')
        lims.y(1) = -1; lims.y(2) = 1;
    end
end



if col==1 
    ylabel(ylab,'Interpreter', 'latex')
else
    set(hout,'YTickLabel',{})
end
if row==1
    xlabel(xlab,'Interpreter', 'latex')
else
    set(hout,'XTickLabel',{})
end
if exist('lims')
    xlim(lims.x);
    ylim(lims.y);
    text(lims.x(1)+.9*(lims.x(2)-lims.x(1)),lims.y(2)-.1*(lims.y(2)-lims.y(1)),panel_lab,'FontSize',fontsize, 'Interpreter', 'latex')
end

end
