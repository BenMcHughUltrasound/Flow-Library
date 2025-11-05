
function figure = STANDARDIZE_FIGURE(fig_comps)
    %========================================================
    % INCLUDE PLOT_STANDARDS
    
    PS = PLOT_STANDARDS();
        
    
    %========================================================
    % GET FIGURE HANDLE
    
    fig = fig_comps;
    
    
    %========================================================
    % DIFFERENT STYLES FOR NON-TILED VS TILED
    
    % Non-Tiled Plot
        numel(findobj(fig, 'type', 'tile')) == 0;
        %========================================================
        % PLOT FULLSCREEN BY DEFAULT
        % Make the plot 'maximized'. The plot opens up in fullscreen mode by default
        %fig.WindowState = 'maximized';
        %========================================================
        % GET CURRENT FIGURE AXES
        % get axes for current figure
        ax = fig.CurrentAxes;
        %========================================================
        % SET PROPERTIES FOR X-AXIS AND Y-AXIS
        % Set Properties for X and Y Axis Numbers
        ax.FontName = PS.AxisNumbersFontName;
        ax.FontSize = PS.AxisNumbersFontSize;
        ax.FontWeight = 'bold';
        % Set Properties for X, Y and Z Labels
        ax.XLabel.FontName = PS.AxisFont;
        ax.YLabel.FontName = PS.AxisFont;
        ax.ZLabel.FontName = PS.AxisFont;
        ax.XLabel.FontSize = PS.AxisFontSize;
        ax.YLabel.FontSize = PS.AxisFontSize;
        ax.ZLabel.FontSize = PS.AxisFontSize;
        ax.XLabel.FontWeight = 'bold';
        ax.YLabel.FontWeight = 'bold';
        ax.ZLabel.FontWeight = 'bold';
        
        if length(ax.XLabel.String) > 4
            if all(ax.XLabel.String([1, 2, end - 1, end]) == '$$$$')
                ax.XLabel.Interpreter = 'latex';
            end
        end
        if length(ax.YLabel.String) > 4
            if all(ax.YLabel.String([1, 2, end - 1, end]) == '$$$$')
                ax.YLabel.Interpreter = 'latex';
            end
        end
        if length(ax.ZLabel.String) > 4
            if all(ax.ZLabel.String([1, 2, end - 1, end]) == '$$$$')
                ax.ZLabel.Interpreter = 'latex';
            end
        end
        %========================================================
        % SET PROPERTIES FOR TITLE
        ax.Title.FontName = PS.TitleFont;
        ax.Title.FontSize = PS.TitleFontSize;
        ax.Title.FontWeight = 'bold';
        if length(ax.Title.String) > 4
            if all(ax.Title.String([1, 2, end - 1, end]) == '$$$$')
                ax.Title.Interpreter = 'latex';
            end
        end
        %========================================================
        % SET PROPERTIES FOR LEGEND
        
        if (numel(ax.Legend) ~= 0)
            ax.Legend.FontName = PS.LegendFont;
            ax.Legend.FontSize = PS.LegendFontSize;
            for i_Legend = 1: length(ax.Legend.String)
                if length(ax.Legend.String{i_Legend}) > 4
                    if all(ax.Legend.String{i_Legend}([1, 2, end - 1, end]) == '$$$$')
                        ax.Legend.Interpreter = 'latex';
                        break
                    end
                end
            end
            ax.Legend.Box = 'on';
            ax.Legend.LineWidth = PS.DefaultLegendBoxLineWidth;
            ax.Legend.AutoUpdate = 'off';

        end
        %========================================================
        % SET PROPERTIES FOR CHILDREN OF AXES -> LINES, TEXT ON PLOT
        axisChildren = ax.Children;
        for i = 1:length(axisChildren)
            if isequal(axisChildren(i).Type, 'text')
                axisChildren(i).FontName = PS.PlotTextFont;
                axisChildren(i).FontSize = PS.PlotTextFontSize;
                if length(axisChildren(i).String) > 4
                    if all(axisChildren(i).String([1, 2, end - 1, end]) == '$$$$')
                        axisChildren(i).Interpreter = 'latex';
                    end
                end
            end
        end
        %========================================================
        % ADJUST AXES PROPERTIES
        
        ax.Box = 'on';
        ax.TickDir = 'out';
        ax.TickLength = [PS.AxisTickLength, PS.AxisTickLength];
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.ZMinorTick = 'on';
        ax.XColor = PS.AxisColor;
        ax.YColor = PS.AxisColor;
        ax.ZColor = PS.AxisColor;
        ax.XLabel.Color = PS.AxisLabelColor;
        ax.YLabel.Color = PS.AxisLabelColor;
        ax.ZLabel.Color = PS.AxisLabelColor;
        ax.LineWidth = PS.DefaultLineWidth;
    
    
        
    %========================================================
% HOW TO CITE THIS TOOLBOX? (Please support :) )
% Refer - https://in.mathworks.com/matlabcentral/answers/1575908-how-to-cite-matlab-file-exchange-file-using-bibtex?s_tid=es_ans_avr_ans_view#answers_list
% OR
% Include the following in your LaTeX .bib file:
% @misc
% {  atharva2021,
%    author = {Atharva Aalok}, 
%    title = {Professional Plots}, 
%    year = 2021
%    howpublished = "\url{https://www.mathworks.com/matlabcentral/fileexchange/100766-professional_plots}",
%    note = "[Online; accessed October 31, 2022]"
% }
% PLEASE SUPPORT ME BY CITING THIS TOOLBOX. THANK YOU!!!
end
