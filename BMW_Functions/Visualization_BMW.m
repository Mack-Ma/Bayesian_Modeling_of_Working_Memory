%% Visualization_BMW
% 
% Visualize the data and statistics of data

function Figure_BMW=Visualization_BMW(data,type,config)

switch type
    case 'ModelSelectionColormap'
        if nargin==2 || ~isfield(config,'Criteria')
            error('Please designate which criterion (e.g. BIC) you are willing to employ...\nType ""BMW(''Manual'')"" for further information')
        end
        if ~isfield(data,config.Criteria)
            error('The criterion you designated was not detected in the data file...')
        end
        if strcmp(config.Criteria,'BIC') || strcmp(config.Criteria,'LME_BS') || strcmp(config.Criteria,'LME_GHM')
            useLME=1;
        else
            useLME=0;
        end
        
        % Load data
        COI=data.(config.Criteria);
        COI_weight=COI.([config.Criteria,'_weight']);
        COI_Group_PP=COI.('Group_PP');
        if useLME==1 && isfield(COI,'EP')
            useRFXBMS=1;
            COI_EP=COI.('EP');
            COI_ModelFreq=COI.('ModelFreq');
        else
            useRFXBMS=0;
        end
        Nsubj=size(COI_weight,1);
        Nmodel=size(COI_weight,2);
        
        % Start constructing the figure
        Figure_BMW=figure;
        % Create axes
        if Nsubj>20
            Position1=[0.12 0.3 0.7 0.6];
            UnitWidth=0.03;
        else
            UnitWidth=0.03;
            Position1=[0.12 0.3 0.7 UnitWidth*Nsubj];
        end
        axes1 = axes('Parent',Figure_BMW,...
            'Position',Position1);
        hold(axes1,'on');
        colormap(parula(256));
        % Create image
        image(COI_weight,'Parent',axes1,'CDataMapping','scaled');
        % Set colormap limits
        caxis([0,1]);
        % Create ylabel
        ylabel('Participant');
        % Set the X-limits of the axes
        xlim(axes1,[0.5 Nmodel+.5]);
        % Set the Y-limits of the axes
        ylim(axes1,[0.5 Nsubj+.5]);
        box(axes1,'on');
        axis(axes1,'ij');
        % Set the remaining axes properties
        if Nsubj>2
            YTick_subj=[1, floor(Nsubj/2), Nsubj];
        else
            YTick_subj=zeros(1,0);
        end
        set(axes1,'Layer','top','XTick',zeros(1,0),'YTick',YTick_subj);
        % Create colorbar
        colorbar('peer',axes1,'Position',[Position1(1)+Position1(3)+0.05, Position1(2), 0.05, Position1(4)]);
        % Create axes
        Position2=[Position1(1) Position1(2)-0.05 Position1(3) UnitWidth];
        axes2 = axes('Parent',Figure_BMW,...
            'Position',Position2);
        hold(axes2,'on');
        % Create image
        image([1 Nmodel],1,COI_Group_PP,'Parent',axes2,'CDataMapping','scaled');
        % Set colormap limits
        caxis([0,1]);
        if useRFXBMS==0
            % Create xlabel
            xlabel('Model');
            XTick2=1:Nmodel;
        else
            XTick2=[];
        end
        % Set the X-limits of the axes
        xlim(axes2,[0.5 Nmodel+.5]);
        % Set the Y-limits of the axes
        ylim(axes2,[0.5 1.5]);
        box(axes2,'on');
        axis(axes2,'ij');
        % Set the remaining axes properties
        set(axes2,'Layer','top','XTick',XTick2,'XTickLabel',...
            num2cell(1:Nmodel),'YTick',1,...
            'YTickLabel',{'GPP'});
        if useRFXBMS==1
            % Model Frequency here
            Position3=[Position2(1) Position2(2)-0.05 Position2(3) UnitWidth];
            axes3 = axes('Parent',Figure_BMW,...
                'Position',Position3);
            hold(axes3,'on');
            % Create image
            image([1 Nmodel],1,COI_ModelFreq,'Parent',axes3,'CDataMapping','scaled');
            % Set colormap limits
            caxis([0,1]);
            % Set the X-limits of the axes
            xlim(axes3,[0.5 Nmodel+.5]);
            % Set the Y-limits of the axes
            ylim(axes3,[0.5 1.5]);
            box(axes3,'on');
            axis(axes3,'ij');
            % Set the remaining axes properties
            set(axes3,'Layer','top','XTick',zeros(1,0),'XTickLabel',...
                zeros(1,0),'YTick',1,...
                'YTickLabel',{'r'});
            % EP here
            % Create axes
            Position4=[Position3(1) Position3(2)-0.05 Position3(3) UnitWidth];
            axes4 = axes('Parent',Figure_BMW,...
                'Position',Position4);
            hold(axes4,'on');
            % Create image
            image([1 Nmodel],1,COI_EP,'Parent',axes4,'CDataMapping','scaled');
            % Set colormap limits
            caxis([0,1]);
            % Create xlabel
            xlabel('Model');
            % Set the X-limits of the axes
            xlim(axes4,[0.5 Nmodel+.5]);
            % Set the Y-limits of the axes
            ylim(axes4,[0.5 1.5]);
            box(axes4,'on');
            axis(axes4,'ij');
            % Set the remaining axes properties
            set(axes4,'Layer','top','XTick',1:Nmodel,'XTickLabel',...
                num2cell(1:Nmodel),'YTick',1,...
                'YTickLabel',{'EP'});
            if Nsubj<=20
                FigHeight=Nsubj*0.03;
                FigWidth=0.3;
            else
                FigHeight=0.6;
                FigWidth=0.3;
            end
        else
            if Nsubj<=20
                FigHeight=Nsubj*0.03-0.05;
                FigWidth=0.3;
            else
                FigHeight=0.55;
                FigWidth=0.3;
            end
        end
        FigSize=[0.2,0.25,FigWidth,FigHeight];
        set(gcf,'unit','normalized','position',FigSize)
        
    case 'HistRecallError'
        
    case 'HistCatError'
        
    case 'MCMCChain'
        
    case 'MCMCDensity'
        
end

end