function [RT, MT] = calc_rt_mt
% cond_name = {'MotorOnly-L','MotorOnly-S','CueOnly-L','CueOnly-S',...
%              'BothRep-L','BothRep-S','NonRep-L','NonRep-S','Non-Interest'};
glm = 2;
k = 1;
for s=[1:3 5 6 8:14]
    R = construct_dsgmat(sprintf('S%02d',s),glm);
    R2 = construct_dsgmat(sprintf('R%02d',s),glm);
    R = combine_behavdata(R,R2);
    for c=1:8

        temp = R.MT(R.cond==c);
        MT(k,c) = mean(temp(temp~=0));
        temp = R.RT(R.cond==c);
        RT(k,c) = mean(temp(temp~=0));
    end
    k=k+1;
end

RT = RT(:,[7 1 3 5 8 2 4 6]);
MT = MT(:,[7 1 3 5 8 2 4 6]);
seqType = [1 1 1 1 2 2 2 2];
cond = [1:8];
repType = [1 2 3 4 1 2 3 4];
CAT.facecolor={
    [0.3 0.3 0.6],...
    [0.3 0.6 0.3],...
    [0.6 0.3 0.3],...
    [0.3 0.3 0.3]
    };

T.RT = reshape(RT,prod(size(RT)),1);
T.MT = reshape(MT,prod(size(MT)),1);
T.seqtype = [ones(48,1);2*ones(48,1)];
T.reptype = reshape(repmat(repType,12,1),prod(size(RT)),1);
figure;
set(gcf,'color','w');
subplot(1,2,1)
legend_text = {'NRep','Seq-Rep','Cue-Rep','Both-Rep'};
title_text = 'Reaction time';
x_coord = barplot([T.seqtype],T.RT,'split',[T.reptype],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                    'leg',legend_text);
% ylim([400 520]);      
        set(gca,'YLim',[300 520],'TickLength',[0.01 0.01],'YTick',[300:100:500],...
            'XTick',(x_coord(2:4:end)+x_coord(3:4:end))/2,'XTickLabel',{'Letter','Spatial'},...
            'XLim',[x_coord(1)-1.2 x_coord(end)+1.2],'FontSize',12,'LineWidth',1,'FontName', 'Arial');
        ylabel('Reaction time (ms)','FontSize',12,'FontName','Arial');
subplot(1,2,2)
legend_text = {'NRep','Seq-Rep','Cue-Rep','Both-Rep'};
% title_text = 'Reaction time';
x_coord = barplot([T.seqtype],T.MT,'split',[T.reptype],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                    'leg',legend_text);
% ylim([1100 1300]);      
        set(gca,'YLim',[1000 1300],'TickLength',[0.01 0.01],'YTick',[1100:100:1300],...
            'XTick',(x_coord(2:4:end)+x_coord(3:4:end))/2,'XTickLabel',{'Letter','Spatial'},...
            'XLim',[x_coord(1)-1.2 x_coord(end)+1.2],'FontSize',12,'LineWidth',1,'FontName', 'Arial');
        ylabel('Movement time (ms)','FontSize',12,'FontName','Arial');

% Create bar plot
figure;
subplot(1,2,1);
bar([1 2 3 4], mean(RT(:,1:4)), 'b');
hold on;
errorbar(mean(RT(:,1:4)), std(RT(:,1:4))/sqrt(size(RT,1)));
bar([7 8 9 10], mean(RT(:,5:8)), 'r');
errorbar(mean(RT(:,5:8)),std(RT(:,5:8))/sqrt(size(RT,1)));
ylim([400 550]);
subplot(1,2,2);
bar([1 2 3 4], mean(MT(:,1:4)), 'b');
hold on;
errorbar(mean(MT(:,1:4)), std(MT(:,1:4))/sqrt(size(MT,1)));
bar([7 8 9 10], mean(MT(:,5:8)), 'r');
errorbar(mean(MT(:,5:8)),std(MT(:,5:8))/sqrt(size(MT,1)));
ylim([1000 1300]);
% 
% % Overlay individual data points
% for i = 1:size(data, 2)
%     scatter(repmat(i, 1, size(data, 1)), data(:, i), 'filled', 'MarkerFaceColor', 'b');
% end
% 
% % Customize the plot
% xlabel('Conditions');
% ylabel('Values');
% title('Bar and Errorbar Plot with Individual Data Points');
% grid on;
% hold off;
% 
% % figure()
% % set(gcf,'color','w');
% % legend_text = {'NRep','Seq-Rep','Cue-Rep','Both-Rep'};
% title_text = 'Repetition effects';
% x_coord = barplot([cond],RT,'split',[seqType],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
%                     'leg',legend_text);