function [fitresult, gof] = createFit(t, u)
%CREATEFIT(T,U)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: t
%      Y 输出: u
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( t, u );

% 设置 fittype 和选项。
ft = fittype( 'poly4' );

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft );

% 绘制数据拟合图。
figure( 'Name', '左边界条件 ' );
h = plot( fitresult, xData, yData );
legend( h, 'u vs. t', '左边界条件', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 't', 'Interpreter', 'none' );
ylabel( 'u', 'Interpreter', 'none' );
grid on


