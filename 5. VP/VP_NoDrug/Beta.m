a=0;
b=100;
e=100;

alpha = 4;  % 设置Beta分布参数alpha为1
beta = 4;   % 设置Beta分布参数beta为1
y=Beta_distribution(a,b,alpha,beta,e);


histogram(y,10)

%plot(x,y)   % 绘制Beta分布曲线
title('Standard Beta Distribution')   % 添加标题
xlabel('x')   % 添加x轴标签
ylabel('f(x)')   % 添加y轴标签
% 
% Z = 1;
% W = 0:0.01:1;
% B = beta(Z,W);
% plot(W,B)