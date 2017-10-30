

randn('seed',0);
%%一维高斯函数
mu=0;
sigma=1;
x=-6:0.1:6;
y=normpdf(x,mu,sigma);
plot(x,y);
figure;