clc
clear
warning off
tic
success=0;

data=xlsread('样本6.xlsx');
S1=data(:,1:4);
Y1=data(:,5:15);%数据导入，包括S-x,y的坐标，Y-坐标对应的值
% 初始化参数
theta = [10 10 10 10]; lob = [1E-1 1E-1 1E-1 1E-1]; upb = [20 20 20 20];
[dmodel, perf] = dacefit(S1, Y1, @regpoly2, @corrgauss, theta, lob, upb);%建立模型
for cs=1:200
%速度更新参数
c1=1.5;
c2=1.5;
w1=0.9;
w2=0.4;
maxgen=100; %迭代次数
sizepop=20; %种群规模

%输入实测频率及位移
f1=13.75924007;
f2=32.36150951;
f3=57.5260883;
f4=61.12926202;
f5=101.8984548;
f6=113.2772432;
d1=0.038937255;
d2=0.0798484;
d3=0.0908588;
d4=0.0719816;
d5=0.033604515;
%个体和速度最大最小值
%popmax=254100000;popmin=161700000;
popmax=1;popmin=0;
Vmax=1;Vmin=-1;
%pop1=zeros(20,16);
for i=1:sizepop
    %随机产生一个种群
    %pop(i,:)=46200000*rands(1,16)+207900000*ones(1,16);
    pop(i,:)=rand(1,4);
    V(i,:)=rands(1,4);
    %生成matlab产生的试验数据
    
    mu=[pop(i,1) pop(i,2) pop(i,3) pop(i,4)];
    yreal=predictor(mu,dmodel);%5*50
    
    fre=yreal(1:6);
    disp=yreal(7:11);
    
    
      
   

        a1=1;
        a2=1;
        a3=1;
        a4=1;
        a5=1;
        a6=1;

        b1=1;
        b2=1;
        b3=1;
        b4=1;
        b5=1;
        fitness(i)=1.*(a1*(((fre(1,1)-f1)./f1)^2)+a2*(((fre(2,1)-f2)./f2)^2)+...
    a3*(((fre(3,1)-f3)./f3)^2)+a4*(((fre(4,1)-f4)./f4)^2)+a5*(((fre(5,1)-f5)./f5)^2)+...
    a6*(((fre(6,1)-f6)./f6)^2))+1.*(b1*(((disp(1,1)+d1)./d1)^2)+b2*(((disp(2,1)+d2)./d2)^2)+...
    b3*(((disp(3,1)+d3)./d3)^2)+b4*(((disp(4,1)+d4)./d4)^2)+b5*(((disp(5,1)+d5)./d5)^2));
end

[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);     %群体极值位置
gbest=pop;     %个体极值位置
fitnessgbest=fitness;     %个体极值适应度值
fitnesszbest=bestfitness;    %群体极值适应度值

%迭代寻优
for i=1:maxgen
    w=w1-((w1-w2)*i)/100;%线性递减惯性权重
    %粒子位置和速度更新
    for j=1:sizepop
        
        %速度更新
        V(j,:)=w*V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %粒子更新
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        mu=[pop(j,1) pop(j,2) pop(j,3) pop(j,4)];
        yreal=predictor(mu,dmodel);%5*50
    
        fre=yreal(1:6);
        disp=yreal(7:11);
        


       fitness(j)=1.*(a1*(((fre(1,1)-f1)./f1)^2)+a2*(((fre(2,1)-f2)./f2)^2)+...
    a3*(((fre(3,1)-f3)./f3)^2)+a4*(((fre(4,1)-f4)./f4)^2)+a5*(((fre(5,1)-f5)./f5)^2)+...
    a6*(((fre(6,1)-f6)./f6)^2))+1.*(b1*(((disp(1,1)+d1)./d1)^2)+b2*(((disp(2,1)+d2)./d2)^2)+...
    b3*(((disp(3,1)+d3)./d3)^2)+b4*(((disp(4,1)+d4)./d4)^2)+b5*(((disp(5,1)+d5)./d5)^2));
    end
    
    %个体极值和群体极值更新
    for j=1:sizepop
        
        %个体极值更新
        if fitness(j)<fitnessgbest(j)
            gbest(j,:)=pop(j,:);
            fitnessgbest(j)=fitness(j);
        end
        
        %群体极值更新
        if fitness(j)<fitnesszbest
            zbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    
    %每代最优值记录到yy数组中
    yy(i)=fitnesszbest;
end
    fitness_set(cs)=fitnesszbest;
    set_of_zbest(cs,:)=zbest;
   
    dis(:,cs)=disp;
    fr(:,cs)=fre;
    if (1.617+1.386*zbest(1,1)-2.079)/2.079<=0.05&&(1.617+1.386*zbest(1,2)-1.848)/1.848<=0.05...
        &&(1.617+1.386*zbest(1,3)-1.9635)/1.9635<=0.05&&(4914+4212*zbest(1,4)-6318)/6318<=0.05;
        success=success+1;
        fitness_set_success(cs)=fitnesszbest;
        
    end
end
%画出每代最优个体适应度值
plot(yy)
title('最优个体适应度值','fontsize',12);
xlabel('进化次数','fontsize',12);
ylabel('适应度值','fontsize',12);
Zbest_average=sum(set_of_zbest)./200;
mu_average=[Zbest_average(1,1) Zbest_average(1,2) Zbest_average(1,3) Zbest_average(1,4)];
yreal_average=predictor(mu_average,dmodel);
S1=std(set_of_zbest(:,1),1);
S2=std(set_of_zbest(:,2),1);
S3=std(set_of_zbest(:,3),1);
S4=std(set_of_zbest(:,4),1);
S=[S1 S2 S3 S4];  
Ddis_average=sum(transpose(dis))./200;
dis_average=transpose(Ddis_average);
Ffr_average=sum(transpose(fr))./200;
fr_average=transpose(Ffr_average);
D1=std(dis(1,:),1);
D2=std(dis(2,:),1);
D3=std(dis(3,:),1);
D4=std(dis(4,:),1);
D5=std(dis(5,:),1);
D=[D1 D2 D3 D4 D5];
DDD=transpose(D);
F1=std(fr(1,:),1);
F2=std(fr(2,:),1);
F3=std(fr(3,:),1);
F4=std(fr(4,:),1);
F5=std(fr(5,:),1);
F6=std(fr(6,:),1);
F=[F1 F2 F3 F4 F5 F6];
FFF=transpose(F);