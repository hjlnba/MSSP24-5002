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
%      mu=[0.5 0.5 0.5 0.5 0.5 0.5 0.5];
%      yreal=predictor(mu,dmodel);

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
f5=106.7507622;
f6=118.6713977;
d1=0.0407914;
d2=0.0798484;
d3=0.0908588;
d4=0.0719816;
d5=0.0352047;
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
    
    mu0=[pop(i,1) pop(i,2) pop(i,3) pop(i,4)];
    yreal0=predictor(mu0,dmodel);%5*50
    
    fre0=yreal0(1:6);
    disp0=yreal0(7:11);
    
   
    %计算粒子适应度
    %fre由ansys模型跑出
        t=1/f1+1/f2+1/f3+1/f4+1/f5+1/f6;
        t1=1/f1;t2=1/f2;t3=1/f3;t4=1/f4;t5=1/f5;t6=1/f6;
        T1=(t1/t);
        T2=(t2/t);
        T3=(t3/t);
        T4=(t4/t);
        T5=(t5/t);
        T6=(t6/t);
       

        r=d1+d2+d3+d4+d5;
        r1=d1;r2=d2;r3=d3;r4=d4;r5=d5;
        R1=r1/r;
        R2=r2/r;
        R3=r3/r;
        R4=r4/r;
        R5=r5/r;
    
        
         fitness(i)=1.*(T1*(((fre0(1,1)-f1)./f1)^2)+T2*(((fre0(2,1)-f2)./f2)^2)+...
    T3*(((fre0(3,1)-f3)./f3)^2)+T4*(((fre0(4,1)-f4)./f4)^2)+T5*(((fre0(5,1)-f5)./f5)^2)+...
    T6*(((fre0(6,1)-f6)./f6)^2))+1.*(R1*(((disp0(1,1)+d1)./d1)^2)+R2*(((disp0(2,1)+d2)./d2)^2)+...
    R3*(((disp0(3,1)+d3)./d3)^2)+R4*(((disp0(4,1)+d4)./d4)^2)+R5*(((disp0(5,1)+d5)./d5)^2));
end

[bestfitness bestindex]=min(fitness);
zbest0=pop(bestindex,:);     %群体极值位置
gbest=pop;     %个体极值位置
fitnessgbest=fitness;     %个体极值适应度值
fitnesszbest=bestfitness;    %群体极值适应度值

%迭代寻优
for i=1:maxgen
    w=w1-((w1-w2)*i)/100;%线性递减惯性权重
    %粒子位置和速度更新
    for j=1:sizepop
        
        %速度更新
        V(j,:)=w*V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest0-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %粒子更新
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
    mu0=[pop(j,1) pop(j,2) pop(j,3) pop(j,4)];
    yreal0=predictor(mu0,dmodel);%5*50
    
    fre0=yreal0(1:6);
    disp0=yreal0(7:11);
        
        fitness(j)=1.*(T1*(((fre0(1,1)-f1)./f1)^2)+T2*(((fre0(2,1)-f2)./f2)^2)+...
    T3*(((fre0(3,1)-f3)./f3)^2)+T4*(((fre0(4,1)-f4)./f4)^2)+T5*(((fre0(5,1)-f5)./f5)^2)+...
    T6*(((fre0(6,1)-f6)./f6)^2))+1.*(R1*(((disp0(1,1)+d1)./d1)^2)+R2*(((disp0(2,1)+d2)./d2)^2)+...
    R3*(((disp0(3,1)+d3)./d3)^2)+R4*(((disp0(4,1)+d4)./d4)^2)+R5*(((disp0(5,1)+d5)./d5)^2));
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
            zbest0=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    yy(i)=fitnesszbest;
end
    dis0(:,cs)=disp0;
    fr0(:,cs)=fre0;
figure(1);
plot(yy);
hold on;
%title('最优个体适应度值','fontsize',12);
xlabel('进化次数','fontsize',12);
ylabel('适应度值','fontsize',12);
% xlabel('迭代次数','Fontname','SimSun','Fontsize',12);
% ylabel('适应度值','Fontname','SimSun','Fontsize',12);
% 设置坐标轴数字的字体为 Times New Roman
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);    


%第二阶段修正
for i=1:sizepop
    %随机产生一个种群
    %pop(i,:)=46200000*rands(1,16)+207900000*ones(1,16);
    pop(i,:)=rand(1,4);
    pop(1,:)=[zbest0(1,1) zbest0(1,2) zbest0(1,3) zbest0(1,4)];
    V(i,:)=rands(1,4);
    %生成matlab产生的试验数据
    
    mu=[pop(i,1) pop(i,2) pop(i,3) pop(i,4)];
    yreal=predictor(mu,dmodel);%5*50
    
    fre=yreal(1:6);
    disp=yreal(7:11);
    
    
      
   
    %计算粒子适应度
    %fre由ansys模型跑出
        t=1/f1+1/f2+1/f3+1/f4+1/f5;
        t=1/f1+1/f2+1/f3+1/f4+1/f5+1/f6;
        t1=1/f1;t2=1/f2;t3=1/f3;t4=1/f4;t5=1/f5;t6=1/f6;
        T1=(t1/t);
        T2=(t2/t);
        T3=(t3/t);
        T4=(t4/t);
        T5=(t5/t);
        T6=(t6/t);

        r=d1+d2+d3+d4+d5;
        r1=d1;r2=d2;r3=d3;r4=d4;r5=d5;
        R1=r1/r;
        R2=r2/r;
        R3=r3/r;
        R4=r4/r;
        R5=r5/r;
    
        E1f=abs((fre(1,1)-f1)./f1);
        E2f=abs((fre(2,1)-f2)./f2);
        E3f=abs((fre(3,1)-f3)./f3);
        E4f=abs((fre(4,1)-f4)./f4);
        E5f=abs((fre(5,1)-f5)./f5);
        E6f=abs((fre(6,1)-f6)./f6);
        Ef=E1f+E2f+E3f+E4f+E5f+E6f;
        a1=T1*(1-E1f/Ef)^2;
        a2=T2*(1-E2f/Ef)^2;
        a3=T3*(1-E3f/Ef)^2;
        a4=T4*(1-E4f/Ef)^2;
        a5=T5*(1-E5f/Ef)^2;
        a6=T6*(1-E6f/Ef)^2;
        A=a1+a2+a3+a4+a5+a6;
        A1=a1/A;A2=a2/A;A3=a3/A;A4=a4/A;A5=a5/A;A6=a6/A;%xiu'gai'mu'biao'han'shu
        E1d=abs((disp(1,1)+d1)./d1);
        E2d=abs((disp(2,1)+d2)./d2);
        E3d=abs((disp(3,1)+d3)./d3);
        E4d=abs((disp(4,1)+d4)./d4);
        E5d=abs((disp(5,1)+d5)./d5);
        Ed=E1d+E2d+E3d+E4d+E5d;
        b1=R1*(1-E1d/Ed)^2;
        b2=R2*(1-E2d/Ed)^2;
        b3=R3*(1-E3d/Ed)^2;
        b4=R4*(1-E4d/Ed)^2;
        b5=R5*(1-E5d/Ed)^2;
        B=b1+b2+b3+b4+b5;
        B1=b1/B;B2=b2/B;B3=b3/B;B4=b4/B;B5=b5/B;
        
        fitness(i)=1.*(A1*(((fre(1,1)-f1)./f1)^2)+A2*(((fre(2,1)-f2)./f2)^2)+...
    A3*(((fre(3,1)-f3)./f3)^2)+A4*(((fre(4,1)-f4)./f4)^2)+A5*(((fre(5,1)-f5)./f5)^2)+...
    A6*(((fre(6,1)-f6)./f6)^2))+1.*(B1*(((disp(1,1)+d1)./d1)^2)+B2*(((disp(2,1)+d2)./d2)^2)+...
   B3*(((disp(3,1)+d3)./d3)^2)+B4*(((disp(4,1)+d4)./d4)^2)+B5*(((disp(5,1)+d5)./d5)^2));
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
        

        E1f=abs((fre(1,1)-f1)./f1);
        E2f=abs((fre(2,1)-f2)./f2);
        E3f=abs((fre(3,1)-f3)./f3);
        E4f=abs((fre(4,1)-f4)./f4);
        E5f=abs((fre(5,1)-f5)./f5);
        E6f=abs((fre(6,1)-f6)./f6);
        Ef=E1f+E2f+E3f+E4f+E5f+E6f;
        a1=T1*(1-E1f/Ef)^2;
        a2=T2*(1-E2f/Ef)^2;
        a3=T3*(1-E3f/Ef)^2;
        a4=T4*(1-E4f/Ef)^2;
        a5=T5*(1-E5f/Ef)^2;
        a6=T6*(1-E6f/Ef)^2;
        A=a1+a2+a3+a4+a5+a6;
        A1=a1/A;A2=a2/A;A3=a3/A;A4=a4/A;A5=a5/A;A6=a6/A;%xiu'gai'mu'biao'han'shu
        E1d=abs((disp(1,1)+d1)./d1);
        E2d=abs((disp(2,1)+d2)./d2);
        E3d=abs((disp(3,1)+d3)./d3);
        E4d=abs((disp(4,1)+d4)./d4);
        E5d=abs((disp(5,1)+d5)./d5);
        Ed=E1d+E2d+E3d+E4d+E5d;
        b1=R1*(1-E1d/Ed)^2;
        b2=R2*(1-E2d/Ed)^2;
        b3=R3*(1-E3d/Ed)^2;
        b4=R4*(1-E4d/Ed)^2;
        b5=R5*(1-E5d/Ed)^2;
        B=b1+b2+b3+b4+b5;
        B1=b1/B;B2=b2/B;B3=b3/B;B4=b4/B;B5=b5/B;
        
        fitness(j)=1.*(A1*(((fre(1,1)-f1)./f1)^2)+A2*(((fre(2,1)-f2)./f2)^2)+...
    A3*(((fre(3,1)-f3)./f3)^2)+A4*(((fre(4,1)-f4)./f4)^2)+A5*(((fre(5,1)-f5)./f5)^2)+...
    A6*(((fre(6,1)-f6)./f6)^2))+1.*(B1*(((disp(1,1)+d1)./d1)^2)+B2*(((disp(2,1)+d2)./d2)^2)+...
   B3*(((disp(3,1)+d3)./d3)^2)+B4*(((disp(4,1)+d4)./d4)^2)+B5*(((disp(5,1)+d5)./d5)^2));
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
    yyy(i)=fitnesszbest;
     
    
    
    
    
end
    fitness_set(cs)=fitnesszbest;
    set_of_zbest(cs,:)=zbest;
    
    dis(:,cs)=disp;
    fr(:,cs)=fre;
    %if (1.617+1.386*zbest(1,1)-2.079)/2.079<=0.05&&(1.617+1.386*zbest(1,2)-1.848)/1.848<=0.05...
    %    &&(1.617+1.386*zbest(1,3)-1.9635)/1.9635<=0.05&&(4914+4212*zbest(1,4)-6318)/6318<=0.05;
    %   success=success+1;
    %    fitness_set_success(cs)=fitnesszbest;
    %end
%画出每代最优个体适应度值
   figure(2)
   plot(yyy)
   title('最优个体适应度值','fontsize',12);
   xlabel('进化次数','fontsize',12);
   ylabel('适应度值','fontsize',12);

%xlabel('迭代次数','Fontname','SimSun','Fontsize',12);
%ylabel('适应度值','Fontname','SimSun','Fontsize',12);
   % 设置坐标轴数字的字体为 Times New Roman
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
hold on;

end
% figure(3)
% plot(yyy)
% title('最优个体适应度值','fontsize',12);
% xlabel('进化次数','fontsize',12);
% ylabel('适应度值','fontsize',12);
hold on;


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

