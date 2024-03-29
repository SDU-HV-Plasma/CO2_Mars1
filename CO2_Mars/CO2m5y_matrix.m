function [B,G,Y]=CO2m5y_matrix(X,nsp)

global B G

% 建立有向关系
% 建立相互关系矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        B{i,j}=intersect(find(X(i,:)>0),find(X(j,:)<0));
    end
    i=i+1;
end

% % 在该矩阵中对有粒子参与碰撞但粒子数目增加或维持不变的反应进行添加
B{2,1}=[2,88,B{2,1}];
B{3,1}=[4,107,B{3,1}];
B{4,1}=[2,4,5,28,47,67,90,106,108,B{4,1}];
B{6,1}=[106,162,B{6,1}];
B{8,1}=[1,B{8,1}];
B{11,1}=[28,46,89,B{11,1}];
B{12,1}=[33,46,47,B{12,1}];
B{14,1}=[67,68,B{14,1}];
B{16,1}=[88,89,90,91,106,162,B{16,1}];
B{17,1}=[107,108,109,B{17,1}];
B{3,2}=[24,B{3,2}];
B{13,2}=[58,B{13,2}];
B{1,3}=[17,41,87,B{1,3}];
B{4,3}=[17,64,B{4,3}];
B{5,3}=[57,87,101,B{5,3}];
B{6,3}=[19,B{6,3}];
B{10,3}=[20,B{10,3}];
B{12,3}=[100,138,B{12,3}];
B{13,3}=[59,B{13,3}];
B{15,3}=[75,77,B{15,3}];
B{19,3}=[120,B{19,3}];
B{3,4}=[105,B{3,4}];
B{1,10}=[22,B{1,10}];
B{3,10}=[60,B{3,10}];
B{4,10}=[22,B{4,10}];
B{5,10}=[14,B{5,10}];
B{21,10}=[128,B{21,10}];
% % M的添加仅为findequation.m提供帮助
B{2,9}=[32,151];
B{3,9}=[34,66,137];
B{4,9}=[40,92,169,170];
B{5,9}=[66,169];
B{7,9}=137;
B{8,9}=158;
B{10,9}=[148,149,158];
B{11,9}=40;
B{14,9}=151;
B{18,9}=119;
B{22,9}=133;
B{24,9}=148;
B{25,9}=149;
B{26,9}=145;

% 拓扑矩阵G
G=B;
G{3,1}=[7,19,75,G{3,1}];
G{9,1}=[34,92];
G{1,2}=[28,G{1,2}];
G{3,2}=[20,G{3,2}];
G{9,2}=[148,149,40];
G{10,2}=[128,G{10,2}];
G{1,3}=[5,162,G{1,3}];
G{9,3}=[133,170];
G{10,3}=[14,G{10,3}];
G{2,4}=[16,24,G{2,4}];
G{3,4}=[19,20,65,100,G{3,4}];
G{9,4}=32;
G{10,4}=[14,60,G{10,4}];
G{1,5}=[4,106,G{1,5}];
G{3,5}=[75,G{3,5}];
G{2,6}=[58,G{2,6}];
G{3,6}=[17,57,59,77,G{3,6}];
G{4,6}=[105,G{4,6}];
G{10,6}=[22,42,G{10,6}];
G{3,7}=[101,120,G{3,7}];
G{9,7}=133;
G{9,8}=145;
G{1,10}=[2,G{1,10}];
G{2,10}=[58,G{2,10}];
G{3,10}=[59,120,G{3,10}];
G{9,10}=[119,145];
G{9,11}=32;
G{10,11}=[128,G{10,11}];
G{3,12}=[57,132,G{3,12}];
G{1,12}=[164,G{1,12}];
G{9,12}=[34,66,119];
G{3,15}=[87,G{3,15}];
G{9,15}=[66,169];
G{3,16}=[100,101,G{3,16}];
G{4,16}=[105,G{4,16}];
G{1,16}=[163,G{1,16}];
G{9,16}=[92,169];
G{1,17}=[86,G{1,17}];
G{9,22}=137;
G{9,25}=148;
G{9,26}=[149,158];
G{9,24}=151;

% 建立有向矩阵
i=1;
while i<nsp+1
    for j=1:nsp
        Y1(i,j)=sum(X(i,B{i,j}));
    end
    Y2(i)=sum(X(i,X(i,:)>0));
    i=i+1;
end
Y=Y1./Y2';

end

