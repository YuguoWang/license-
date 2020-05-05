clear all;
clc



file_path =  'D:\课件\大三下\综合课程设计\Part2设计图片库\保存图片\';% 图像文件夹路径
img_path_list = dir(strcat(file_path,'*.jpg'));%获取该文件夹中所有.jpg格式的图像
img_num = length(img_path_list);%获取图像总数
chepai=[];
%iii=1;
if img_num > 0 %有满足条件的图像
        for pn = 1:img_num %逐一读取图像
            image_name = img_path_list(pn).name;% 图像名
            dw =  imread(strcat(file_path,image_name));%读取图像
            fprintf('%d %s\n',pn,strcat(file_path,image_name));% 显示正在处理的图像名
            


%[filename, pathname]=uigetfile('*jpg','select a picture');
%图像数据读入 
%I=imread([pathname,filename]);%采集图象

[heng, shu] = size(dw);
%求输入图像经retinex算法补偿后的图像
subplot(4, 7, 1),imshow(dw);title('1.原图定位')


%letters = zeros(40,20,7);
I2 = rgb2gray(dw);
%subplot(232); imshow(I2); title("2.灰度图像");
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1.边缘算子、直线检测和旋转
I3 = edge(I2,'Sobel','horizontal');
%subplot(233);imshow(I3);title("边缘检测");
se = [1 1 1;1 1 1;1 1 1];
I4 = imdilate(I3,se);
%subplot(234);imshow(I4);title("膨胀");
 
[H,T,R] = hough(I4,'Theta',-89:89);
ccc = max(H);
[value, rot_theta] = max(ccc);

if rot_theta<=90
     img_correction = imrotate(dw ,rot_theta,'bilinear', 'loose');
     cut_theta = rot_theta;
else
     if rot_theta<=180 && rot_theta>90
    img_correction = imrotate(dw ,180+rot_theta,'bilinear', 'loose');
    cut_theta = 180-rot_theta;
     end
end
subplot(472), imshow(img_correction);title('2.矫正后的图像');
[xinheng,xinshu]=size(img_correction);
cut_length = ceil((shu*sin(cut_theta*pi/180))/3);
img_cut=imcrop(img_correction,[1, cut_length, xinshu/3, (xinheng-2*cut_length)]);
img_cut = imresize(img_cut,[150 450],'bilinear');
subplot(473), imshow(img_cut);title('3.裁剪后的图像');

Gray=rgb2gray(img_cut); %double型
Gray=imadjust(Gray,[0.2,0.7],[]);

subplot(4, 7, 4),imshow(Gray),title('4.车牌灰度图像')

%Gray = medfilt2(Gray);
l=graythresh(Gray); %利用Ostu法获取阈值
Bin=imbinarize(Gray,l); %logical型
Med= medfilt2(Bin);

[m,n]=size(Bin);
sumd = sum(sum(Med(floor(m/4):floor(m*3/4),:)));

if sumd/(m*n*0.5) > 0.4
    Gray=rgb2gray(img_correction);
    Gray=imadjust(Gray,[0.4,0.8],[]);
    subplot(4, 7, 4),imshow(Gray),title('4.车牌灰度图像')
    l=graythresh(Gray); %利用Ostu法获取阈值
    Bin=imbinarize(Gray,l); %logical型
    Med= medfilt2(Bin);
end

%Med = medfilt2(Gray); %中值滤波
subplot(4, 7, 5),imshow(Med),title('5.中值滤波后')

SE=strel('line',2,90);%创建一个线条形状
%Ero = imclose(Med,SE);
Ero=imerode(Med,SE);
subplot(4, 7, 6),imshow(Ero),title('6.腐蚀图像');


%pic=imerode(pic,SE);
Close = imclose(Ero,SE);%闭运算, 连接各字符之间的缝隙
%Close = imerode(pic,SE);
[m1, n1] = size(Close);
Del2=bwareaopen(Close,floor(m1*n1/800));%删除面积小于一定大小的区域
subplot(4, 7, 7),imshow(Del2),title('7.车牌二值图像');

%Del3=chubuqiege(Del2);
subplot(4, 7, 8),imshow(Del2),title('8.初步切割后图片')
Del2 = imresize(Del2,[140 440],'bilinear');
d=qiege(Del2);
subplot(4, 7, 9),imshow(d),title('9.切割后图片')
d = imresize(d,[90 410],'bilinear');
%d = bwareaopen(d,90);%删除面积小于一定大小的区域
subplot(4, 7, 10),imshow(d),title('10.统一大小后的图片')
if(sum(sum(d(:,371:380)))>870)
    pos=390;
    while(sum(d(:,pos))>2)
        pos= pos-1;
    end
    d(:,pos:410) = 0;
end

d_adjust = chuizhiqingxiejiaozheng(d);
d_adjust = chubuqiege(d_adjust);
subplot(4, 7, 11),imshow(d_adjust),title('11.校正后图片')
d_adjust = imresize(d_adjust,[90 410],'bilinear');
d2=d_adjust(:,37:410);
SE1=strel('line',2,90);%创建一个线条形状
for i=1:2
d2 = imerode(d2,SE1);
end
d_adjust(:,37:410)=d2;
L = bwlabel(d_adjust);
nummax=max(max(L));
num=nummax;
liantongyu = zeros(1,nummax); %创建保存连通域的和向量
liantongyuyou = zeros(1,nummax); %创建保存连通域最右点向量
liantongyuzuo = zeros(1,nummax); %创建保存连通域最左点向量
liantongyuheight = zeros(1,nummax); %创建保存连通域高度的向量
liantongyutop = zeros(1,nummax); %创建保存连通域最高点和向量
liantongyubut = zeros(1,nummax); %创建保存连通域最低点向量
for i =1:nummax  %初始化连通域最低坐标
    liantongyubut(i) = 90;
end
for i=1:90      %计算各连通域的像素点数和
    for j=1:410
        if(L(i,j)>0)
            liantongyu(L(i,j)) = liantongyu(L(i,j))+1;
        end
    end
end
for j = 1:410   %计算各连通域最右像素点的纵坐标
    for i=1:90
        if(L(i,j)>0)
            liantongyuyou(L(i,j)) = j;
        end
    end
end
j=410;
while(j>=1)   %计算各连通域最左像素点的纵坐标
    for i=1:90
        if(L(i,j)>0)
            liantongyuzuo(L(i,j)) = j;
        end
    end
    j=j-1;
end

for i = 1:90  %计算各连通域最高和最低点坐标
    for j=1:410
        if(L(i,j)>0)
            if(i<liantongyubut(L(i,j)))
                liantongyubut(L(i,j)) = i;
            end
            if(i>liantongyutop(L(i,j)))
                liantongyutop(L(i,j)) = i;
            end
        end
    end
end
for i = 1:nummax  %计算各连通域高度
    liantongyuheight(i) = liantongyutop(i) - liantongyubut(i);
end
i=1;
while(i<=nummax)   %去除边框
    if(liantongyu(i)<300 && (liantongyuzuo(i)<30 || liantongyuzuo(i)>380) && liantongyuheight(i)>55)
        for k = 1:90
            for j=1:410
                if(L(k,j)==i)
                    L(k,j) = 0;
                end
                if(L(k,j)>i)
                    L(k,j) = L(k,j)-1;
                end
            end
        end
        liantongyu(i)=[];
        liantongyuzuo(i)=[];
        liantongyuyou(i)=[];
        liantongyuheight(i)=[];
        liantongyutop(i)=[];
        liantongyubut(i)=[];
        nummax = nummax - 1;
    end
    i=i+1;
end
i=1;
while(i<=nummax)   %去除小的连通域
    if(liantongyu(i)<300 && liantongyuzuo(i)>80)
        for k = 1:90
            for j=1:410
                if(L(k,j)==i)
                    L(k,j) = 0;
                end
                if(L(k,j)>i)
                    L(k,j) = L(k,j)-1;
                end
            end
        end
        liantongyu(i)=[];
        liantongyuzuo(i)=[];
        liantongyuyou(i)=[];
        liantongyuheight(i)=[];
        liantongyutop(i)=[];
        liantongyubut(i)=[];
        nummax = nummax - 1;
    end
    i=i+1;
end
num = nummax;




subplot(4, 7, 12),imshow(L),title('12.连通域图')

for i = 1:6  %判断是字符并进行切割
    while(~(liantongyu(num)>695 && liantongyuheight(num)>50 && (liantongyuyou(num)-liantongyuzuo(num))>=7 && (liantongyuyou(num)-liantongyuzuo(num))<70) && ~(liantongyu(num)>580 && liantongyuheight(num)>30 && (liantongyuyou(num)-liantongyuzuo(num))>=30 && (liantongyuyou(num)-liantongyuzuo(num))<70) )
        num=num-1;
    end
    Ifin = imcrop(L,[liantongyuzuo(num),1,liantongyuyou(num)-liantongyuzuo(num),90]);
    Qie(:,:,8-i) = imresize(Ifin,[90 50],'bilinear');
    num=num-1;
end
Ifin = imcrop(L,[liantongyuzuo(num+1)-58,1,50,90]);
Qie(:,:,1) = imresize(Ifin,[90 50],'bilinear');
for i=1:7  %把图像重新变回二值图像
    for j=1:90
        for k=1:50
            if(Qie(j,k,i)>=1)
                Qie(j,k,i) = 1;
            else
                Qie(j,k,i) = 0;
            end
        end
    end
end
for i=3:7
    max1(i) = max(sum(Qie(:,:,i)));
end
max2 = max(max1);
for i=3:7     %找出1，把1的宽度恢复
    sumd=sum(sum(Qie(:,:,i)));
    maxd=max(sum(Qie(:,:,i)));
    sumoflie=sum(Qie(:,:,i));
    sum1 = 0;sum2 = 0;
    for j = 1:50
        if(sumoflie(j)>5)
            sum1=sum1+1;
        end
        if(sumoflie(j)>(0.6*maxd))
            sum2=sum2+1;
        end
    end
    
    if(sumd/4500>0.4 && (sum2/sum1)>=0.55 && (maxd/max2)>0.88 && sumoflie(25)>(0.85*maxd))
        s1 = imresize(Qie(:,:,i),[90 14],'bilinear');
        Qie(:,:,i) = 0;
        Qie(:,(19:32),i) = s1(:,:);
    end
end


        
for i = 1:7
    subplot(4, 7, (21+i)),imshow(Qie(:,:,i));
end
pic = zeros(90,350);
for(i=1:7)
    pic(:,(1+(i-1)*50):i*50) = Qie(:,:,i);
end

% for iii=1:7
% savepathname = sprintf('D:%s课件%s大三下%s综合课程设计%sPart2设计图片库%s分割图片%s%d图.jpg','\','\','\','\','\','\',pn*7-(7-iii));
% imwrite(Qie(:,:,iii),savepathname);
% end

subplot(4, 7, 13),imshow(pic),title('13.分割后图片粘连图') 

Qie1=imresize(Qie(:,:,1),[90,42]);
Qie1=[zeros(90,8) Qie1(:,:,1)];
Qie(:,:,1)=Qie1;
for i=1:7
    cha{i} = [Qie(:,:,i)];
end
cha{1} = imresize(cha{1}, [32, 16]);
for i = 2:7
%    [a,b] = size(cha{i});
%    if a > 3 * b
%        tmp = round((a - b)/2);
%        tmpI = cha{i};
%        cha{i} = [zeros(a,tmp), tmpI, zeros(a,tmp)];
%    end
%   cha{i} = imresize(cha{i}, [32, 16]);
cha{i} = imresize(cha{i}, [40, 32]);
end
character_im=cha;
characters = total_shibie(character_im);


chepai=[chepai;characters];
%iii=iii+1;
%D:\课件\大三下\综合课程设计\Part2设计图片库\分割图片
% for iii=1:7
% savepathname = sprintf('D:%s课件%s大三下%s综合课程设计%sPart2设计图片库%s分割图片1%s%d图.jpg','\','\','\','\','\','\',pn*7-(7-iii));
% imwrite(dw,savepathname);
% end

end
end


%初步切割
    function pic = chubuqiege(d)
    [m,n]=size(d);
    top=1;bottom=m;left=1;right=n;  

    while sum(d(top,:))<=3 && top<m     %切割出白色区域（横切）
        top=top+1;
    end
    while sum(d(bottom,:))<=3 && bottom>1   %同上
        bottom=bottom-1;
    end
    while sum(d(:,left))<=3 && left<n        %切割出白区域（纵切）
        left=left+1;
    end
    while sum(d(:,right))<=3 && right>1
        right=right-1;
    end
    dd=right-left;
    hh=bottom-top;
    cut=imcrop(d,[left top dd hh]);%裁剪图片
    cut = imresize(cut, [m,n], 'bilinear');%保持图片大小一致
    %pic1=imclearborder(cut);%去除与边界相连的部分
    [m1, n1] = size(cut);
    pic=cut;%删除面积较小的区域
    end

    %切割
    function cut=qiege(d)
    [m,n]=size(d);
    bottom1= floor(m/2) ;top1= floor(m/2) ;left = 1; right= n ;  

    while sum(d(bottom1,:))>20 && bottom1<m     %切割出白色区域（横切）
        bottom1=bottom1+1;
    end
    while sum(d(top1,:))>20 && top1>1   %同上
        top1=top1-1;
    end
    while sum(d(:,left))<=0 && left<n        %切割出白区域（纵切）
        left=left+1;
    end
    while sum(d(:,right))<=0 && right>1
        right=right-1;
    end
    dd=right-left;
    hh=bottom1-top1;
    cut=imcrop(d,[left, top1, dd, hh]);%裁剪图片
    end

    %分割字符
    function Qie = Fenge(d)
        m_1 = zeros(1,90);%创建空的横投影向量
        n_1 = zeros(1,410);%创建空的竖投影向量
        for i = 1:410
        n_1(i) = sum(d(:, i)); %竖投影
        end
        for i = 1:90
            m_1(i) = sum(d(i, :)); %横投影
        end
        
        left1 = zeros(1,6);
        left = 1;
        flag = 0;
        while flag == 0
            while n_1(left)<1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            left1(1) = left;
            while n_1(left)>=1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            left1(2) = left;
            while n_1(left)<1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            while n_1(left)>=1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            left1(3) = left;
            while n_1(left)<1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            while n_1(left)>=1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            left1(4) = left;
            while n_1(left)<1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            while n_1(left)>=1
                left = left +1 ;
            end
            
            if left >200
                flag = 1;
                break;
            end
            
            left1(5) = left;
            left1(6) = left1(5) - left1(1);
            %subplot(4,4,1),plot(1:6, left1);
            if(left1(6) < 72 )
                d=imcrop(d,[left1(2) 1 (410 -left1(2))  90]);
                d = imresize(d, [90,410], 'bilinear');                
                flag = 1;
            end
        end

        for i = 1:410
                n_1(i) = sum(d(:, i)); %竖投影
        end
        for i = 1:90
                m_1(i) = sum(d(i, :)); %横投影
        end
        
        right = 1;
        while(n_1(410-right)<1)
            right = right+1;
        end
        right1 = right;
        while(n_1(410-right)>1)
            right = right +1;
        end
        right2 = right;
        if((right1 - right2)<18)
            n_1(right2:right1) = 0;
        end
        
        subplot(4,5,9),plot(1:410, n_1),title('竖投影')
        subplot(4,5,10),plot(1:90, m_1),title('横投影')
        subplot(4, 7, 11),imshow(d),title('11.最终切割图')
        
        lettermid = zeros(1,7);
        lettermid(1) = 20;
        letterright = 390;
        MINnum = min(n_1(20:390))+2;
        for i = 1:6
            while(n_1(letterright)<=MINnum)
                letterright = letterright - 1;
            end
            r = letterright;
            while(n_1(letterright)>MINnum)
                letterright = letterright - 1;
            end
            l = letterright;
            lettermid(8-i) = floor((r+l)/2);
        end
        subplot(4,4,13),plot(1:7, lettermid),title('字的中心')
        cutsite = zeros(1,7);
        cutsite(7) = 410;
        cutsite(1:6) = lettermid(1:6);
        for i = 1:20
            if((n_1(410-i)<=n_1(cutsite(7))) && n_1(410-i)<4)
                cutsite(7) = 410-i;
            end
        end
        for j = 1:6
            for i = lettermid(j):lettermid(j+1)
                if(n_1(i)<=n_1(cutsite(j)))
                    cutsite(j) = i;
                end
            end
        end
        subplot(4,4,14),plot(1:7, cutsite),title('初始切割点')
        cutsiteright = cutsite;
        while((n_1(cutsiteright(7))<=n_1(cutsite(7))) && cutsiteright(7)<410)   %对于最右切点的修正
            cutsiteright(7) = cutsiteright(7) + 1;
        end
        cutsite(7) = floor((cutsite(7)+cutsiteright(7))/2);
        
        for i = 1:6
            while(n_1(cutsiteright(i))<=n_1(cutsite(i)))   %对于其他切点的修正
                cutsiteright(i) = cutsiteright(i) - 1;
            end
            cutsite(i) = floor((cutsite(i)+cutsiteright(i))/2);
        end
        
        subplot(4,4,15),plot(1:7, cutsite),title('修正切割点')
        lengthofcut = zeros(1,7);
        lengthofcut(1)= cutsite(1);
        for i = 2:7
            lengthofcut(i) = cutsite(i) - cutsite(i-1);
        end
        subplot(4,4,16),plot(1:7, lengthofcut),title('各字符长度')
        Ifin = imcrop(d,[1,1,lengthofcut(1),90]);
        Qie(:,:,1) = imresize(Ifin,[90 50],'bilinear');
        for i = 1:6
            Ifin = imcrop(d,[cutsite(i),1,lengthofcut(i+1),90]);
            Qie(:,:,i+1) = imresize(Ifin,[90 50],'bilinear');
        end
    end

   %倾斜矫正
    function d_adjust = chuizhiqingxiejiaozheng(I)
        [m,n] = size(I);
        kk = 1;
        for i = 1:n
            for j = 1:m
                if I(j,i) ~= 0
                    X1(kk) = i;
                    Y1(kk) = j;
                    kk = kk+1;
                end
            end
        end
        Zaa = [X1 ;Y1];
        
        Eu = zeros(1,n);
        Ed = zeros(1,n);
        Rk = zeros(1,2*n-1);
        thatas = 0 ;
        thatat = 0.013;
        k = 1;
        N = 10;
        thata(1) = 1;
        flag = 0;
    while(flag ~=1 && k<N)
        
        kk = 1;
        for i = 1:n
            for j = 1:m
                if I(j,i) ~= 0
                    X(kk) = i;
                    Y(kk) = j;
                    kk = kk+1;
                end
            end
        end
        Z = [X ;Y];   %生成坐标矩阵
        Z = Z';
        
        for i = 1:n    %计算上包络
            j = 1;
            while(I(j,i)<1 && j < m)
                j = j+1;
            end
            Eu(i) = j;         
        end

        for i =1:n     %计算下包络
            j = 1;
            while(I(m-j+1,i)<1 && j < m)
                j = j+1;
            end
            Ed(i) = j;
        end
        
        subplot(4,4,9),plot(1:n, Eu),title('10.上包络')
        subplot(4,4,11),plot(1:n, Ed),title('11.下包络')
        
        fft1 = fft(Eu);fft11 = fft(Ed);
        subplot(4,5,16),plot(1:length(fft1),fft1),title('15.上包络频谱')
        subplot(4,5,18),plot(1:length(fft1),fft1),title('15.下包络频谱')
        lpf = myfilter;
        Eu2 = filter(lpf,Eu);Ed2 = filter(lpf,Ed);
        fft2 = fft(Eu2);fft22 = fft(Ed2);
        subplot(4,4,10),plot(1:n, Eu2),title('10.滤波后上包络')
        subplot(4,4,12),plot(1:n, Ed2),title('10.滤波后下包络')
        subplot(4,5,17),plot(1:length(fft2),fft2),title('15.滤波后上包络频谱')
        subplot(4,5,19),plot(1:length(fft2),fft2),title('15.滤波后上包络频谱')
        
        
        Euu = zeros(1,3*n-2);
        Edd = zeros(1,3*n-2);
        for i =1:n
        Euu(n-1+i) = Eu2(i);
        end
        for i =1:n
        Edd(n-1+i) = Ed2(i);
        end
        for i = 1:2*n-1   %计算互相关
            for j = 1:n
                Rk(i) = Rk(i) + Euu(j+n-1)*Edd(j+i-1);
            end
        end
        
        %subplot(4,7,12),plot(-n+1:n-1, Rk),title('12.互相关')    
        
        
        [x ,n0] = max(Rk);    %找出互相关最大值坐标
        thata(k) = atan((n0-n)/m); %计算错切角
        if abs(thata(k)) > thatat    %如果thtak大于设定的thatat，说明仍然倾斜，应该继续
            if k<N   %判断是否仍小于最大迭代次数
                thatas = sum(thata);    %计算总
                Z = Z*[1 0;-tan(thata(k)) 1];
                I = zeros(m,n);
                for i = 1:length(X)
                    if ceil(Z(i,2))<=0      %排除矫正完坐标超出索引值
                        Z(i,2) = 1;
                    end
                    if ceil(Z(i,1))<=0      %排除矫正完坐标超出索引值
                        Z(i,1) = 1;
                    end
                    %if ceil(Z(i,2))>n      %排除矫正完坐标超出索引值
                       % Z(i,2) = n;
                    %end
                    %if ceil(Z(i,1))>n      %排除矫正完坐标超出索引值
                        %Z(i,1) = n;
                   % end
                    I(ceil(Z(i,2)),ceil(Z(i,1))) = 1;
                end
                k = k+1;
            else
                thatas = sum(thata);
                flag = 1;
            end
        else
            thatas = sum(thata);
            flag = 1;
        end
    end
 Zaa = Zaa'; 
 Zaa = Zaa*[1 0;-tan(thatas) 1];
 Zaa = Zaa';
 
 I_ad = zeros(m,n);
 for i = 1:length(X)
                    if ceil(Zaa(2,i))<=0      %排除矫正完坐标超出索引值
                        Zaa(2,i) = 1;
                    end
                    if ceil(Zaa(1,i))<=0      %排除矫正完坐标超出索引值
                        Zaa(1,i) = 1;
                    end
                   % if ceil(Zaa(2,i))>n      %排除矫正完坐标超出索引值
                       % Zaa(2,i) = n;
                   % end
                   % if ceil(Zaa(1,i))>n      %排除矫正完坐标超出索引值
                      %  Zaa(1,i) = n;
                   % end
     I_ad(ceil(Zaa(2,i)),ceil(Zaa(1,i))) = 1;
 end
  
  %subplot(4, 4, 13),imshow(I_ad),title("倾斜矫正后图像");
        
        d_adjust= I_ad;
    end 
