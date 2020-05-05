clear all;
clc



file_path =  'D:\�μ�\������\�ۺϿγ����\Part2���ͼƬ��\����ͼƬ\';% ͼ���ļ���·��
img_path_list = dir(strcat(file_path,'*.jpg'));%��ȡ���ļ���������.jpg��ʽ��ͼ��
img_num = length(img_path_list);%��ȡͼ������
chepai=[];
%iii=1;
if img_num > 0 %������������ͼ��
        for pn = 1:img_num %��һ��ȡͼ��
            image_name = img_path_list(pn).name;% ͼ����
            dw =  imread(strcat(file_path,image_name));%��ȡͼ��
            fprintf('%d %s\n',pn,strcat(file_path,image_name));% ��ʾ���ڴ����ͼ����
            


%[filename, pathname]=uigetfile('*jpg','select a picture');
%ͼ�����ݶ��� 
%I=imread([pathname,filename]);%�ɼ�ͼ��

[heng, shu] = size(dw);
%������ͼ��retinex�㷨�������ͼ��
subplot(4, 7, 1),imshow(dw);title('1.ԭͼ��λ')


%letters = zeros(40,20,7);
I2 = rgb2gray(dw);
%subplot(232); imshow(I2); title("2.�Ҷ�ͼ��");
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1.��Ե���ӡ�ֱ�߼�����ת
I3 = edge(I2,'Sobel','horizontal');
%subplot(233);imshow(I3);title("��Ե���");
se = [1 1 1;1 1 1;1 1 1];
I4 = imdilate(I3,se);
%subplot(234);imshow(I4);title("����");
 
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
subplot(472), imshow(img_correction);title('2.�������ͼ��');
[xinheng,xinshu]=size(img_correction);
cut_length = ceil((shu*sin(cut_theta*pi/180))/3);
img_cut=imcrop(img_correction,[1, cut_length, xinshu/3, (xinheng-2*cut_length)]);
img_cut = imresize(img_cut,[150 450],'bilinear');
subplot(473), imshow(img_cut);title('3.�ü����ͼ��');

Gray=rgb2gray(img_cut); %double��
Gray=imadjust(Gray,[0.2,0.7],[]);

subplot(4, 7, 4),imshow(Gray),title('4.���ƻҶ�ͼ��')

%Gray = medfilt2(Gray);
l=graythresh(Gray); %����Ostu����ȡ��ֵ
Bin=imbinarize(Gray,l); %logical��
Med= medfilt2(Bin);

[m,n]=size(Bin);
sumd = sum(sum(Med(floor(m/4):floor(m*3/4),:)));

if sumd/(m*n*0.5) > 0.4
    Gray=rgb2gray(img_correction);
    Gray=imadjust(Gray,[0.4,0.8],[]);
    subplot(4, 7, 4),imshow(Gray),title('4.���ƻҶ�ͼ��')
    l=graythresh(Gray); %����Ostu����ȡ��ֵ
    Bin=imbinarize(Gray,l); %logical��
    Med= medfilt2(Bin);
end

%Med = medfilt2(Gray); %��ֵ�˲�
subplot(4, 7, 5),imshow(Med),title('5.��ֵ�˲���')

SE=strel('line',2,90);%����һ��������״
%Ero = imclose(Med,SE);
Ero=imerode(Med,SE);
subplot(4, 7, 6),imshow(Ero),title('6.��ʴͼ��');


%pic=imerode(pic,SE);
Close = imclose(Ero,SE);%������, ���Ӹ��ַ�֮��ķ�϶
%Close = imerode(pic,SE);
[m1, n1] = size(Close);
Del2=bwareaopen(Close,floor(m1*n1/800));%ɾ�����С��һ����С������
subplot(4, 7, 7),imshow(Del2),title('7.���ƶ�ֵͼ��');

%Del3=chubuqiege(Del2);
subplot(4, 7, 8),imshow(Del2),title('8.�����и��ͼƬ')
Del2 = imresize(Del2,[140 440],'bilinear');
d=qiege(Del2);
subplot(4, 7, 9),imshow(d),title('9.�и��ͼƬ')
d = imresize(d,[90 410],'bilinear');
%d = bwareaopen(d,90);%ɾ�����С��һ����С������
subplot(4, 7, 10),imshow(d),title('10.ͳһ��С���ͼƬ')
if(sum(sum(d(:,371:380)))>870)
    pos=390;
    while(sum(d(:,pos))>2)
        pos= pos-1;
    end
    d(:,pos:410) = 0;
end

d_adjust = chuizhiqingxiejiaozheng(d);
d_adjust = chubuqiege(d_adjust);
subplot(4, 7, 11),imshow(d_adjust),title('11.У����ͼƬ')
d_adjust = imresize(d_adjust,[90 410],'bilinear');
d2=d_adjust(:,37:410);
SE1=strel('line',2,90);%����һ��������״
for i=1:2
d2 = imerode(d2,SE1);
end
d_adjust(:,37:410)=d2;
L = bwlabel(d_adjust);
nummax=max(max(L));
num=nummax;
liantongyu = zeros(1,nummax); %����������ͨ��ĺ�����
liantongyuyou = zeros(1,nummax); %����������ͨ�����ҵ�����
liantongyuzuo = zeros(1,nummax); %����������ͨ�����������
liantongyuheight = zeros(1,nummax); %����������ͨ��߶ȵ�����
liantongyutop = zeros(1,nummax); %����������ͨ����ߵ������
liantongyubut = zeros(1,nummax); %����������ͨ����͵�����
for i =1:nummax  %��ʼ����ͨ���������
    liantongyubut(i) = 90;
end
for i=1:90      %�������ͨ������ص�����
    for j=1:410
        if(L(i,j)>0)
            liantongyu(L(i,j)) = liantongyu(L(i,j))+1;
        end
    end
end
for j = 1:410   %�������ͨ���������ص��������
    for i=1:90
        if(L(i,j)>0)
            liantongyuyou(L(i,j)) = j;
        end
    end
end
j=410;
while(j>=1)   %�������ͨ���������ص��������
    for i=1:90
        if(L(i,j)>0)
            liantongyuzuo(L(i,j)) = j;
        end
    end
    j=j-1;
end

for i = 1:90  %�������ͨ����ߺ���͵�����
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
for i = 1:nummax  %�������ͨ��߶�
    liantongyuheight(i) = liantongyutop(i) - liantongyubut(i);
end
i=1;
while(i<=nummax)   %ȥ���߿�
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
while(i<=nummax)   %ȥ��С����ͨ��
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




subplot(4, 7, 12),imshow(L),title('12.��ͨ��ͼ')

for i = 1:6  %�ж����ַ��������и�
    while(~(liantongyu(num)>695 && liantongyuheight(num)>50 && (liantongyuyou(num)-liantongyuzuo(num))>=7 && (liantongyuyou(num)-liantongyuzuo(num))<70) && ~(liantongyu(num)>580 && liantongyuheight(num)>30 && (liantongyuyou(num)-liantongyuzuo(num))>=30 && (liantongyuyou(num)-liantongyuzuo(num))<70) )
        num=num-1;
    end
    Ifin = imcrop(L,[liantongyuzuo(num),1,liantongyuyou(num)-liantongyuzuo(num),90]);
    Qie(:,:,8-i) = imresize(Ifin,[90 50],'bilinear');
    num=num-1;
end
Ifin = imcrop(L,[liantongyuzuo(num+1)-58,1,50,90]);
Qie(:,:,1) = imresize(Ifin,[90 50],'bilinear');
for i=1:7  %��ͼ�����±�ض�ֵͼ��
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
for i=3:7     %�ҳ�1����1�Ŀ�Ȼָ�
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
% savepathname = sprintf('D:%s�μ�%s������%s�ۺϿγ����%sPart2���ͼƬ��%s�ָ�ͼƬ%s%dͼ.jpg','\','\','\','\','\','\',pn*7-(7-iii));
% imwrite(Qie(:,:,iii),savepathname);
% end

subplot(4, 7, 13),imshow(pic),title('13.�ָ��ͼƬճ��ͼ') 

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
%D:\�μ�\������\�ۺϿγ����\Part2���ͼƬ��\�ָ�ͼƬ
% for iii=1:7
% savepathname = sprintf('D:%s�μ�%s������%s�ۺϿγ����%sPart2���ͼƬ��%s�ָ�ͼƬ1%s%dͼ.jpg','\','\','\','\','\','\',pn*7-(7-iii));
% imwrite(dw,savepathname);
% end

end
end


%�����и�
    function pic = chubuqiege(d)
    [m,n]=size(d);
    top=1;bottom=m;left=1;right=n;  

    while sum(d(top,:))<=3 && top<m     %�и����ɫ���򣨺��У�
        top=top+1;
    end
    while sum(d(bottom,:))<=3 && bottom>1   %ͬ��
        bottom=bottom-1;
    end
    while sum(d(:,left))<=3 && left<n        %�и�����������У�
        left=left+1;
    end
    while sum(d(:,right))<=3 && right>1
        right=right-1;
    end
    dd=right-left;
    hh=bottom-top;
    cut=imcrop(d,[left top dd hh]);%�ü�ͼƬ
    cut = imresize(cut, [m,n], 'bilinear');%����ͼƬ��Сһ��
    %pic1=imclearborder(cut);%ȥ����߽������Ĳ���
    [m1, n1] = size(cut);
    pic=cut;%ɾ�������С������
    end

    %�и�
    function cut=qiege(d)
    [m,n]=size(d);
    bottom1= floor(m/2) ;top1= floor(m/2) ;left = 1; right= n ;  

    while sum(d(bottom1,:))>20 && bottom1<m     %�и����ɫ���򣨺��У�
        bottom1=bottom1+1;
    end
    while sum(d(top1,:))>20 && top1>1   %ͬ��
        top1=top1-1;
    end
    while sum(d(:,left))<=0 && left<n        %�и�����������У�
        left=left+1;
    end
    while sum(d(:,right))<=0 && right>1
        right=right-1;
    end
    dd=right-left;
    hh=bottom1-top1;
    cut=imcrop(d,[left, top1, dd, hh]);%�ü�ͼƬ
    end

    %�ָ��ַ�
    function Qie = Fenge(d)
        m_1 = zeros(1,90);%�����յĺ�ͶӰ����
        n_1 = zeros(1,410);%�����յ���ͶӰ����
        for i = 1:410
        n_1(i) = sum(d(:, i)); %��ͶӰ
        end
        for i = 1:90
            m_1(i) = sum(d(i, :)); %��ͶӰ
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
                n_1(i) = sum(d(:, i)); %��ͶӰ
        end
        for i = 1:90
                m_1(i) = sum(d(i, :)); %��ͶӰ
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
        
        subplot(4,5,9),plot(1:410, n_1),title('��ͶӰ')
        subplot(4,5,10),plot(1:90, m_1),title('��ͶӰ')
        subplot(4, 7, 11),imshow(d),title('11.�����и�ͼ')
        
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
        subplot(4,4,13),plot(1:7, lettermid),title('�ֵ�����')
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
        subplot(4,4,14),plot(1:7, cutsite),title('��ʼ�и��')
        cutsiteright = cutsite;
        while((n_1(cutsiteright(7))<=n_1(cutsite(7))) && cutsiteright(7)<410)   %���������е������
            cutsiteright(7) = cutsiteright(7) + 1;
        end
        cutsite(7) = floor((cutsite(7)+cutsiteright(7))/2);
        
        for i = 1:6
            while(n_1(cutsiteright(i))<=n_1(cutsite(i)))   %���������е������
                cutsiteright(i) = cutsiteright(i) - 1;
            end
            cutsite(i) = floor((cutsite(i)+cutsiteright(i))/2);
        end
        
        subplot(4,4,15),plot(1:7, cutsite),title('�����и��')
        lengthofcut = zeros(1,7);
        lengthofcut(1)= cutsite(1);
        for i = 2:7
            lengthofcut(i) = cutsite(i) - cutsite(i-1);
        end
        subplot(4,4,16),plot(1:7, lengthofcut),title('���ַ�����')
        Ifin = imcrop(d,[1,1,lengthofcut(1),90]);
        Qie(:,:,1) = imresize(Ifin,[90 50],'bilinear');
        for i = 1:6
            Ifin = imcrop(d,[cutsite(i),1,lengthofcut(i+1),90]);
            Qie(:,:,i+1) = imresize(Ifin,[90 50],'bilinear');
        end
    end

   %��б����
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
        Z = [X ;Y];   %�����������
        Z = Z';
        
        for i = 1:n    %�����ϰ���
            j = 1;
            while(I(j,i)<1 && j < m)
                j = j+1;
            end
            Eu(i) = j;         
        end

        for i =1:n     %�����°���
            j = 1;
            while(I(m-j+1,i)<1 && j < m)
                j = j+1;
            end
            Ed(i) = j;
        end
        
        subplot(4,4,9),plot(1:n, Eu),title('10.�ϰ���')
        subplot(4,4,11),plot(1:n, Ed),title('11.�°���')
        
        fft1 = fft(Eu);fft11 = fft(Ed);
        subplot(4,5,16),plot(1:length(fft1),fft1),title('15.�ϰ���Ƶ��')
        subplot(4,5,18),plot(1:length(fft1),fft1),title('15.�°���Ƶ��')
        lpf = myfilter;
        Eu2 = filter(lpf,Eu);Ed2 = filter(lpf,Ed);
        fft2 = fft(Eu2);fft22 = fft(Ed2);
        subplot(4,4,10),plot(1:n, Eu2),title('10.�˲����ϰ���')
        subplot(4,4,12),plot(1:n, Ed2),title('10.�˲����°���')
        subplot(4,5,17),plot(1:length(fft2),fft2),title('15.�˲����ϰ���Ƶ��')
        subplot(4,5,19),plot(1:length(fft2),fft2),title('15.�˲����ϰ���Ƶ��')
        
        
        Euu = zeros(1,3*n-2);
        Edd = zeros(1,3*n-2);
        for i =1:n
        Euu(n-1+i) = Eu2(i);
        end
        for i =1:n
        Edd(n-1+i) = Ed2(i);
        end
        for i = 1:2*n-1   %���㻥���
            for j = 1:n
                Rk(i) = Rk(i) + Euu(j+n-1)*Edd(j+i-1);
            end
        end
        
        %subplot(4,7,12),plot(-n+1:n-1, Rk),title('12.�����')    
        
        
        [x ,n0] = max(Rk);    %�ҳ���������ֵ����
        thata(k) = atan((n0-n)/m); %������н�
        if abs(thata(k)) > thatat    %���thtak�����趨��thatat��˵����Ȼ��б��Ӧ�ü���
            if k<N   %�ж��Ƿ���С������������
                thatas = sum(thata);    %������
                Z = Z*[1 0;-tan(thata(k)) 1];
                I = zeros(m,n);
                for i = 1:length(X)
                    if ceil(Z(i,2))<=0      %�ų����������곬������ֵ
                        Z(i,2) = 1;
                    end
                    if ceil(Z(i,1))<=0      %�ų����������곬������ֵ
                        Z(i,1) = 1;
                    end
                    %if ceil(Z(i,2))>n      %�ų����������곬������ֵ
                       % Z(i,2) = n;
                    %end
                    %if ceil(Z(i,1))>n      %�ų����������곬������ֵ
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
                    if ceil(Zaa(2,i))<=0      %�ų����������곬������ֵ
                        Zaa(2,i) = 1;
                    end
                    if ceil(Zaa(1,i))<=0      %�ų����������곬������ֵ
                        Zaa(1,i) = 1;
                    end
                   % if ceil(Zaa(2,i))>n      %�ų����������곬������ֵ
                       % Zaa(2,i) = n;
                   % end
                   % if ceil(Zaa(1,i))>n      %�ų����������곬������ֵ
                      %  Zaa(1,i) = n;
                   % end
     I_ad(ceil(Zaa(2,i)),ceil(Zaa(1,i))) = 1;
 end
  
  %subplot(4, 4, 13),imshow(I_ad),title("��б������ͼ��");
        
        d_adjust= I_ad;
    end 
