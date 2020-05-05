function characters = LicPlateRec(character_image)
% µ÷ÓÃº¯Êýrecognise Ê¶±ð×Ö·û£¬·µ»Ø³µÅÆ×Ö·û´®
%% ºº×ÖÊ¶±ð
characters = '';
lib1 = '¾©½ò¼½½úÃÉÁÉ¼ªºÚ»¦ËÕÕãÍîÃö¸ÓÂ³Ô¥¶õÏæÔÁ¹ðÇíÓå´¨¹óÔÆ²ØÉÂ¸ÊÇàÄþÐÂ';
lib2 = '1234567890ABCDEFGHJKLMNPQRSTUVWXYZ';

temp_char = character_image{1};
temp_char = temp_char(:)';
load('hanzi_net1.mat');
load('hanzi_net2.mat');
characters = strcat(characters, hanzishibie(hanzi_theta1, hanzi_theta2, temp_char, lib1));
if characters=='¸Ó'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(1:10,13:16)))<=4
        characters='´¨';
    end
end
if characters=='ÔÁ'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(30:32,5:11)))<=5
        characters='´¨';
    end
end
if characters=='½ú'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(1:2,1:16)))<=3
        characters='¹ó';
    end
end
if characters=='¹ó'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(:,6)))<=1
        characters='´¨';
    end
end
if characters=='Õã'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(1:22,6)))>=3
        characters='ÐÂ';
    end
end
if characters=='Íî'
    temp_temp_char=character_image{1};
    if sum(sum(temp_temp_char(20:26,1:4)))<=1
        characters='´¨';
    end
end
% if characters=='ËÕ'
%     temp_temp_char=character_image{1};
%     if sum(sum(temp_temp_char(11:17,1:3)))>=1
%         characters='Ô¥';
%     end
% end
switch characters
    case '¹ð'
        characters='´¨';
    case 'ºÚ'
        characters='´¨';
    case '¸Ê'
        characters='´¨';
    case 'Äþ'
        characters='´¨';
    case 'Â³'
        characters='¶õ';
    case 'Çà'
        characters='ËÕ';
        temp_temp_char=character_image{1};
        if sum(sum(temp_temp_char(11:17,1:3)))>=1
            characters='Ô¥';
        end
    case '½ò'
        characters='Ô¥';
    case 'ÃÉ'
        characters='¸Ó';
    case 'Çí'
        characters='ÔÆ';
    case '¼ª'
        characters='Íî';
    case 'ÁÉ'
        characters='´¨';
    otherwise
            characters=characters;
end
%% Êý×ÖºÍ×ÖÄ¸Ê¶±ð
finalout=[];
%finalout(1)=characters;

for i=2:7
    bw_4032=character_image{i};
        for cnt= 1:4 %40/4=10
            for cnt2=1:4 %32/4=8
                Atemp=sum(bw_4032(((cnt*10-9):(cnt*10)),((cnt2*8-7):(cnt2*8))));
                lett((cnt-1)*4+cnt2)=sum(Atemp);
            end
        end
        lett=((100-lett)/100);
        %lett=lett/320;
    lett=lett';%ÐÐÏòÁ¿±äÁÐÏòÁ¿
    feature(:,i)=lett;
end



load zimu.mat; 
a2=sim(net,feature(:,2));
load zimuandshuzi.mat; 
for i=3:7
    a3(:,i)=sim(net,feature(:,i));
end

%a1=sim('hanzi',feature(:,1));
%a2=sim('zimu',feature(:,2));
%a3=sim('zimuandshuzi',feature(:,3:7));

%out(1)=find(a1(:,1)==max(a1(:,1)));
out(2)=find(a2(:,1)==max(a2(:,1)));
for i=1:5
    [f out(i+2)]=max(a3(:,i+2));%==max(a3(:,i+2)));
end

switch out(2)
    case 1        
        finalout(2)='A';
    case 2        
        finalout(2)='B';
    case 3        
        finalout(2)='C';
    case 4        
        finalout(2)='D';    
    case 5        
        finalout(2)='E'; 
    case 6        
        finalout(2)='F';
    case 7        
        finalout(2)='G';
    case 8        
        finalout(2)='H';    
    case 9        
        finalout(2)='J';
    case 10        
        finalout(2)='K';
    case 11       
        finalout(2)='L';
    temp_temp_char=character_image{2};
    if sum(sum(temp_temp_char(1:20,16:32)))>=5
        finalout(2)='A';
    end
    case 12       
        finalout(2)='M';
    case 13       
        finalout(2)='N';
    case 14       
        finalout(2)='P';
    case 15       
        finalout(2)='Q';
    case 16       
        finalout(2)='R';
    case 17       
        finalout(2)='S';
    case 18       
        finalout(2)='T';
    case 19       
        finalout(2)='U';
    case 20       
        finalout(2)='V';    
    case 21       
        finalout(2)='W'; 
    case 22       
        finalout(2)='X';
    case 23       
        finalout(2)='Y';
    case 24       
        finalout(2)='Z';    
end
for i=3:7
    switch out(i)
        case 1        
            finalout(i)='0';
        case 2       
            finalout(i)='1';
        case 3        
            finalout(i)='2';
        case 4        
            finalout(i)='3';    
        case 5        
            finalout(i)='4'; 
        case 6        
            finalout(i)='5';
        case 7        
            finalout(i)='6';
        case 8        
            finalout(i)='7';    
        case 9        
            finalout(i)='8';
        case 10        
            finalout(i)='9';
        case 11       
            finalout(i)='A';
        case 12       
            finalout(i)='B';
        case 13       
            finalout(i)='C';
        case 14       
            finalout(i)='D';
        case 15       
            finalout(i)='E';
        case 16       
            finalout(i)='F';
        case 17       
            finalout(i)='G';
        case 18       
            finalout(i)='H';
        case 19       
            finalout(i)='J';
        case 20       
            finalout(i)='K';    
        case 21       
            finalout(i)='L'; 
        case 22       
            finalout(i)='M';
        case 23       
            finalout(i)='N';
        case 24       
            finalout(i)='P';    
        case 25       
            finalout(i)='Q';
        case 26       
            finalout(i)='R';
        case 27       
            finalout(i)='S';
        case 28       
            finalout(i)='T';
        case 29       
            finalout(i)='U';
        case 30       
            finalout(i)='V';
        case 31       
            finalout(i)='W';
        case 32      
            finalout(i)='X';
        case 33       
            finalout(i)='Y';
        case 34       
            finalout(i)='Z';
    end
end
characters=[characters finalout(2:7)]; 
