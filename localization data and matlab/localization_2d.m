%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code is suit for paper:
% RFID 3D Indoor localization for tag and tag-free target based on interference
% and the data :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
cc=299792458;% speed of light  
wl=cc/920625000*100;% wavelength(cm)
tag=-1000:2000;
at=42;%reference tag number (a reference tag is used in both x-axis and y-axis)
%% Ideal target coordinate:(xt，yt，zt)
xt=60;yt=60;
tpoint=0:15:300;%  interval of reference tag is 15cm
%%initial state data
%ideal distance and ideal unwrapped phase from antenna to reference tags
row_ideal_antenna_dist1=sqrt(270^2+150^2+(fliplr(tpoint)-150).^2);
row_ideal_ant1=row_ideal_antenna_dist1/wl*2*pi;
colomn_ideal_antenna_dist1=sqrt(270^2+150^2+(tpoint-150).^2);
colomn_ideal_ant1=colomn_ideal_antenna_dist1/wl*2*pi;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measured data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial state ,the data is list as: EPC, Phase,RSSI,Antenna Number
p3=importdata('D:\\实验数据\\第一篇\\2d_16点结果新\\data\\instant noodles_16points\\(60,60)initial1.txt');
[m3,n3]=size(p3);
p3_1=[];
p3_2=[];
for i=1:m3
    if p3(i,4)==1
        p3_1=[p3_1;p3(i,1:3)];
    end
    if p3(i,4)==2
        p3_2=[p3_2;p3(i,1:3)];
    end
end
phase_gan_mean3=[];
rssi_gan_mean3=[];
for j=1:41;
    phase_gan3{j}=[];  
    rssi_gan3{j}=[];
    for i=1:length(p3_1);    
        if p3_1(i,1)==j+240;
            phase_gan3{j}=[phase_gan3{j},p3_1(i,2)]; 
            rssi_gan3{j}=[rssi_gan3{j},p3_1(i,3)];
        end
    end
    if max(phase_gan3{j})-min(phase_gan3{j})>5;
        fprintf('phase jumpping, tag number:%d\n',j);
        for k=1:length(phase_gan3{j});
            if phase_gan3{j}(k)>5.5;
                phase_gan3{j}(k)=phase_gan3{j}(k)-2*pi;
            end
            fprintf('solving phase jumpping complete%d\n',j);
        end
    end    
    phase_gan_mean3=[phase_gan_mean3,mean(phase_gan3{j})]; 
    rssi_gan_mean3=[rssi_gan_mean3,mean(10^6*10.^(rssi_gan3{j}./10))];% convert the RSSI unit to mw
end
phase_gan_mean3=2*pi-phase_gan_mean3;
%% Data measured after target was presented
p4=importdata('D:\\实验数据\\第一篇\\2d_16点结果新\\data\\instant noodles_16points\\(60,60)noodles1.txt');
[m4,n4]=size(p4);
p4_1=[];
p4_2=[];
for i=1:m4
    if p4(i,4)==1
        p4_1=[p4_1;p4(i,1:3)];
    end
    if p4(i,4)==2
        p4_2=[p4_2;p4(i,1:3)];
    end
end
phase_tag_mean3=[];
rssi_tag_mean3=[];
for j=1:41;
    phase_tag3{j}=[];  
    rssi_tag3{j}=[];
    for i=1:length(p4_1);    
        if p4_1(i,1)==j+240;
            phase_tag3{j}=[phase_tag3{j},p4_1(i,2)];  
            rssi_tag3{j}=[rssi_tag3{j},p4_1(i,3)];
        end
    end
    if max(phase_tag3{j})-min(phase_tag3{j})>5;
        fprintf('phase jumpping,tag number:%d\r',j);
        for k=1:length(phase_tag3{j});
            if phase_tag3{j}(k)>5.5;
                phase_tag3{j}(k)=phase_tag3{j}(k)-2*pi;
            end
            fprintf('solving phase jumping complete%d\n',j);
        end
    end    
end

error_x=[];
error_y=[];
x_all=[];
y_all=[];
% in rare case, one reference tag sample number is less than,this is just in case
 min_tag=length(phase_tag3{1});min_tag_num=1;
    for i =1:41
        if length(phase_tag3{i}) < min_tag
            min_tag=length(phase_tag3{i});
            min_tag_num=i;
        end
    end
while min_tag<80
    min_tag=length(phase_tag3{1});
    for i =1:41
        if length(phase_tag3{i}) < min_tag
            min_tag=length(phase_tag3{i});
            min_tag_num=i;
        end
    end
    if min_tag==0
        fprintf('sample number is less than 5,stop, tag number：%d\r\n',min_tag_num);
        phase_tag3{min_tag_num}=[ phase_tag3{min_tag_num-1}];
        rssi_tag3{min_tag_num}=[ rssi_tag3{min_tag_num-1}];
        for i =1:41
            if length(phase_tag3{i}) < min_tag
                min_tag=length(phase_tag3{i});
                min_tag_num=i;
            end
        end
        break;
    end
    if min_tag<80 && min_tag>0

        fprintf('sample number is less than 5，stop,tag number: %d\r\n',min_tag_num);
        while length(phase_tag3{min_tag_num})<80
           phase_tag3{min_tag_num}=[phase_tag3{min_tag_num},phase_tag3{min_tag_num}];
           rssi_tag3{min_tag_num}=[rssi_tag3{min_tag_num},rssi_tag3{min_tag_num}];
        end
%         break;
    end
end

% use mean of 5 values as a measurement.
for ii=1:floor(min_tag/5)
phase_tag_single=[];
rssi_tag_single=[];  
for jj=1:41
    phase_tag_5(jj,:)=phase_tag3{jj}((ii-1)*5+1:(ii-1)*5+5);
    phase_tag_single=[phase_tag_single,mean(phase_tag3{jj}((ii-1)*5+1:(ii-1)*5+5))];
    rssi_tag_single=[rssi_tag_single,10^6*10^(mean(rssi_tag3{jj}((ii-1)*5+1:(ii-1)*5+5))/10)];
    phase_tag_mean3=2*pi-phase_tag_single;
    rssi_tag_mean3=rssi_tag_single;
end
% vector operation 
i=sqrt(-1);
vec_gan3=[];
vec_tag3=[];
vec_gan3=rssi_gan_mean3.*exp(i*phase_gan_mean3);
vec_tag3=rssi_tag_mean3.*exp(i*phase_tag_mean3);
vec_dif3=vec_tag3 - vec_gan3;
real_vec_dif3=real(vec_dif3);
image_vec_dif3=imag(vec_dif3);
row_vec=flipud([real_vec_dif3(1:21)',image_vec_dif3(1:21)']);
colomn_vec=[real_vec_dif3(22:41)',image_vec_dif3(22:41)'];
vec_21_colomn=rssi_tag_mean3(21).*exp(i*(phase_tag_mean3(21)+pi)) - rssi_gan_mean3(21).*exp(i*(phase_gan_mean3(21)+pi));
colomn_vec=[[real(vec_21_colomn),imag(vec_21_colomn)];colomn_vec];

dif_phase_vec3=angle(vec_dif3);
dif_abs_vec3=abs(vec_dif3);
row_ori3=phase_gan_mean3(1:at/2);
colomn_ori3=phase_gan_mean3(21:41)';

colomn_ori3(1)=colomn_ori3(1)+pi;
if colomn_ori3(1)>2*pi
    colomn_ori3(1)=colomn_ori3(1)-3*pi;
end
row_final3=dif_phase_vec3(1:at/2);
row_final3=fliplr(row_final3);
colomn_final3=dif_phase_vec3(21:41)';
colomn_final3(1)=colomn_final3(1)+pi;% this reference tag is used both in x-axis and y-axis,the rotation performence of tag 
if colomn_final3(1)>2*pi
    colomn_final3(1)=colomn_final3(1)-3*pi;
end
row_rssi3=fliplr(dif_abs_vec3(1:at/2));
colomn_rssi3=dif_abs_vec3(21:41);
% finding initial state wrapped times,row_roll_ori3 and colomn_roll_ori3
row_roll_ori3=zeros(1,at/2);
colomn_roll_ori3=zeros(21,1);
for i=(at/2+1)/2:-1:2    
    while (row_ori3(i-1)+row_roll_ori3(i-1)*2*pi) - (row_ori3(i)+row_roll_ori3(i)*2*pi)<-0.5
        row_roll_ori3(i-1)=row_roll_ori3(i-1)+1;
    end
end
for i=(at/2+1)/2:at/2-1
    while (row_ori3(i+1)+row_roll_ori3(i+1)*2*pi) - (row_ori3(i)+row_roll_ori3(i)*2*pi)<-0.5
        row_roll_ori3(i+1)=row_roll_ori3(i+1)+1;
    end
end
for i=10:-1:2
    while (colomn_ori3(i-1,1)+colomn_roll_ori3(i-1,1)*2*pi) - (colomn_ori3(i,1)+colomn_roll_ori3(i,1)*2*pi)<-0.5
        colomn_roll_ori3(i-1,1)=colomn_roll_ori3(i-1,1)+1;
    end
end
for i=10:21-1
    while (colomn_ori3(i+1,1)+colomn_roll_ori3(i+1,1)*2*pi) - (colomn_ori3(i,1)+colomn_roll_ori3(i,1)*2*pi)<-0.5
        colomn_roll_ori3(i+1,1)=colomn_roll_ori3(i+1,1)+1;
    end
end
row_ori_adj3=row_ori3+row_roll_ori3*2*pi;
colomn_ori_adj3= colomn_ori3+colomn_roll_ori3*2*pi;

%% Finding x-axis target projection point： fr3
angle_row_fir=[];
angle_row_four=[];
for i=1:length(row_vec)-1
    angle_cos_row(i)=dot(row_vec(i+1,:),row_vec(i,:))/(norm(row_vec(i+1,:))*norm(row_vec(i,:)));
    row_acos(i)=acos(angle_cos_row(i));
    
    cross_rowi=cross([row_vec(i+1,:),0],[row_vec(i,:),0]);
    row_sign(i)=cross_rowi(3)/abs(cross_rowi(3));
    angle_sin_row(i)=cross_rowi(3)/(norm(row_vec(i+1,:))*norm(row_vec(i,:)));
    row_asin(i)=asin(angle_sin_row(i));
%  
end
row_ia_sel=find(row_acos<pi/2);
%expand the range
row_ia_sel_add=row_ia_sel(1);
for i=1:length(row_ia_sel)-1
    if row_ia_sel(i+1)-row_ia_sel(i)==2
        row_ia_sel_add=[row_ia_sel_add,row_ia_sel(i+1)-1,row_ia_sel(i+1)];
    else
        row_ia_sel_add=[row_ia_sel_add,row_ia_sel(i+1)];
    end
end
var_row_seper={};
if length(row_ia_sel_add)==1
    sec_row_num=1;
     var_row_seper{1}(1)=row_ia_sel_add(1);
else
    for i=1:length(row_ia_sel_add)
        if i==1
            sec_row_num=1;
            var_row_seper{1}(1)=row_ia_sel_add(1);
        else
        if row_ia_sel_add(i)-row_ia_sel_add(i-1)==1
            var_row_seper{sec_row_num}=[var_row_seper{sec_row_num},row_ia_sel_add(i)];
        else

             var_row_seper{sec_row_num+1}=row_ia_sel_add(i);
             sec_row_num=sec_row_num+1;
        end
        end  
    end 
end
max_third=0;
max_sec=0;max_secnum=0;
max_vec_row=0;

for i=1:sec_row_num
    if i==1
        max_vec_row=length(var_row_seper{1});
        max_num=1;
    else
        if length(var_row_seper{i})>=max_vec_row

            max_third=max_sec;
            max_thirdnum=max_secnum;
            max_sec=max_vec_row;
            max_secnum=max_num; 
            max_vec_row=length(var_row_seper{i});
            max_num=i;         
        else
             if length(var_row_seper{i})>=max_sec
                max_third=max_sec;
                max_thirdnum=max_secnum;
                max_sec=length(var_row_seper{i});
                max_secnum=i;
             else
                 if length(var_row_seper{i})>=max_sec
                        max_third=length(var_row_seper{i});
                        max_thirdnum=i;
                 end
             end
         end
     end
end

%in rare large noise case,delete the point that rssi <max(rssi),these point is regarded as too far to projection point
row_rssi_min=find(row_rssi3<0.1*max(row_rssi3));
for i=length(row_ia_sel):-1:1
    if ismember(row_ia_sel(i),row_rssi_min)
        row_ia_sel(i)=[];
    end
end
%in rare large noise case,delete the point that included angle less than pi/2,but is far from projection point
%using asin of two vector
flag_zhengfu=0;%0 normal，1 all plus，2 all minus
if row_asin(row_ia_sel(1))>0
    for i=1:length(row_ia_sel)
        if row_asin(row_ia_sel(i))<0
            flag_zhengfu=0;
            break;
        elseif i==length(row_ia_sel)
            flag_zhengfu=1;
        end
    end
elseif row_asin(row_ia_sel(1))<0
    for i=1:length(row_ia_sel)
        if row_asin(row_ia_sel(i))>0
            flag_zhengfu=0;
            break;
        elseif i==length(row_ia_sel)
            flag_zhengfu=2;
        end
    end
end
fr3_pr=round(median(row_ia_sel));
if flag_zhengfu==0  %
    %search from left side 
    row_left_record=[];
    for i=1:length(row_ia_sel)
        if row_acos(row_ia_sel(i))<pi/2 && row_asin(row_ia_sel(i))<0
            row_left_record=[row_left_record,i];
        else 
            break;
        end
    end
    if ~isempty(row_left_record)
        for i=length(row_left_record):-1:1
            if abs(row_ia_sel(row_left_record(i))-fr3_pr)>5
                row_ia_sel(row_left_record(i))=[];
            end
        end
    end
    %search from right side 
    row_right_record=[];
    for i=length(row_ia_sel):-1:1
        if row_acos(row_ia_sel(i))<pi/2 && row_asin(row_ia_sel(i))>0
            row_right_record=[row_right_record,i];
        else 
            break;
        end
    end
    if ~isempty(row_right_record)
        for i=1:length(row_right_record)
            if abs(row_ia_sel(row_right_record(i))-fr3_pr)>5 
            row_ia_sel(row_right_record(i))=[];
            end
        end
    end
end
fr3_pr=round(median(row_ia_sel));
for i=length(row_ia_sel):-1:1
    if abs(row_ia_sel(i)-fr3_pr) >=8
        row_ia_sel(i)=[];
    end
end
fr3_pr=round(median(row_ia_sel));% coarse projection point
flag_median_row=0; %using to judge the phase jumping between adjacent tags near projection point
for i=max(2,fr3_pr-5):min(fr3_pr+5,at/2-1)
if (abs(row_final3(i)-row_final3(i-1))>pi && abs(row_final3(i)-row_final3(i+1))>pi) || (abs(row_final3(i)-row_final3(i-1))>pi &&abs(row_final3(i)-row_final3(i+1))<pi/4&& abs(row_final3(i+1)-row_final3(i+2))>pi)
    flag_median_row=1;
end
end
if ~flag_median_row 
    if ismember(fr3_pr,var_row_seper{max_num})
        fr3=ceil((var_row_seper{max_num}(1)+ var_row_seper{max_num}(end))/2);
    elseif max_sec>=3 && ismember(fr3_pr,var_row_seper{max_secnum})
         
            fr3=ceil((var_row_seper{max_secnum}(1)+ var_row_seper{max_secnum}(end))/2);
    else
             fr3=fr3_pr;
    end

else
             fr3=fr3_pr;
end



%% Finding y-axis target wrapped times： fc3
for i=1:length(colomn_vec)-1
    angle_cos_colomn(i)=dot(colomn_vec(i+1,:),colomn_vec(i,:))/(norm(colomn_vec(i+1,:))*norm(colomn_vec(i,:)));
    colomn_acos(i)=acos(angle_cos_colomn(i));
    cross_colomni=cross([colomn_vec(i+1,:),0],[colomn_vec(i,:),0]);
    colomn_sign(i)=cross_colomni(3)/abs(cross_colomni(3));
    colomn_asin(i)=asin(cross_colomni(3)/(norm(colomn_vec(i+1,:))*norm(colomn_vec(i,:))));
end
colomn_ia_sel=find(colomn_acos<pi/2);
colomn_ia_sel_add=colomn_ia_sel(1);
%expand the range
for i=1:length(colomn_ia_sel)-1
    if colomn_ia_sel(i+1)-colomn_ia_sel(i)==2
        colomn_ia_sel_add=[colomn_ia_sel_add,colomn_ia_sel(i+1)-1,colomn_ia_sel(i+1)];
    else
        colomn_ia_sel_add=[colomn_ia_sel_add,colomn_ia_sel(i+1)];
    end
end
var_colomn_seper={};
if length(colomn_ia_sel_add)==1
    sec_colomn_num=1;
     var_colomn_seper{1}(1)=colomn_ia_sel_add(1);
else
for i=1:length(colomn_ia_sel_add)
    if i==1
        sec_colomn_num=1;
        var_colomn_seper{1}(1)=colomn_ia_sel_add(1);
    else
    if colomn_ia_sel_add(i)-colomn_ia_sel_add(i-1)==1
        var_colomn_seper{sec_colomn_num}=[var_colomn_seper{sec_colomn_num},colomn_ia_sel_add(i)];
    else
         var_colomn_seper{sec_colomn_num+1}=colomn_ia_sel_add(i);
         sec_colomn_num=sec_colomn_num+1;
    end
    end
        
end 
end

max_third=0;
max_sec=0;max_secnum=0;
max_vec_colomn=0;
for i=1:sec_colomn_num
    if i==1
        max_vec_colomn=length(var_colomn_seper{1});
        max_num=1;
    else
    if length(var_colomn_seper{i})>=max_vec_colomn
        
        max_third=max_sec;
        max_thirdnum=max_secnum;
       
        max_sec=max_vec_colomn;
        max_secnum=max_num;
        max_vec_colomn=length(var_colomn_seper{i});
        max_num=i;         
    else
        if length(var_colomn_seper{i})>=max_sec
        max_third=max_sec;
        max_thirdnum=max_secnum;
       
        max_sec=length(var_colomn_seper{i});
        max_secnum=i;
        else
            if  length(var_colomn_seper{i})>=max_third
            max_third=length(var_colomn_seper{i});
            max_thirdnum=i;
            end
        end
    end
    end
        
end

colomn_ia_sel=find(colomn_acos<pi/2);
% delete the rssi <0.1 max(rssi) points,these point is regarded as too far to projection point
colomn_rssi_min=find(colomn_rssi3<0.1*max(colomn_rssi3));
for i=length(colomn_ia_sel):-1:1
    if ismember(colomn_ia_sel(i),colomn_rssi_min)
        colomn_ia_sel(i)=[];
    end
end

% delete the point that included angle less than pi/2,but is far from projection point
%using asin of two vector

flag_zhengfu=0;%0 normal ，1 all plus，2 all minus
if colomn_asin(colomn_ia_sel(1))>0
    for i=1:length(colomn_ia_sel)
        if colomn_asin(colomn_ia_sel(i))<0
            flag_zhengfu=0;
            break;
        elseif i==length(colomn_ia_sel)
            flag_zhengfu=1;
        end
    end
elseif colomn_asin(colomn_ia_sel(1))<0
    for i=1:length(colomn_ia_sel)
        if colomn_asin(colomn_ia_sel(i))>0
            flag_zhengfu=0;
            break;
        elseif i==length(colomn_ia_sel)
            flag_zhengfu=2;
        end
    end
end
fc3_pr=round(median(colomn_ia_sel));
if flag_zhengfu==0    
     %search from left
    colomn_left_record=[];
    for i=1:length(colomn_ia_sel)
        if colomn_acos(colomn_ia_sel(i))<pi/2 && colomn_asin(colomn_ia_sel(i))<0
            colomn_left_record=[colomn_left_record,i];
        else 
            break;
        end
    end
    if ~isempty(colomn_left_record)
        for i=length(colomn_left_record):-1:1
            if abs(colomn_ia_sel(colomn_left_record(i))-fc3_pr)>5
                colomn_ia_sel(colomn_left_record(i))=[];
            end
        end
    end
    %search from right
    colomn_right_record=[];
    for i=length(colomn_ia_sel):-1:1
        if colomn_acos(colomn_ia_sel(i))<pi/2 && colomn_asin(colomn_ia_sel(i))>0
            colomn_right_record=[colomn_right_record,i];
        else 
            break;
        end
    end
    if ~isempty(colomn_right_record)
        for i=1:length(colomn_right_record)
            if abs(colomn_ia_sel(colomn_right_record(i))-fc3_pr)>5
                colomn_ia_sel(colomn_right_record(i))=[];
            end
        end
    end
end
fc3_pr=round(median(colomn_ia_sel));
for i=length(colomn_ia_sel):-1:1
    if abs(colomn_ia_sel(i)-fc3_pr) >=8
        colomn_ia_sel(i)=[];
    end
end
fc3_pr=round(median(colomn_ia_sel));% coarse projection points
flag_median_colomn=0;
for i=max(2,fc3_pr-5):min(fc3_pr+5,at/2-1)
    if (abs(colomn_final3(i)-colomn_final3(i-1))>pi && abs(colomn_final3(i)-colomn_final3(i+1))> pi)||(abs(colomn_final3(i)-colomn_final3(i-1))>pi && abs(colomn_final3(i)-colomn_final3(i+1))<pi/4 && abs(colomn_final3(i+1)-colomn_final3(i+2))> pi) %在中间位置有相位循环
        flag_median_colomn=1;
    end
end
if ~flag_median_colomn  
   
    if ismember(fc3_pr,var_colomn_seper{max_num})
        fc3=ceil((var_colomn_seper{max_num}(1)+ var_colomn_seper{max_num}(end))/2);
    elseif max_sec>=3 && ismember(fc3_pr,var_colomn_seper{max_secnum})
        
            fc3=ceil((var_colomn_seper{max_secnum}(1)+ var_colomn_seper{max_secnum}(end))/2);
    else
        fc3=fc3_pr;
    end
else
    fc3=fc3_pr;
end

%% Calculate the relative phase wrapping times 
row_roll_final_ref3=zeros(1,at/2);
colomn_roll_final_ref3=zeros(21,1);
for i=fr3:-1:2    
    while (row_final3(i-1)+row_roll_final_ref3(i-1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi)<-0.2
        row_roll_final_ref3(i-1)=row_roll_final_ref3(i-1)+1;
    end
    if abs(i-1-fr3)<=3
        if (row_final3(i-1)+row_roll_final_ref3(i-1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi) >4
            row_roll_final_ref3(i-1)=row_roll_final_ref3(i-1)-1;
        end
    elseif abs(i-1-fr3)>=4 && abs(i-fr3)<=6
        if (row_final3(i-1)+row_roll_final_ref3(i-1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi) >5
            row_roll_final_ref3(i-1)=row_roll_final_ref3(i-1)-1;
        end
    elseif abs(i-1-fr3)>=7
        if (row_final3(i-1)+row_roll_final_ref3(i-1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi) <0.5
            row_roll_final_ref3(i-1)=row_roll_final_ref3(i-1)+1;
        end
    end
end
for i=fr3:at/2-1
   
    while (row_final3(i+1)+row_roll_final_ref3(i+1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi)<-0.2
        row_roll_final_ref3(i+1)=row_roll_final_ref3(i+1)+1;
    end
    if abs(i+1-fr3)<=3
        if (row_final3(i+1)+row_roll_final_ref3(i+1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi) >4
             row_roll_final_ref3(i+1)=row_roll_final_ref3(i+1)-1;
        end
    elseif abs(i+1-fr3)>=4 && abs(i+1-fr3)<=6
         if (row_final3(i+1)+row_roll_final_ref3(i+1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi) >5
             row_roll_final_ref3(i+1)=row_roll_final_ref3(i+1)-1;
         end
    elseif abs(i+1-fr3)>=7
        if (row_final3(i+1)+row_roll_final_ref3(i+1)*2*pi) - (row_final3(i)+row_roll_final_ref3(i)*2*pi)<0.5
             row_roll_final_ref3(i+1)=row_roll_final_ref3(i+1)+1;
        end
    end
end
for i=fc3:-1:2
    
    while (colomn_final3(i-1,1)+colomn_roll_final_ref3(i-1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi)<-0.2
        colomn_roll_final_ref3(i-1,1)=colomn_roll_final_ref3(i-1,1)+1;
    end
    if abs(i-1-fc3)<=3
        if (colomn_final3(i-1,1)+colomn_roll_final_ref3(i-1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) >4
            colomn_roll_final_ref3(i-1,1)=colomn_roll_final_ref3(i-1,1)-1;
        end
    elseif abs(i-1-fc3)>=4 && abs(i-1-fc3)<=6
        if (colomn_final3(i-1,1)+colomn_roll_final_ref3(i-1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) >5
            colomn_roll_final_ref3(i-1,1)=colomn_roll_final_ref3(i-1,1)-1;
        end
    elseif abs(i-1-fc3)>=7
        if (colomn_final3(i-1,1)+colomn_roll_final_ref3(i-1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) <0.5
            colomn_roll_final_ref3(i-1,1)=colomn_roll_final_ref3(i-1,1)+1;
        end
    end
end
for i=fc3:21-1
     
    while (colomn_final3(i+1,1)+colomn_roll_final_ref3(i+1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi)<-0.2
        colomn_roll_final_ref3(i+1,1)=colomn_roll_final_ref3(i+1,1)+1;
    end
    if abs(i+1-fc3)<=3
        if (colomn_final3(i+1,1)+colomn_roll_final_ref3(i+1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) >4
            colomn_roll_final_ref3(i+1,1)=colomn_roll_final_ref3(i+1,1)-1;
        end
    elseif abs(i+1-fc3)>=4 && abs(i+1-fc3)<=6
        if (colomn_final3(i+1,1)+colomn_roll_final_ref3(i+1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) >5
            colomn_roll_final_ref3(i+1,1)=colomn_roll_final_ref3(i+1,1)-1;
        end
    elseif abs(i+1-fc3)>=7
         if (colomn_final3(i+1,1)+colomn_roll_final_ref3(i+1,1)*2*pi) - (colomn_final3(i,1)+colomn_roll_final_ref3(i,1)*2*pi) <0.5
            colomn_roll_final_ref3(i+1,1)=colomn_roll_final_ref3(i+1,1)+1;
        end
    end
end
row_final_adj3=row_final3+row_roll_final_ref3*2*pi;
colomn_final_adj3= colomn_final3+colomn_roll_final_ref3*2*pi;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the phase difference and distance difference. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:at/2-1
    dif_row_aij1(1,i)    = (row_final_adj3(i+1) - row_final_adj3(i))  -  (row_ideal_ant1(i+1) - row_ideal_ant1(i));
end
for i=1:21-1
    dif_colomn_aij1(i,1) = (colomn_final_adj3(i+1) - colomn_final_adj3(i))  -  (colomn_ideal_ant1(i+1) - colomn_ideal_ant1(i));
end
% getting x-coordinate
xdata_row1=0:15:300-15;
ydata_row1=dif_row_aij1;
empty1=isnan(ydata_row1);
for i=length(empty1):-1:1
    if empty1(i)==1
        xdata_row1(i)=[];ydata_row1(i)=[];
    end
end
x01=[1,1];
options = optimset('MaxFunEvals',20000,'MaxIter',20000,'TolFun',1e-10,'Display','iter');

[x1,resnorm_x1,residual_row1]=lsqcurvefit('nodao',x01,xdata_row1,ydata_row1,[],[],options);
x_result=x1;
% in rare case, resnorm >3 means the point has great interference effect and delete the point and recalculation.
if max(abs(residual_row1))>3
    xdata_row1_again=xdata_row1;ydata_row1_again=ydata_row1;
    part_row1=find(abs(residual_row1)>3);
    part_row2=part_row1;
    for i=1:length(part_row1)
        if part_row1(i)+1 <= length(residual_row1) &&  part_row1(i)-1 >= 1
            if abs(residual_row1(part_row1(i)+1)) > abs(residual_row1(part_row1(i)-1)) 
                part_row2=[part_row2,part_row1(i)+1];
            else part_row2=[part_row2,part_row1(i)-1];
            end
        else if part_row1(i)==length(residual_row1)
                if abs(residual_row1(part_row1(i)-1)) > 1.5
                     part_row2=[part_row2,part_row1(i)-1];
                end
            else if part_row1(i)==1
                    if abs(residual_row1(part_row1(i)+1)) > 1.5
                         part_row2=[part_row2,part_row1(i)+1];
                    end   
                end
            end
        end
                    
    end
    part_row2=sort(unique(part_row2));
    for i=length(part_row2):-1:1
         ydata_row1_again(part_row2(i))=[];xdata_row1_again(part_row2(i))=[];
    end
    x01again=x1;
    options = optimset('MaxFunEvals',20000,'MaxIter',20000,'TolFun',1e-10,'Display','iter');
    [x1_again,resnorm_x2,residual_row1_again]=lsqcurvefit('nodao',x01again,xdata_row1_again,ydata_row1_again,[],[],options);
    x_result=x1_again;
   
else
    resnorm_x2=0;
end
if ~isreal(x_result)
    continue;
end
% getting y-coordinate
xdata_colomn1=0:15:300-15;

ydata_colomn1=dif_colomn_aij1';
empty2=isnan(ydata_colomn1);
for i=length(empty2):-1:1
    if empty2(i)==1
        xdata_colomn1(i)=[];ydata_colomn1(i)=[];%dif_colomn1_e(i)=[];
    end
end
y01=[1,1];
options = optimset('MaxFunEvals',20000,'MaxIter',20000,'TolFun',1e-10,'Display','iter');

[y1,resnorm_y1,residual_colomn1]=lsqcurvefit('nodao',y01,xdata_colomn1,ydata_colomn1,[],[],options);
y_result=y1;



if max(abs(residual_colomn1))>3
    xdata_colomn1_again=xdata_colomn1;ydata_colomn1_again=ydata_colomn1;%dif_colomn1_again=dif_colomn1_e;
    part_colomn1=find(abs(residual_colomn1)>3);
    part_colomn2=part_colomn1;
    for i=1:length(part_colomn1)
        if part_colomn1(i)+1 <= length(residual_colomn1) &&  part_colomn1(i)-1 >= 1
            if abs(residual_colomn1(part_colomn1(i)+1)) > abs(residual_colomn1(part_colomn1(i)-1)) %&& abs(residual_colomn1(part_colomn1(i)+1))>1.5
                part_colomn2=[part_colomn2,part_colomn1(i),part_colomn1(i)+1];
            else part_colomn2=[part_colomn2,part_colomn1(i),part_colomn1(i)-1];
            end
        else if part_colomn1(i)==length(residual_colomn1)
                if abs(residual_colomn1(part_colomn1(i)-1)) > 1.5
                     part_colomn2=[part_colomn2,part_colomn1(i)-1];
                end
            else if part_colomn1(i)==1
                    if abs(residual_colomn1(part_colomn1(i)+1)) > 1.5
                         part_colomn2=[part_colomn2,part_colomn1(i)+1];
                    end   
                end
            end
        end
                    
    end
    part_colomn2=sort(unique(part_colomn2));
    for i=length(part_colomn2):-1
         ydata_colomn1_again(part_colomn2(i))=[];xdata_colomn1_again(part_colomn2(i))=[];dif_colomn1_again(part_colomn2(i))=[];
    end
    y01_again=y1';
    options = optimset('MaxFunEvals',20000,'MaxIter',20000,'TolFun',1e-10,'Display','iter');

    [y1_again,resnorm_y2,residual_colomn1_again]=lsqcurvefit('nodao',y01_again,xdata_colomn1_again,ydata_colomn1_again,[],[],options);
    y_result=y1_again';
else
    resnorm_y2=0;
end
if ~isreal(y_result)
    continue;
end
if x_result(1)>300 || x_result(1)<0 || y_result(1)>300 || y_result(1)<0 || ~isreal(x_result(1)) || ~isreal(y_result(1))
    continue;
end

x_all=[x_all;x_result];
y_all=[y_all;y_result];
error_x=[error_x,x_result(1)-xt];
error_y=[error_y,y_result(1)-yt];
end


figure;plot(error_x);hold on;plot(error_y,'r');
%% Save the result
% fxyresult=fopen('D:\\实验数据\\2d_16点结果新\\可乐结果\\（240,120）基准1可乐5结果.txt','w');
% for i=1:length(error_x)
%     fprintf(fxyresult,'%f \t%f \t%f \t%f\r\n',x_all(i,1),y_all(i,1),error_x(i),error_y(i));
% end
% fclose(fxyresult);

