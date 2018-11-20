function left_bottom_r=left_bottom_main(I,r,meanx,meany)
%% 左侧半圆环积分，步长为5 
    [m,n]=size(I);
    img_incline=zeros(m,n);
     for i=meany:m
         for j=1:meanx  % meanx:n %1:n  
               if atan((i-meany)/(j-meanx))*180/pi >-30 && atan((i-meany)/(j-meanx))*180/pi<-20
                 img_incline(i,j)=I(i,j);
               end
         end
     end
     %% 显示部分
     %{
     sun=1;
     for i=3:m-1
         for j=3:n-1   % meanx:n %1:n  
               if img_incline(i,j)==0 && img_incline(i+1,j)~=0 || img_incline(i+1,j)==0 && img_incline(i,j)~=0
                    img_incline(i,j)=1; 
               end
             if img_incline(i,j)==0 && img_incline(i+1,j)~=0 && sun<19000  || img_incline(i+1,j)==0 && img_incline(i,j)~=0 && sun<19000
                    img_incline(i-1,j)=1; 
                    sun=sun+1;
             end
             if img_incline(i,j)==0 && img_incline(i+1,j)~=0 && sun<19000 || img_incline(i+1,j)==0 && img_incline(i,j)~=0 && sun<19000
                    img_incline(i-1,j)=1;  img_incline(i-1,j-1)=1; img_incline(i-2,j-2)=1; img_incline(i,j-1)=1;img_incline(i,j)=1; img_incline(i-1,j-2)=1;img_incline(i-2,j-1)=1;img_incline(i,j-2)=1;
                    sun=sun+1;
              end
         end
     end
     for i=1:m
         for j=1:meanx
             if i==450  || i==449 ||  i==448 ||i==447 ||i==446 
                img_incline(i-1,j)=1; 
             end
         end
     end
%}
%       figure,imshow(img_incline);   
     %%left_r
    I_sum(1)=0; step=1;
    for k=1:5:round(n-meanx) % round((m-meany)/sin(2*pi/9))
       I_sum=0; number=0; I_test=ones(m,n);
        for i=meany:m
            for j= 1:meanx %1:n %meanx:n    1:n
                if  sqrt((i-meany)^2+(j-meanx)^2)<k+5 && k<=sqrt((i-meany)^2+(j-meanx)^2)  &&  atan((i-meany)/(j-meanx))*180/pi>-30 && atan((i-meany)/(j-meanx))*180/pi<-20
                      I_sum=I_sum+img_incline(i,j);
                      number=number+1;
                end
            end
        end
         I_SUM(k)= I_sum/number ;  % 单位面积灰度值
    end
    n=length(I_SUM);
    x=1:5:n;
    v=I_SUM(x);
    xq=1:1:n;
    vq2 = interp1(x,v,xq,'linear');
%     figure,plot(xq,vq2);
%     xlim([1 n]); 
    % title('linear Interpolation'); 
    
    %%
    
    for i=round(1.3*r): length(vq2)
        data(i)=vq2(i);
    end
%    figure,plot(data);
    %%
    IndMin=find(diff(sign(diff(data)))>0)+1;%,Find local minima， data(IndMin)对应峰值，IndMin对应峰值在x轴的位置
%      figure; hold on; box on;  
%     plot(1:length(data),data);  
%     plot(IndMin,data(IndMin),'r^');
%     legend('曲线','波谷点');
%     title('计算离散节点的波谷信息', 'FontWeight', 'Bold');
    %%
    [pks_max,locs]=findpeaks(data,'minpeakdistance',1);% Find local maxima，pks_max 对应峰值，locs 对应峰值的横坐标位置
%     figure; hold on; box on;  
%     plot(1:length(data),data);  
%     plot(locs,data(locs),'r^');
%     legend('曲线','波峰点');
%     title('计算离散节点的波峰信息', 'FontWeight', 'Bold');
    [~ ,y1]=max(pks_max); % y1极大值在极值中的位置（有限个从1开始的正整数），locs(y1)是极大值对应的横坐标
    if data(locs(y1))<data(length(data))
        locs(y1)=length(data);
    end
    
    
    
    %%
    trough=zeros(1,length(IndMin)-1);
    for i=2:length(IndMin)
        trough(i-1)=IndMin(i); %IndMin对应峰值在x轴的位置，trough去除了第一个点
    end
    [min_trough,loc_trough]=min(data(trough)); %几个波谷点比较最小值，min_trough是极小值，loc_trough极小值在极值中的位置（有限个从1开始的正整数）
                                               % IndMin(loc_trough+1)是极小值在原波形中的位置(横坐标)，data(IndMin(loc_trough+1))是极小值对应的大小（纵坐标）
                                               % min_trough 与 data(IndMin(loc_trough+1))一样大小
%      figure; hold on; box on; 
%      plot(1:length(data),data);  
%     plot(IndMin(loc_trough+1),data(IndMin(loc_trough+1)),'r^') ;
    if locs(y1)<IndMin(loc_trough+1) % 保证locs(y1)必须在极小值点右侧
        locs(y1)=length(vq2);
    end
    %{
   if loc_trough==length(IndMin)-1  %如果该点是最后一个点,那么就从最后一个点开始到最后找相邻两个点梯度最大的地方就是边界
        for i=IndMin(loc_trough+1): locs(y1)-1% length(data)-1  % 根据波峰的位置，确定边界范围
            precise_point(i)=data(i+1)-data(i);
            [maxr,index]=max(precise_point);
        end
   else    % precise_point=zeros(1,locs(max_num+1)-1-locs(max_num)); %如果该点是不是最后一个点,那么就从该点到下一个波谷点找相邻两个点梯度最大的地方就是边界
        for i=IndMin(loc_trough+1): locs(y1)-1% IndMin(loc_trough+2)  % 根据波峰的位置，确定边界范围
            precise_point(i)=data(i+1)-data(i);
            [maxr,index]=max(precise_point);
        end
    
   end
    %}
    for i=IndMin(loc_trough+1): locs(y1)-1% length(data)-1  % 根据波峰的位置，确定边界范围
            precise_point(i)=data(i+1)-data(i);
            [maxr,index]=max(precise_point);
    end
%     figure; hold on; box on;  
%     plot(xq,vq2);
%     plot(index,data(index),'r^');
    left_bottom_r=index+5;