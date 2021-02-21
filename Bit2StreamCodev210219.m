clc;
clear;
close all; 
% Notice that the white part will be etched, while the dark part will be left.

DwellTime_Max = 100;  % 单点最大停留时间10us
xx_res_part = 80;     %最大4095 （0~4095）                                                                 % % % % 如果设置的坐标和原图坐标非等比例怎么办？
%yy_res_part = 50;
N_Process_Repeating = 20;     %扫描次数

N_Array_Row = 6 ; % 阵列重复行数
N_Array_Col = 4 ;  % 阵列重复列数
P_Array_Row = 100; % 阵列行间距
P_Array_Col = 80; % 阵列列间距

%% Step1 读取数据至矩阵
    imagedata = imread('Logo_CAS.png'); % % % 导入图片
    % imagedata = imread('logo_Color.png');
    % imdata1 = imagedata(:,:,1); imdata2 = imagedata(:,:,2); imdata3 = imagedata(:,:,3); %三色矩阵
    imdata_Gay = (1/3 *imagedata(:,:,1) + 1/3 *imagedata(:,:,2) + 1/3 *imagedata(:,:,3));  %三色混合为灰度
    imdata_Gay = im2double(imdata_Gay);
    imdata_Gay = flip(imdata_Gay);                                                                       % % % % 此操作是否有必要？
    [numRows,numCols] = size(imdata_Gay);    %读取图的横纵尺寸  [numRows/行数, numCols/行数]
    yy_res_part = round(xx_res_part./numCols.* numRows);   %等比例构建插值图
    % All_White = ones(numRows,numCols) ; 
    % imdata_Gay = All_White - imdata_Gay;  % 取反
          figure('Name','imdata_Gay','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(imdata_Gay); shading interp%伪彩色图
          print(['Step1_imdata_Gay','.tif'] ,'-dtiffn','-r300');
 
%% Step2 确定dwell time/ 停留时间 （此处我们默认单点的最大停留时间为10us）
    Dwell_Time =  floor(imdata_Gay * DwellTime_Max);
          figure('Name','Dwell time','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time); shading interp%伪彩色图
          print(['Step2_Dwell time','.tif'] ,'-dtiffn','-r300');
     
%% Step3 根据尺寸重排
  % Step3.1 从二维矩阵构建xyz散点坐标
      yy_region_Gay = linspace(1, numRows, numRows); %定义x轴坐标
      xx_region_Gay = linspace(1, numCols, numCols); %定义y轴坐标   
      yy_region_Gay = yy_region_Gay';
      for y_jj_Gay = 1:numRows
         xx_region_Gay(y_jj_Gay,:) = xx_region_Gay(1,:);
      end
      for x_ii_Gay = 1:numCols
         yy_region_Gay(:,x_ii_Gay) = yy_region_Gay(:,1);
      end
      xx_region_Gay_2 = reshape(xx_region_Gay,[numRows*numCols,1]);  %构建xyz坐标的数据y列
      yy_region_Gay_2 = reshape(yy_region_Gay,[numRows*numCols,1]);  %构建xyz坐标的数据x列
      Dwell_Time_2 = reshape(Dwell_Time,[numRows*numCols,1]);        %构建xyz坐标的数据z列
      Data_3_1(:,:,1) = xx_region_Gay_2; Data_3_1(:,:,2) = yy_region_Gay_2; Data_3_1(:,:,3) = Dwell_Time_2;
       % Dwell_Time_xyz = [yy_region_Gay_2; xx_region_Gay_2; Dwell_Time_2]   % xyz散点坐标 
  % Step3.2 构建插值用的新坐标矩阵
      xx_region_part = linspace(0 , numRows, xx_res_part+1); %定义x轴坐标     %注意此处加入了0位置的坐标以贴合矩阵   
       %xx_region_part = xx_region_part(:,1:xx_res_part+1);               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 此种插值方式会导致插值的具体位置发生偏移，如何优化？
      yy_region_part = linspace(0 , numCols, yy_res_part+1); %定义y轴坐标
       %yy_region_part = yy_region_part(:,1:yy_res_part+1);     
      yy_region_part = yy_region_part';
      for y_jj_part = 1:yy_res_part+1
         xx_region_part(y_jj_part,:) = xx_region_part(1,:);
      end
         xx_region_part = xx_region_part';
      for x_ii_part = 1:xx_res_part+1
         yy_region_part(:,x_ii_part) = yy_region_part(:,1);
      end
         yy_region_part = yy_region_part';
  % tep3.3 从图片坐标出发，对应实际加工坐标格点进行坐标插值
      Dwell_Time_gridded = griddata(Data_3_1(:,:,2),Data_3_1(:,:,1),Data_3_1(:,:,3),xx_region_part,yy_region_part);  % 插值
      xx_region_part = xx_region_part./ numRows .* xx_res_part;
      yy_region_part = yy_region_part./ numCols .* yy_res_part;
          figure('Name','Dwelltime_Grid','NumberTitle','off');  % 画图
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%伪彩色图
          print(['Step3_Dwell_Time_gridded','.tif'] ,'-dtiffn','-r300');
        
%% Step4 矩阵转换结构尺寸 （根据所要加工的阵列尺寸修改）     % xx_res_part = 80; yy_res_part = 64; 
  % Step4.1 修剪多余部分
      xx_region_part = xx_region_part(2:xx_res_part+1,2:yy_res_part+1);     
      yy_region_part = yy_region_part(2:xx_res_part+1,2:yy_res_part+1);   
      Dwell_Time_gridded = Dwell_Time_gridded(2:xx_res_part+1,2:yy_res_part+1);
      Data_4_1(:,:,1) = yy_region_part; Data_4_1(:,:,2) = xx_region_part; Data_4_1(:,:,3) = Dwell_Time_gridded;
  % Step4.2 加工顺序为首尾相接           【奇数列 A(:,1:2:end)  偶数列  A(:,2:2:end)】    B = flip(A)；   C(1:2:2n)=A   C(2:2:2n)=B;
      SizeData_4_1 = size(Data_4_1(:,:,1));      
      Dwell_Time_gridded_1 = zeros(SizeData_4_1);       
      Dwell_Time_gridded_1(:,1:2:end) = Data_4_1(:,1:2:yy_res_part,1);            %Dwell_Time_gridded_Flipped_1       
      Dwell_Time_gridded_1(:,2:2:end) = flip(Data_4_1(:,2:2:yy_res_part,1));     %Dwell_Time_gridded_UnFlipped_1
      Dwell_Time_gridded_2 = zeros(SizeData_4_1);       
      Dwell_Time_gridded_2(:,1:2:end) = Data_4_1(:,1:2:yy_res_part,2);            %Dwell_Time_gridded_Flipped_2
      Dwell_Time_gridded_2(:,2:2:end) = flip(Data_4_1(:,2:2:yy_res_part,2));     %Dwell_Time_gridded_UnFlipped_2      
      Dwell_Time_gridded_3 = zeros(SizeData_4_1);       
      Dwell_Time_gridded_3(:,1:2:end) = Data_4_1(:,1:2:yy_res_part,3);            %Dwell_Time_gridded_Flipped_3
      Dwell_Time_gridded_3(:,2:2:end) = flip(Data_4_1(:,2:2:yy_res_part,3));     %Dwell_Time_gridded_UnFlipped_3
          figure('Name','Dwelltime_TailHead','NumberTitle','off');  % 画图
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%伪彩色图
          print(['Step4_Dwell_Time_TailHead','.tif'] ,'-dtiffn','-r300');
  % Step4.2 加工顺序为首尾相接
      % 此处可选择加工顺序首位相接，或者螺旋取点
  
  
%% Step5 矩阵顺序重排      
      xx_region_part_N = reshape(Dwell_Time_gridded_2,[xx_res_part*yy_res_part,1]);  %构建xyz坐标的数据y列
      yy_region_part_N = reshape(Dwell_Time_gridded_1,[xx_res_part*yy_res_part,1]);  %构建xyz坐标的数据x列
      Dwell_Time_gridded_N = reshape(Dwell_Time_gridded_3,[xx_res_part*yy_res_part,1]); %构建xyz坐标的数据z列
      Dwell_Time_xyz = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];   % xyz散点坐标 
      
%% Step6 单元阵列化          
      if  xx_res_part > P_Array_Col
          P_Array_Col = xx_res_part; %设置color bar 显示范围
      end    
      if  yy_res_part > P_Array_Row
          P_Array_Row = xx_res_part; %设置color bar 显示范围
      end
      Data_6_1 = [];
      for y_ii_Array = 1 : N_Array_Col % 6
          for x_jj_Array = 1 : N_Array_Row %4
              Dwell_Time_xyz_ij = [];
              Dwell_Time_xyz_ij(:,1:3) = Dwell_Time_xyz ;
              Dwell_Time_xyz_ij(:,4)= ones(xx_res_part * yy_res_part,1).* (y_ii_Array - 1);
              Dwell_Time_xyz_ij(:,5)= ones(xx_res_part * yy_res_part,1).* (x_jj_Array - 1);
              TransferMatrix_S6 =[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 P_Array_Row 0 1 0;0 0 P_Array_Col 0 1]; % 构建平移矩阵
              Dwell_Time_xyz_ij = Dwell_Time_xyz_ij * TransferMatrix_S6;
              Data_6_1 = [Data_6_1 ; Dwell_Time_xyz_ij]; % 保存平移后的矩阵
          end
      end;
      Dwell_Time_gridded_N = Data_6_1(:,1); yy_region_part_N = Data_6_1(:,2); xx_region_part_N = Data_6_1(:,3);
      
          %二维散点图
          S = 5; %坐标点的大小/尺寸
          scatter(Data_6_1(:,2),Data_6_1(:,3),S,Data_6_1(:,1),'filled') %filled表示点是实心点，缺省则为空心点
          grid on
          % h = colorbar;
          set(gcf,'Position',[50 50 400 300]);
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step6_Dwell_Time_Array','.tif'] ,'-dtiffn','-r300');
          
%% Step7 按StreamFile格式输出txt文件
  
  % 此处有必要加入一个判别 非0 数值的筛，以减少加工的数据量

  % Step7.1 构建数据，补充格式
       %Data_7_1_1 = "s"; Data_6_1_2 = N_Process_Repeating; Data_6_1_3 = numel(Dwell_Time_gridded_N);
       % Data_7_1 = [Data_6_1_1; Data_6_1_2; Data_6_1_3];
      Data_7_1 = ["s"," "," "; N_Process_Repeating," "," "; numel(Dwell_Time_gridded_N)," "," "];
      Data_7_2 = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];
      Data_7_3 = (round(Data_7_2))';
       % Data_7 = [Data_7_1;round(Data_7_2)];
  % Step7.2 按标准格式输出文件   
      fid = fopen('Data_7.txt', 'w');           % Open a file for writing
      fprintf(fid, 's\n');                      % Stream File Header
      fprintf(fid, '%d\n',N_Process_Repeating); 
      fprintf(fid, '%d\n',numel(Dwell_Time_gridded_N)); 
      fprintf(fid, '%d %d %d\n', Data_7_3);     % Print values in column order    % Three values appear on each row of the file
      fclose(fid);







