clc;
clear;
close all; 
% Notice that the white part will be etched, while the dark part will be left.

DwellTime_Max = 50;  % 单点最大停留时间10us
xx_res_part = 500;     %最大4095 （0~4095）                                                                 % % % % 如果设置的坐标和原图坐标非等比例怎么办？
%yy_res_part = 50;
N_Process_Repeating = 50;     %扫描次数
x_Process_Size_max = 4096 ; y_Process_Size_max = 4096 ;  % x范围为4096，y范围为4096. Step7中使用
N_enlarged = 1;   %坐标扩大倍数

N_Array_Row = 1 ; % 阵列重复行数
N_Array_Col = 1;  % 阵列重复列数
P_Array_Row = 120; % 阵列行间距
P_Array_Col = 120; % 阵列列间距

%% Step1 读取数据至矩阵
    imagedata = imread('S26-Hole-D150P330.bmp'); % % % 导入图片
    % imagedata = imread('logo_Color.png');
    % imdata1 = imagedata(:,:,1); imdata2 = imagedata(:,:,2); imdata3 = imagedata(:,:,3); %三色矩阵
    imdata_Gay = (1/3 *imagedata(:,:,1) + 1/3 *imagedata(:,:,2) + 1/3 *imagedata(:,:,3));  %三色混合为灰度
    imdata_Gay = im2double(imdata_Gay);
    imdata_Gay = flip(imdata_Gay);                                                                       % % % % 此操作是否有必要？
    [numRows,numCols] = size(imdata_Gay);    %读取图的横纵尺寸  [numRows/行数, numCols/行数]
    yy_res_part = round(xx_res_part./numCols.* numRows);   %等比例构建插值图
    
    All_White = ones(numRows,numCols) ; 
    % imdata_Gay = All_White - imdata_Gay;  % 取反
          figure('Name','Step1_imdata_Gay','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[0 700 400 300]);  % Figure_1 [距离屏幕左边的像素 距离屏幕下边的像素 横向像素 纵向像素]
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(imdata_Gay); shading interp%伪彩色图
          print(['Step1_imdata_Gay','.tif'] ,'-dtiffn','-r300');
 
%% Step2 确定dwell time/ 停留时间 （此处我们默认单点的最大停留时间为10us）
    Dwell_Time =  floor(imdata_Gay * DwellTime_Max);
          figure('Name','Step2_Dwell time','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[350 700 400 300]); % Figure_2
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time); shading interp %伪彩色图
          print(['Step2_Dwell time','.tif'] ,'-dtiffn','-r300');
     
%% Step3 根据尺寸重排        % 横纵比例问题需
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
          figure('Name','Step3_Dwelltime_Grid','NumberTitle','off');  % 画图
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[700 700 400 300]); % Figure_3
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp %伪彩色图
          print(['Step3_Dwell_Time_gridded','.tif'] ,'-dtiffn','-r300');
        
%% Step4 矩阵转换结构尺寸 （根据所要加工的阵列尺寸修改）     % xx_res_part = 80; yy_res_part = 64; 
  % Step4.1 修剪多余部分
  if (1)    
      xx_region_part = xx_region_part(2:xx_res_part+1,2:yy_res_part+1);     
      yy_region_part = yy_region_part(2:xx_res_part+1,2:yy_res_part+1);   
      Dwell_Time_gridded = Dwell_Time_gridded(2:xx_res_part+1,2:yy_res_part+1);
      Data_4_1(:,:,1) = yy_region_part; Data_4_1(:,:,2) = xx_region_part; Data_4_1(:,:,3) = Dwell_Time_gridded;
  end
  % Step4.2_1 加工顺序为首尾相接           【奇数列 A(:,1:2:end)  偶数列  A(:,2:2:end)】    B = flip(A)；   C(1:2:2n)=A   C(2:2:2n)=B;
      if (1)   % 此处与Step5 同时使用，注意修改 （0）to（1）
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
          figure('Name','Step4_Dwelltime_TailHead','NumberTitle','off');  % 画图
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[1050 700 400 300]); %Figure_4
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%伪彩色图
          print(['Step4_Dwell_Time_TailHead','.tif'] ,'-dtiffn','-r300');
      end
  % Step4.2_2 加工顺序为首尾相接
      % 此处可选择加工顺序首位相接，或者螺旋取点
      if (0)
          %  Data_4_1(:,:,1) = yy_region_part; Data_4_1(:,:,2) = xx_region_part; Data_4_1(:,:,3) = Dwell_Time_gridded;
          Result_Dwell_Time = [];
          
          Result_Data_4_X_position = [];
          Result_Data_4_Y_position = [];
%           Res_X_Line = linspace(1,Line_X,Line_X);
%           Res_Y_Line = linspace(1,Line_Y,Line_Y);
%           % b=repmat(a,10,1)
%           Res_X_position = repmat(Res_X_Line,Line_Y,1);
%           Res_Y_position = (repmat(Res_Y_Line,Line_X,1))';
                    

                    
          while(yy_region_part)       %y坐标矩阵螺旋读取
              for column=1:size(yy_region_part,2)
                  Result_Data_4_X_position(1,size(Result_Data_4_X_position,2)+1) = yy_region_part(1,column);
              end
              yy_region_part(1,:)=[];
              yy_region_part = rot90(yy_region_part);
          end
          
          while(xx_region_part)       %x坐标矩阵螺旋读取
              for column=1:size(xx_region_part,2)
                  Result_Data_4_Y_position(1,size(Result_Data_4_Y_position,2)+1) = xx_region_part(1,column);
              end
              xx_region_part(1,:)=[];
              xx_region_part = rot90(xx_region_part);
          end
          

          
         %Dwell_Time_gridded = Dwell_Time_gridded./ 100;
          while(Dwell_Time_gridded)   %加工时间矩阵螺旋读取
              for column=1:size(Dwell_Time_gridded,2)
                  Result_Dwell_Time(1,size(Result_Dwell_Time,2)+1) = Dwell_Time_gridded(1,column);
              end
              Dwell_Time_gridded(1,:)=[];
              Dwell_Time_gridded = rot90(Dwell_Time_gridded);
          end
          
         Data_Spiral_Result = [ Result_Dwell_Time' Result_Data_4_X_position' Result_Data_4_Y_position'];
         % Dwell_Time_xyz = Data_Spiral_Result;        % 由外向内加工
         Dwell_Time_xyz = flipud(Data_Spiral_Result);  % 由内向外加工
      end
       
  
%% Step5 矩阵顺序重排      
if (1)      
      xx_region_part_N = reshape(Dwell_Time_gridded_2,[xx_res_part*yy_res_part,1]);  %构建xyz坐标的数据y列
      yy_region_part_N = reshape(Dwell_Time_gridded_1,[xx_res_part*yy_res_part,1]);  %构建xyz坐标的数据x列
      Dwell_Time_gridded_N = reshape(Dwell_Time_gridded_3,[xx_res_part*yy_res_part,1]); %构建xyz坐标的数据z列
      Dwell_Time_xyz = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];   % xyz散点坐标 
end     
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
      % Dwell_Time_gridded_N = Data_6_1(:,1); yy_region_part_N = Data_6_1(:,2); xx_region_part_N = Data_6_1(:,3); % % % % % % % % % %
      
          %二维散点图
          figure('Name','Step6_Dwell_Time_Array','NumberTitle','off');
          S = 10; %坐标点的大小/尺寸
          set(gcf,'Colormap',bone);
          scatter(Data_6_1(:,2),Data_6_1(:,3),S,Data_6_1(:,1),'filled') %filled表示点是实心点，缺省则为空心点
          grid on
          % h = colorbar;
          set(gcf,'Position',[0 300 400 300]);  %Figure_6
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step6_Dwell_Time_Array','.tif'] ,'-dtiffn','-r300');
          
%% Step7 补偿xy轴位置并居中   
      % x_Process_Size_max = 4096
      %Data_7_L1 = ones (size (Data_6_1(:,2))) * DwellTime_Max - (Data_6_1(:,1)) ;
      Data_7_L1 = (Data_6_1(:,1)) ;  % 补偿xy位置实现加工区域的居中
      if(1)    %系统控制加工坐标点的扩散倍率， 1 为自动计算，0 为手动输入；
          x_Process_Size = max (Data_6_1(:,2));  
          y_Process_Size = max (Data_6_1(:,3)); 
          N_x_Process_Size = x_Process_Size_max *(5/8) / x_Process_Size;   %扩大至屏幕中
          N_y_Process_Size = y_Process_Size_max *(5/8) / y_Process_Size;
          N_Process_Size = max ( [N_x_Process_Size N_y_Process_Size] );
          N_enlarged = floor (N_Process_Size);                             %坐标扩大倍数
      end   
          Data_6_2_x = Data_6_1(:,2) * floor (N_Process_Size);
          Data_6_2_y = Data_6_1(:,3) * floor (N_Process_Size); 
      
      Data_7_L2 = (Data_6_2_x) + ones (size (Data_6_2_x)) * ((x_Process_Size_max/2) - max (Data_6_2_x)/2) ;
      Data_7_L3 = (Data_6_2_y) + ones (size (Data_6_2_y)) * ((y_Process_Size_max/2) - max (Data_6_2_y)/2) ;
      Data_7_New = [Data_7_L1 , Data_7_L2 , Data_7_L3];
          
          figure('Name','Step7_Dwell_Time_Array_ReArranged','NumberTitle','off');
          S = 10; %坐标点的大小/尺寸
          set(gcf,'Colormap',bone);
          scatter(Data_7_New(:,2),Data_7_New(:,3),S,Data_7_New(:,1),'filled') %filled表示点是实心点，缺省则为空心点
          grid on
          % h = colorbar;
          set(gcf,'Position',[350 300 400 300]);  %Figure_6
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step7_Dwell_Time_Array_ReArranged','.tif'] ,'-dtiffn','-r300');
          
%% Step8 按StreamFile格式输出txt文件
  
  % Step8.0 此处加入一个判别 非0 数值的筛，以减少加工的数据量
        % Data_8_a=[0,0,0;1,0,0;1,2,3;0,2,3];   %判别删除0的行
             
        Data_8_a = Data_7_New ; 
        Data_8_b = [];               % Line selected 
        for i = 1:size(Data_8_a,1)
            b = Data_8_a(i,1); 
            if b==0                  % 删除第一列即强度为0的数据
                Data_8_b(end+1) = i;
            end
        end
        Data_8_a(Data_8_b,:) = [];
        Data_8_a_L1 = floor (Data_8_a(:,1) / max (Data_8_a(:,1)) * DwellTime_Max);
        Data_8_a = [ Data_8_a_L1 , Data_8_a(:,2) , Data_8_a(:,3) ];
        % disp(Data_7_a)
        
      Dwell_Time_gridded_N = Data_8_a(:,1); yy_region_part_N = Data_8_a(:,2); xx_region_part_N = Data_8_a(:,3);  %Line 142
        
 if(1) 
  % Step8.1 构建数据，补充格式

      Data_8_1 = ["s"," "," "; N_Process_Repeating," "," "; numel(Dwell_Time_gridded_N)," "," "];
      Data_8_2 = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];
      Data_8_3 = (round(Data_8_2));
      %%Data_8_3 = (round(Data_8_2))';
       % Data_8 = [Data_8_1;round(Data_8_2)];
  % Step8.2 按标准格式输出文件   
      out2=['Data_8_2.txt'];
      dlmwrite(out2,Data_8_3,'delimiter',' ','newline','pc');
      fid = fopen('Data_8_1.txt', 'w');           % Open a file for writing
      fprintf(fid, 's\n');                      % Stream File Header
      fprintf(fid, '%d\n',N_Process_Repeating); 
      fprintf(fid, '%d\n',numel(Dwell_Time_gridded_N)); 
      fclose(fid); 
      
      %下面是融合两个文件
      fidA=fopen('Data_8_1.txt','r');
      fidB=fopen('Data_8_2.txt','r');
      DataA=fread(fidA);
      DataB=fread(fidB);
      fidC=fopen('Data_8.txt','w');
      fwrite(fidC,DataA);
      fwrite(fidC,DataB);
      fclose(fidA);
      fclose(fidB);
      fclose(fidC);
      
end







