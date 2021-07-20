clc;
clear;
close all; 
% Notice that the white part will be etched, while the dark part will be left.

DwellTime_Max = 50;  % �������ͣ��ʱ��10us
xx_res_part = 500;     %���4095 ��0~4095��                                                                 % % % % ������õ������ԭͼ����ǵȱ�����ô�죿
%yy_res_part = 50;
N_Process_Repeating = 50;     %ɨ�����
x_Process_Size_max = 4096 ; y_Process_Size_max = 4096 ;  % x��ΧΪ4096��y��ΧΪ4096. Step7��ʹ��
N_enlarged = 1;   %����������

N_Array_Row = 1 ; % �����ظ�����
N_Array_Col = 1;  % �����ظ�����
P_Array_Row = 120; % �����м��
P_Array_Col = 120; % �����м��

%% Step1 ��ȡ����������
    imagedata = imread('S26-Hole-D150P330.bmp'); % % % ����ͼƬ
    % imagedata = imread('logo_Color.png');
    % imdata1 = imagedata(:,:,1); imdata2 = imagedata(:,:,2); imdata3 = imagedata(:,:,3); %��ɫ����
    imdata_Gay = (1/3 *imagedata(:,:,1) + 1/3 *imagedata(:,:,2) + 1/3 *imagedata(:,:,3));  %��ɫ���Ϊ�Ҷ�
    imdata_Gay = im2double(imdata_Gay);
    imdata_Gay = flip(imdata_Gay);                                                                       % % % % �˲����Ƿ��б�Ҫ��
    [numRows,numCols] = size(imdata_Gay);    %��ȡͼ�ĺ��ݳߴ�  [numRows/����, numCols/����]
    yy_res_part = round(xx_res_part./numCols.* numRows);   %�ȱ���������ֵͼ
    
    All_White = ones(numRows,numCols) ; 
    % imdata_Gay = All_White - imdata_Gay;  % ȡ��
          figure('Name','Step1_imdata_Gay','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[0 700 400 300]);  % Figure_1 [������Ļ��ߵ����� ������Ļ�±ߵ����� �������� ��������]
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(imdata_Gay); shading interp%α��ɫͼ
          print(['Step1_imdata_Gay','.tif'] ,'-dtiffn','-r300');
 
%% Step2 ȷ��dwell time/ ͣ��ʱ�� ���˴�����Ĭ�ϵ�������ͣ��ʱ��Ϊ10us��
    Dwell_Time =  floor(imdata_Gay * DwellTime_Max);
          figure('Name','Step2_Dwell time','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[350 700 400 300]); % Figure_2
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time); shading interp %α��ɫͼ
          print(['Step2_Dwell time','.tif'] ,'-dtiffn','-r300');
     
%% Step3 ���ݳߴ�����        % ���ݱ���������
  % Step3.1 �Ӷ�ά���󹹽�xyzɢ������
      yy_region_Gay = linspace(1, numRows, numRows); %����x������
      xx_region_Gay = linspace(1, numCols, numCols); %����y������   
      yy_region_Gay = yy_region_Gay';
      for y_jj_Gay = 1:numRows
         xx_region_Gay(y_jj_Gay,:) = xx_region_Gay(1,:);
      end
      for x_ii_Gay = 1:numCols
         yy_region_Gay(:,x_ii_Gay) = yy_region_Gay(:,1);
      end
      xx_region_Gay_2 = reshape(xx_region_Gay,[numRows*numCols,1]);  %����xyz���������y��
      yy_region_Gay_2 = reshape(yy_region_Gay,[numRows*numCols,1]);  %����xyz���������x��
      Dwell_Time_2 = reshape(Dwell_Time,[numRows*numCols,1]);        %����xyz���������z��
      Data_3_1(:,:,1) = xx_region_Gay_2; Data_3_1(:,:,2) = yy_region_Gay_2; Data_3_1(:,:,3) = Dwell_Time_2;
       % Dwell_Time_xyz = [yy_region_Gay_2; xx_region_Gay_2; Dwell_Time_2]   % xyzɢ������ 
  % Step3.2 ������ֵ�õ����������
      xx_region_part = linspace(0 , numRows, xx_res_part+1); %����x������     %ע��˴�������0λ�õ����������Ͼ���   
       %xx_region_part = xx_region_part(:,1:xx_res_part+1);               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ֲ�ֵ��ʽ�ᵼ�²�ֵ�ľ���λ�÷���ƫ�ƣ�����Ż���
      yy_region_part = linspace(0 , numCols, yy_res_part+1); %����y������
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
  % tep3.3 ��ͼƬ�����������Ӧʵ�ʼӹ���������������ֵ
      Dwell_Time_gridded = griddata(Data_3_1(:,:,2),Data_3_1(:,:,1),Data_3_1(:,:,3),xx_region_part,yy_region_part);  % ��ֵ
      xx_region_part = xx_region_part./ numRows .* xx_res_part;
      yy_region_part = yy_region_part./ numCols .* yy_res_part;
          figure('Name','Step3_Dwelltime_Grid','NumberTitle','off');  % ��ͼ
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[700 700 400 300]); % Figure_3
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp %α��ɫͼ
          print(['Step3_Dwell_Time_gridded','.tif'] ,'-dtiffn','-r300');
        
%% Step4 ����ת���ṹ�ߴ� ��������Ҫ�ӹ������гߴ��޸ģ�     % xx_res_part = 80; yy_res_part = 64; 
  % Step4.1 �޼����ಿ��
  if (1)    
      xx_region_part = xx_region_part(2:xx_res_part+1,2:yy_res_part+1);     
      yy_region_part = yy_region_part(2:xx_res_part+1,2:yy_res_part+1);   
      Dwell_Time_gridded = Dwell_Time_gridded(2:xx_res_part+1,2:yy_res_part+1);
      Data_4_1(:,:,1) = yy_region_part; Data_4_1(:,:,2) = xx_region_part; Data_4_1(:,:,3) = Dwell_Time_gridded;
  end
  % Step4.2_1 �ӹ�˳��Ϊ��β���           �������� A(:,1:2:end)  ż����  A(:,2:2:end)��    B = flip(A)��   C(1:2:2n)=A   C(2:2:2n)=B;
      if (1)   % �˴���Step5 ͬʱʹ�ã�ע���޸� ��0��to��1��
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
          figure('Name','Step4_Dwelltime_TailHead','NumberTitle','off');  % ��ͼ
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[1050 700 400 300]); %Figure_4
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%α��ɫͼ
          print(['Step4_Dwell_Time_TailHead','.tif'] ,'-dtiffn','-r300');
      end
  % Step4.2_2 �ӹ�˳��Ϊ��β���
      % �˴���ѡ��ӹ�˳����λ��ӣ���������ȡ��
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
                    

                    
          while(yy_region_part)       %y�������������ȡ
              for column=1:size(yy_region_part,2)
                  Result_Data_4_X_position(1,size(Result_Data_4_X_position,2)+1) = yy_region_part(1,column);
              end
              yy_region_part(1,:)=[];
              yy_region_part = rot90(yy_region_part);
          end
          
          while(xx_region_part)       %x�������������ȡ
              for column=1:size(xx_region_part,2)
                  Result_Data_4_Y_position(1,size(Result_Data_4_Y_position,2)+1) = xx_region_part(1,column);
              end
              xx_region_part(1,:)=[];
              xx_region_part = rot90(xx_region_part);
          end
          

          
         %Dwell_Time_gridded = Dwell_Time_gridded./ 100;
          while(Dwell_Time_gridded)   %�ӹ�ʱ�����������ȡ
              for column=1:size(Dwell_Time_gridded,2)
                  Result_Dwell_Time(1,size(Result_Dwell_Time,2)+1) = Dwell_Time_gridded(1,column);
              end
              Dwell_Time_gridded(1,:)=[];
              Dwell_Time_gridded = rot90(Dwell_Time_gridded);
          end
          
         Data_Spiral_Result = [ Result_Dwell_Time' Result_Data_4_X_position' Result_Data_4_Y_position'];
         % Dwell_Time_xyz = Data_Spiral_Result;        % �������ڼӹ�
         Dwell_Time_xyz = flipud(Data_Spiral_Result);  % ��������ӹ�
      end
       
  
%% Step5 ����˳������      
if (1)      
      xx_region_part_N = reshape(Dwell_Time_gridded_2,[xx_res_part*yy_res_part,1]);  %����xyz���������y��
      yy_region_part_N = reshape(Dwell_Time_gridded_1,[xx_res_part*yy_res_part,1]);  %����xyz���������x��
      Dwell_Time_gridded_N = reshape(Dwell_Time_gridded_3,[xx_res_part*yy_res_part,1]); %����xyz���������z��
      Dwell_Time_xyz = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];   % xyzɢ������ 
end     
%% Step6 ��Ԫ���л�          
      if  xx_res_part > P_Array_Col
          P_Array_Col = xx_res_part; %����color bar ��ʾ��Χ
      end    
      if  yy_res_part > P_Array_Row
          P_Array_Row = xx_res_part; %����color bar ��ʾ��Χ
      end
      Data_6_1 = [];
      for y_ii_Array = 1 : N_Array_Col % 6
          for x_jj_Array = 1 : N_Array_Row %4
              Dwell_Time_xyz_ij = [];
              Dwell_Time_xyz_ij(:,1:3) = Dwell_Time_xyz ;
              Dwell_Time_xyz_ij(:,4)= ones(xx_res_part * yy_res_part,1).* (y_ii_Array - 1);
              Dwell_Time_xyz_ij(:,5)= ones(xx_res_part * yy_res_part,1).* (x_jj_Array - 1);
              TransferMatrix_S6 =[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 P_Array_Row 0 1 0;0 0 P_Array_Col 0 1]; % ����ƽ�ƾ���
              Dwell_Time_xyz_ij = Dwell_Time_xyz_ij * TransferMatrix_S6;
              Data_6_1 = [Data_6_1 ; Dwell_Time_xyz_ij]; % ����ƽ�ƺ�ľ���
          end
      end;
      % Dwell_Time_gridded_N = Data_6_1(:,1); yy_region_part_N = Data_6_1(:,2); xx_region_part_N = Data_6_1(:,3); % % % % % % % % % %
      
          %��άɢ��ͼ
          figure('Name','Step6_Dwell_Time_Array','NumberTitle','off');
          S = 10; %�����Ĵ�С/�ߴ�
          set(gcf,'Colormap',bone);
          scatter(Data_6_1(:,2),Data_6_1(:,3),S,Data_6_1(:,1),'filled') %filled��ʾ����ʵ�ĵ㣬ȱʡ��Ϊ���ĵ�
          grid on
          % h = colorbar;
          set(gcf,'Position',[0 300 400 300]);  %Figure_6
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step6_Dwell_Time_Array','.tif'] ,'-dtiffn','-r300');
          
%% Step7 ����xy��λ�ò�����   
      % x_Process_Size_max = 4096
      %Data_7_L1 = ones (size (Data_6_1(:,2))) * DwellTime_Max - (Data_6_1(:,1)) ;
      Data_7_L1 = (Data_6_1(:,1)) ;  % ����xyλ��ʵ�ּӹ�����ľ���
      if(1)    %ϵͳ���Ƽӹ���������ɢ���ʣ� 1 Ϊ�Զ����㣬0 Ϊ�ֶ����룻
          x_Process_Size = max (Data_6_1(:,2));  
          y_Process_Size = max (Data_6_1(:,3)); 
          N_x_Process_Size = x_Process_Size_max *(5/8) / x_Process_Size;   %��������Ļ��
          N_y_Process_Size = y_Process_Size_max *(5/8) / y_Process_Size;
          N_Process_Size = max ( [N_x_Process_Size N_y_Process_Size] );
          N_enlarged = floor (N_Process_Size);                             %����������
      end   
          Data_6_2_x = Data_6_1(:,2) * floor (N_Process_Size);
          Data_6_2_y = Data_6_1(:,3) * floor (N_Process_Size); 
      
      Data_7_L2 = (Data_6_2_x) + ones (size (Data_6_2_x)) * ((x_Process_Size_max/2) - max (Data_6_2_x)/2) ;
      Data_7_L3 = (Data_6_2_y) + ones (size (Data_6_2_y)) * ((y_Process_Size_max/2) - max (Data_6_2_y)/2) ;
      Data_7_New = [Data_7_L1 , Data_7_L2 , Data_7_L3];
          
          figure('Name','Step7_Dwell_Time_Array_ReArranged','NumberTitle','off');
          S = 10; %�����Ĵ�С/�ߴ�
          set(gcf,'Colormap',bone);
          scatter(Data_7_New(:,2),Data_7_New(:,3),S,Data_7_New(:,1),'filled') %filled��ʾ����ʵ�ĵ㣬ȱʡ��Ϊ���ĵ�
          grid on
          % h = colorbar;
          set(gcf,'Position',[350 300 400 300]);  %Figure_6
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step7_Dwell_Time_Array_ReArranged','.tif'] ,'-dtiffn','-r300');
          
%% Step8 ��StreamFile��ʽ���txt�ļ�
  
  % Step8.0 �˴�����һ���б� ��0 ��ֵ��ɸ���Լ��ټӹ���������
        % Data_8_a=[0,0,0;1,0,0;1,2,3;0,2,3];   %�б�ɾ��0����
             
        Data_8_a = Data_7_New ; 
        Data_8_b = [];               % Line selected 
        for i = 1:size(Data_8_a,1)
            b = Data_8_a(i,1); 
            if b==0                  % ɾ����һ�м�ǿ��Ϊ0������
                Data_8_b(end+1) = i;
            end
        end
        Data_8_a(Data_8_b,:) = [];
        Data_8_a_L1 = floor (Data_8_a(:,1) / max (Data_8_a(:,1)) * DwellTime_Max);
        Data_8_a = [ Data_8_a_L1 , Data_8_a(:,2) , Data_8_a(:,3) ];
        % disp(Data_7_a)
        
      Dwell_Time_gridded_N = Data_8_a(:,1); yy_region_part_N = Data_8_a(:,2); xx_region_part_N = Data_8_a(:,3);  %Line 142
        
 if(1) 
  % Step8.1 �������ݣ������ʽ

      Data_8_1 = ["s"," "," "; N_Process_Repeating," "," "; numel(Dwell_Time_gridded_N)," "," "];
      Data_8_2 = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];
      Data_8_3 = (round(Data_8_2));
      %%Data_8_3 = (round(Data_8_2))';
       % Data_8 = [Data_8_1;round(Data_8_2)];
  % Step8.2 ����׼��ʽ����ļ�   
      out2=['Data_8_2.txt'];
      dlmwrite(out2,Data_8_3,'delimiter',' ','newline','pc');
      fid = fopen('Data_8_1.txt', 'w');           % Open a file for writing
      fprintf(fid, 's\n');                      % Stream File Header
      fprintf(fid, '%d\n',N_Process_Repeating); 
      fprintf(fid, '%d\n',numel(Dwell_Time_gridded_N)); 
      fclose(fid); 
      
      %�������ں������ļ�
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







