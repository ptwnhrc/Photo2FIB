clc;
clear;
close all; 
% Notice that the white part will be etched, while the dark part will be left.

DwellTime_Max = 100;  % �������ͣ��ʱ��10us
xx_res_part = 80;     %���4095 ��0~4095��                                                                 % % % % ������õ������ԭͼ����ǵȱ�����ô�죿
%yy_res_part = 50;
N_Process_Repeating = 20;     %ɨ�����

N_Array_Row = 6 ; % �����ظ�����
N_Array_Col = 4 ;  % �����ظ�����
P_Array_Row = 100; % �����м��
P_Array_Col = 80; % �����м��

%% Step1 ��ȡ����������
    imagedata = imread('Logo_CAS.png'); % % % ����ͼƬ
    % imagedata = imread('logo_Color.png');
    % imdata1 = imagedata(:,:,1); imdata2 = imagedata(:,:,2); imdata3 = imagedata(:,:,3); %��ɫ����
    imdata_Gay = (1/3 *imagedata(:,:,1) + 1/3 *imagedata(:,:,2) + 1/3 *imagedata(:,:,3));  %��ɫ���Ϊ�Ҷ�
    imdata_Gay = im2double(imdata_Gay);
    imdata_Gay = flip(imdata_Gay);                                                                       % % % % �˲����Ƿ��б�Ҫ��
    [numRows,numCols] = size(imdata_Gay);    %��ȡͼ�ĺ��ݳߴ�  [numRows/����, numCols/����]
    yy_res_part = round(xx_res_part./numCols.* numRows);   %�ȱ���������ֵͼ
    % All_White = ones(numRows,numCols) ; 
    % imdata_Gay = All_White - imdata_Gay;  % ȡ��
          figure('Name','imdata_Gay','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(imdata_Gay); shading interp%α��ɫͼ
          print(['Step1_imdata_Gay','.tif'] ,'-dtiffn','-r300');
 
%% Step2 ȷ��dwell time/ ͣ��ʱ�� ���˴�����Ĭ�ϵ�������ͣ��ʱ��Ϊ10us��
    Dwell_Time =  floor(imdata_Gay * DwellTime_Max);
          figure('Name','Dwell time','NumberTitle','off');
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time); shading interp%α��ɫͼ
          print(['Step2_Dwell time','.tif'] ,'-dtiffn','-r300');
     
%% Step3 ���ݳߴ�����
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
          figure('Name','Dwelltime_Grid','NumberTitle','off');  % ��ͼ
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%α��ɫͼ
          print(['Step3_Dwell_Time_gridded','.tif'] ,'-dtiffn','-r300');
        
%% Step4 ����ת���ṹ�ߴ� ��������Ҫ�ӹ������гߴ��޸ģ�     % xx_res_part = 80; yy_res_part = 64; 
  % Step4.1 �޼����ಿ��
      xx_region_part = xx_region_part(2:xx_res_part+1,2:yy_res_part+1);     
      yy_region_part = yy_region_part(2:xx_res_part+1,2:yy_res_part+1);   
      Dwell_Time_gridded = Dwell_Time_gridded(2:xx_res_part+1,2:yy_res_part+1);
      Data_4_1(:,:,1) = yy_region_part; Data_4_1(:,:,2) = xx_region_part; Data_4_1(:,:,3) = Dwell_Time_gridded;
  % Step4.2 �ӹ�˳��Ϊ��β���           �������� A(:,1:2:end)  ż����  A(:,2:2:end)��    B = flip(A)��   C(1:2:2n)=A   C(2:2:2n)=B;
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
          figure('Name','Dwelltime_TailHead','NumberTitle','off');  % ��ͼ
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[100 100 400 300]);
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          pcolor(Dwell_Time_gridded); shading interp%α��ɫͼ
          print(['Step4_Dwell_Time_TailHead','.tif'] ,'-dtiffn','-r300');
  % Step4.2 �ӹ�˳��Ϊ��β���
      % �˴���ѡ��ӹ�˳����λ��ӣ���������ȡ��
  
  
%% Step5 ����˳������      
      xx_region_part_N = reshape(Dwell_Time_gridded_2,[xx_res_part*yy_res_part,1]);  %����xyz���������y��
      yy_region_part_N = reshape(Dwell_Time_gridded_1,[xx_res_part*yy_res_part,1]);  %����xyz���������x��
      Dwell_Time_gridded_N = reshape(Dwell_Time_gridded_3,[xx_res_part*yy_res_part,1]); %����xyz���������z��
      Dwell_Time_xyz = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];   % xyzɢ������ 
      
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
      Dwell_Time_gridded_N = Data_6_1(:,1); yy_region_part_N = Data_6_1(:,2); xx_region_part_N = Data_6_1(:,3);
      
          %��άɢ��ͼ
          S = 5; %�����Ĵ�С/�ߴ�
          scatter(Data_6_1(:,2),Data_6_1(:,3),S,Data_6_1(:,1),'filled') %filled��ʾ����ʵ�ĵ㣬ȱʡ��Ϊ���ĵ�
          grid on
          % h = colorbar;
          set(gcf,'Position',[50 50 400 300]);
          %set(gca,'Position',[.18 .17 .6 .67]);
          %colorbar('position',[0.82 0.17 0.03 0.73]); 
          print(['Step6_Dwell_Time_Array','.tif'] ,'-dtiffn','-r300');
          
%% Step7 ��StreamFile��ʽ���txt�ļ�
  
  % �˴��б�Ҫ����һ���б� ��0 ��ֵ��ɸ���Լ��ټӹ���������

  % Step7.1 �������ݣ������ʽ
       %Data_7_1_1 = "s"; Data_6_1_2 = N_Process_Repeating; Data_6_1_3 = numel(Dwell_Time_gridded_N);
       % Data_7_1 = [Data_6_1_1; Data_6_1_2; Data_6_1_3];
      Data_7_1 = ["s"," "," "; N_Process_Repeating," "," "; numel(Dwell_Time_gridded_N)," "," "];
      Data_7_2 = [Dwell_Time_gridded_N, yy_region_part_N, xx_region_part_N];
      Data_7_3 = (round(Data_7_2))';
       % Data_7 = [Data_7_1;round(Data_7_2)];
  % Step7.2 ����׼��ʽ����ļ�   
      fid = fopen('Data_7.txt', 'w');           % Open a file for writing
      fprintf(fid, 's\n');                      % Stream File Header
      fprintf(fid, '%d\n',N_Process_Repeating); 
      fprintf(fid, '%d\n',numel(Dwell_Time_gridded_N)); 
      fprintf(fid, '%d %d %d\n', Data_7_3);     % Print values in column order    % Three values appear on each row of the file
      fclose(fid);







