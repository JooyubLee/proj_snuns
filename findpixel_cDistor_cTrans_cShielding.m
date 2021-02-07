%% Distortion problem

theta_R = 0:5:90 ;       % theta
phi_R = 0:5:360 ;      % phi
a=1;                  % z 축기준으로 떨어진 각 대입


cnt_det = zeros(6,1);
cnt_IM = zeros(6,360);
cnt_OM = zeros(6,360);
cnt_IMOM = zeros(6,360);
cnt = zeros(6,360);

Result=ones(360,1);
Coordinates_string = strings(1,19*73);
weightX = ones(6,1);
weightY = ones(6,1);
weightZ = ones(6,1);


tht_obj = [pi*10/180,pi*20/180,pi*30/180,pi*40/180,pi*50/180,pi*60/180];
%V= sin(tht_v1)*cos(phi_v),sin(tht_v1)*sin(phi_v),cos(tht_v1)
tht_v1 = pi/6;
tht_v2 = pi/3;
phi_v = 0;

vx=[1,0,0]; vy=[0,1,0]; vz=[0,0,1];
v1=[sin(tht_obj(1))*cos(phi_v),sin(tht_obj(1))*sin(phi_v),cos(tht_obj(1))]; %% theta = 10 deg
v2=[sin(tht_obj(2))*cos(phi_v),sin(tht_obj(2))*sin(phi_v),cos(tht_obj(2))]; %% theta = 20 deg
v3=[sin(tht_obj(3))*cos(phi_v),sin(tht_obj(3))*sin(phi_v),cos(tht_obj(3))]; %% theta = 30 deg
v4=[sin(tht_obj(4))*cos(phi_v),sin(tht_obj(4))*sin(phi_v),cos(tht_obj(4))]; %% theta = 40 deg
v5=[sin(tht_obj(5))*cos(phi_v),sin(tht_obj(5))*sin(phi_v),cos(tht_obj(5))]; %% theta = 50 deg
v6=[sin(tht_obj(6))*cos(phi_v),sin(tht_obj(6))*sin(phi_v),cos(tht_obj(6))]; %% theta = 60 deg
V= [v1; v2; v3; v4; v5; v6];
rootfolder=pwd;

%%
r1 = 1;
r2 = 3;
r3 = 9;
%% raw Detector Count
for k=2:2:12
    v = V(k/2,:);
    weightX = (v(1)*vx(1)+v(2)*vx(2)+v(3)*vx(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vx(1)^2+vx(2)^2+vx(3)^2);
    weightY = (v(1)*vy(1)+v(2)*vy(2)+v(3)*vy(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vy(1)^2+vy(2)^2+vy(3)^2);
    weightZ = (v(1)*vz(1)+v(2)*vz(2)+v(3)*vz(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vz(1)^2+vz(2)^2+vz(3)^2);
       
for p=0:0
    foldername_det = sprintf('center_Slice_view_slat13_onlydet(%d,%d)',theta_R(k+1),phi_R(p+1));
    branchfolder_det = strcat(rootfolder,'\',foldername_det); 
    cd(branchfolder_det);               
    filename = sprintf('center_Slice_view(%d,%d)_1 deg.png',theta_R(k+1),phi_R(p+1));  % - 대상 파일 이름으로 바꾸기
    src = imread(filename);  
    redChannel = src(:,:, 1);
    greenChannel = src(:,:, 2);
    blueChannel = src(:,:, 3);
    row = size(redChannel,1);
    column = size(redChannel,2);
       
      for l=1:row
          for j=1:column
            if redChannel(l,j)==204   % (0,1,0)
                if greenChannel(l,j)==225
                    cnt_det(a,1) = cnt_det(a)+r1*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==165
                    cnt_det(a,1) = cnt_det(a)+r2*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==105
                    cnt_det(a,1) = cnt_det(a)+r3*abs(1/weightY(1));                   
                end
                
            elseif redChannel(l,j)==230   % (1,0,0)
                if greenChannel(l,j)==225
                    cnt_det(a,1) = cnt_det(a)+r1*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==165
                    cnt_det(a,1) = cnt_det(a)+r2*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==105
                    cnt_det(a,1) = cnt_det(a)+r3*abs(1/weightX(1));                   
                end

            elseif redChannel(l,j)==255   % (0,0,1)
                 
                if greenChannel(l,j)==225
                    cnt_det(a,1) = cnt_det(a)+r1*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==165
                    cnt_det(a,1) = cnt_det(a)+r2*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==105
                    cnt_det(a,1) = cnt_det(a)+r3*abs(1/weightZ(1));                 
                end
            end
          end 
      end  % 1도씩 돌아간 그림 한개에서의 Total count            
end
a=a+1;
end
a=1;

%% IM 만의 영향
for k=2:2:12
    v = V(k/2,:);
    weightX = (v(1)*vx(1)+v(2)*vx(2)+v(3)*vx(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vx(1)^2+vx(2)^2+vx(3)^2);
    weightY = (v(1)*vy(1)+v(2)*vy(2)+v(3)*vy(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vy(1)^2+vy(2)^2+vy(3)^2);
    weightZ = (v(1)*vz(1)+v(2)*vz(2)+v(3)*vz(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vz(1)^2+vz(2)^2+vz(3)^2);

    for p=0:0
        foldername_IM = sprintf('center_Slice_view_slat13_ConIM(%d,%d)',theta_R(k+1),phi_R(p+1));
        branchfolder_IM = strcat(rootfolder,'\',foldername_IM); 
        cd(branchfolder_IM);    
        for i = 1:360 
        filename = sprintf('center_Slice_view(%d,%d)_%d deg.png',theta_R(k+1),phi_R(p+1),i);  % - 대상 파일 이름으로 바꾸기
        src = imread(filename);  
        redChannel = src(:,:, 1);
        greenChannel = src(:,:, 2);
        blueChannel = src(:,:, 3);
        row = size(redChannel,1);
        column = size(redChannel,2);
        
      for l=1:row
          for j=1:column
            if redChannel(l,j)==204   % (0,1,0)
                if greenChannel(l,j)==225
                    cnt_IM(a,i) = cnt_IM(a,i)+r1*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==165
                    cnt_IM(a,i) = cnt_IM(a,i)+r2*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IM(a,i) = cnt_IM(a,i)+r3*abs(1/weightY(1));                   
                end
                
            elseif redChannel(l,j)==230   % (1,0,0)
                if greenChannel(l,j)==225
                    cnt_IM(a,i) = cnt_IM(a,i)+r1*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==165
                    cnt_IM(a,i) = cnt_IM(a,i)+r2*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IM(a,i) = cnt_IM(a,i)+r3*abs(1/weightX(1));                   
                end

            elseif redChannel(l,j)==255   % (0,0,1)
                 
                if greenChannel(l,j)==225
                    cnt_IM(a,i) = cnt_IM(a,i)+r1*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==165
                    cnt_IM(a,i) = cnt_IM(a,i)+r2*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IM(a,i) = cnt_IM(a,i)+r3*abs(1/weightZ(1));                 
                end
            end
          end 
      end  % 1도씩 돌아간 그림 한개에서의 Total count        
    end
    a=a+1; 
    end
end
a=1;

%% OM 만의 영향
for k=2:2:12
    v = V(k/2,:);
    weightX = (v(1)*vx(1)+v(2)*vx(2)+v(3)*vx(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vx(1)^2+vx(2)^2+vx(3)^2);
    weightY = (v(1)*vy(1)+v(2)*vy(2)+v(3)*vy(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vy(1)^2+vy(2)^2+vy(3)^2);
    weightZ = (v(1)*vz(1)+v(2)*vz(2)+v(3)*vz(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vz(1)^2+vz(2)^2+vz(3)^2);
for p=0:0
    foldername_OM = sprintf('center_Slice_view_slat13_ConOO(%d,%d)',theta_R(k+1),phi_R(p+1));
    branchfolder_OM = strcat(rootfolder,'\',foldername_OM); 
    cd(branchfolder_OM);
        for i = 1:360         
        filename = sprintf('center_Slice_view(%d,%d)_%d deg.png',theta_R(k+1),phi_R(p+1),i);  % - 대상 파일 이름으로 바꾸기
        src = imread(filename);  
        redChannel = src(:,:, 1);
        greenChannel = src(:,:, 2);
        blueChannel = src(:,:, 3);
        row = size(redChannel,1);
        column = size(redChannel,2);
       
      for l=1:row
          for j=1:column
            if redChannel(l,j)==204   % (0,1,0)
                if greenChannel(l,j)==225
                    cnt_OM(a,i) = cnt_OM(a,i)+r1*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==165
                    cnt_OM(a,i) = cnt_OM(a,i)+r2*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==105
                    cnt_OM(a,i) = cnt_OM(a,i)+r3*abs(1/weightY(1));                   
                end
                
            elseif redChannel(l,j)==230   % (1,0,0)
                if greenChannel(l,j)==225
                    cnt_OM(a,i) = cnt_OM(a,i)+r1*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==165
                    cnt_OM(a,i) = cnt_OM(a,i)+r2*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==105
                    cnt_OM(a,i) = cnt_OM(a,i)+r3*abs(1/weightX(1));                   
                end

            elseif redChannel(l,j)==255   % (0,0,1)
                 
                if greenChannel(l,j)==225
                    cnt_OM(a,i) = cnt_OM(a,i)+r1*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==165
                    cnt_OM(a,i) = cnt_OM(a,i)+r2*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==105
                    cnt_OM(a,i) = cnt_OM(a,i)+r3*abs(1/weightZ(1));                 
                end
            end
          end 
      end  % 1도씩 돌아간 그림 한개에서의 Total count        
    end

end
a=a+1;
end
a=1;


%% IM 또는 OM의 영향 
for k = 2:2:12
    v = V(k/2,:);
    weightX = (v(1)*vx(1)+v(2)*vx(2)+v(3)*vx(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vx(1)^2+vx(2)^2+vx(3)^2);
    weightY = (v(1)*vy(1)+v(2)*vy(2)+v(3)*vy(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vy(1)^2+vy(2)^2+vy(3)^2);
    weightZ = (v(1)*vz(1)+v(2)*vz(2)+v(3)*vz(3))/sqrt(v(1)^2+v(2)^2+v(3)^2) * sqrt(vz(1)^2+vz(2)^2+vz(3)^2);
    
for p = 0:0
    Coordinates_string(a) = sprintf("(%d , %d)",theta_R(k+1),phi_R(p+1));    
    cd(rootfolder);

    foldernameIMOM = sprintf('center_Slice_view_slat13_DetC(%d,%d)',theta_R(k+1),phi_R(p+1));
    branchfolderIMOM = strcat(rootfolder,'\',foldernameIMOM);
    cd(branchfolderIMOM);
        
    for i = 1:360    
        filename = sprintf('center_Slice_view(%d,%d)_%d deg.png',theta_R(k+1),phi_R(p+1),i);  % - 대상 파일 이름으로 바꾸기
        src = imread(filename);  
        redChannel = src(:,:, 1);
        greenChannel = src(:,:, 2);
        blueChannel = src(:,:, 3);
        row = size(redChannel,1);
        column = size(redChannel,2);
       
      for l=1:row 
          for j=1:column
            if redChannel(l,j)==204   % (0,1,0
                if greenChannel(l,j)==225
                   cnt_IMOM(a,i) = cnt_IMOM(a,i)+r1*abs(1/weightY(1));
                aa
                elseif greenChannel(l,j)==165
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r2*abs(1/weightY(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r3*abs(1/weightY(1));                   
                end
                
            elseif redChannel(l,j)==230   % (1,0,0)
                if greenChannel(l,j)==225
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r1*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==165
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r2*abs(1/weightX(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r3*abs(1/weightX(1));                   
                end

            elseif redChannel(l,j)==255   % (0,0,1)
                 
                if greenChannel(l,j)==225
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r1*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==165
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r2*abs(1/weightZ(1));
                
                elseif greenChannel(l,j)==105
                    cnt_IMOM(a,i) = cnt_IMOM(a,i)+r3*abs(1/weightZ(1));                 
                end
            end
          end 
      end  % 1도씩 돌아간 그림 한개에서의 Total count 
    cnt(a,i) = cnt_IMOM(a,i)+ 0.1920*(cnt_det(a)-cnt_OM(a,i)-cnt_IM(a,i)+cnt_IMOM(a,i)) + 0.4382*(cnt_OM(a,i)-cnt_IMOM(a,i)) + 0.4382*(cnt_IM(a,i)-cnt_IMOM(a,i));
    Result(i)=cnt(a,i);     
    end
end    
    a=a+1;
    Resultfilename = sprintf('Result_rev_slat13_%d_r300_TC_SC139.xlsx',theta_R(k+1));
    xlswrite(Resultfilename,Result,'ModResult');    
end
cd(rootfolder);



