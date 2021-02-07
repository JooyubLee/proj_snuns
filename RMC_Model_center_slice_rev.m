n = 11; % slat °³¼ö
p = 30;
%x->z y->x z->y
%
theta = linspace(0,pi,1+p).';
phi_base = linspace(0,9.5*pi/180,2);
phi_base2 = linspace(170.5*pi/180,pi,2);
phi_Main = linspace(16.5*pi/180,163.5*pi/180,2*n);
phi=[phi_base,phi_Main,phi_base2];
%
    
    x = sin(theta)*cos(phi)*100;
    y = sin(theta)*sin(phi)*100;
    z = cos(theta)*ones(1,2*(n+2))*100;
        
V = [x(:),z(:),y(:)];
F = bsxfun(@plus,(1:p).',0:2*(p+1):size(V,1)-1);
F = F(:);
F(:,2:4) = [1+F,F+2+p,F+1+p];

%
aa=patch('Faces',F, 'Vertices',V, 'FaceColor','green', 'EdgeColor','none');
%aa.FaceVertexAlphaData = 0;    % Set constant transparency 
%aa.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
hold on;

    x = sin(theta)*cos(phi)*95;
    y = sin(theta)*sin(phi)*95;
    z = cos(theta)*ones(1,2*(n+2))*95;

V = [x(:),z(:),y(:)];
F = bsxfun(@plus,(1:p).',0:2*(p+1):size(V,1)-1);
F = F(:);
F(:,2:4) = [1+F,F+2+p,F+1+p];

%
bb=patch('Faces',F, 'Vertices',V, 'FaceColor','green', 'EdgeColor','none');
%bb.FaceVertexAlphaData = 0;    % Set constant transparency 
%bb.FaceAlpha = 'flat' ;
hold on;

    x = sin(theta)*cos(phi)*30;
    y = sin(theta)*sin(phi)*30;
    z = cos(theta)*ones(1,2*(n+2))*30;

V = [x(:),z(:),y(:)];
F = bsxfun(@plus,(1:p).',0:2*(p+1):size(V,1)-1);
F = F(:);
F(:,2:4) = [1+F,F+2+p,F+1+p];

%
cc=patch('Faces',F, 'Vertices',V, 'FaceColor','blue', 'EdgeColor','none');
%cc.FaceVertexAlphaData = 0;    % Set constant transparency 
%cc.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
hold on;


    x = sin(theta)*cos(phi)*25;
    y = sin(theta)*sin(phi)*25;
    z = cos(theta)*ones(1,2*(n+2))*25;

V = [x(:),z(:),y(:)];
F = bsxfun(@plus,(1:p).',0:2*(p+1):size(V,1)-1);
F = F(:);
F(:,2:4) = [1+F,F+2+p,F+1+p];

%
dd=patch('Faces',F, 'Vertices',V, 'FaceColor','blue', 'EdgeColor','none');
%dd.FaceVertexAlphaData = 0;    % Set constant transparency 
%dd.FaceAlpha = 'flat' ;          % Interpolate to find face transparency
hold on;
%%%%% detector 1 x 1 x 1 cm^3
%{
[X Y] = meshgrid(linspace(-5,5,2),linspace(-5,5,2)); % 
Z = 5*ones(length(X),length(Y));
surf(X,Y,Z+5,'edgecolor','black','FaceColor',[1 0 0]); %¶Ñ²± V(0,0,1) 0
hold on;
surf(Y,Z,X+5,'edgecolor','black','FaceColor',[1 0 0.1]); % ¶Ñ²± V(0,1,0) 26
surf(Z,X,Y+5,'edgecolor','black','FaceColor',[1 0 0.2]); % ¶Ñ²± V(1,0,0) 51
surf(-Y,-Z,-X+5,'edgecolor','black','FaceColor',[1 0 0.1]); % ¶Ñ²± V(0,1,0) 26
surf(-Z,-X,-Y+5,'edgecolor','black','FaceColor',[1 0 0.2]); % ¶Ñ²± V(1,0,0) 51
%}

%{
[X,Y] = meshgrid(linspace(-5,5,200),linspace(-5,5,200));
Z = 5*ones(length(X),length(Y));
surf(X,Y,Z+5/3,'edgecolor','red');
hold on;
surf(Y,Z,X+5/3,'edgecolor','red');
surf(Z,X,Y+5/3,'edgecolor','red');
surf(-Y,-Z,-X+5/3,'edgecolor','red');
surf(-Z,-X,-Y+5/3,'edgecolor','red');
surf(-Y,-Z,-X+5/3,'edgecolor','red');

surf(X,Y,Z+5/3,'edgecolor','yellow');
surf(Y,Z,X+5/3,'edgecolor','yellow');
surf(Z,X,Y+5/3,'edgecolor','yellow');
surf(-Y,-Z,-X+5/3,'edgecolor','yellow');
surf(-Z,-X,-Y+5/3,'edgecolor','yellow');
surf(-Y,-Z,-X+5/3,'edgecolor','yellow');

surf(X,Y,Z+5/3,'edgecolor','green');
surf(Y,Z,X+5/3,'edgecolor','green');
surf(Z,X,Y+5/3,'edgecolor','green');
surf(-Y,-Z,-X+5/3,'edgecolor','green');
surf(-Z,-X,-Y+5/3,'edgecolor','green');
surf(-Y,-Z,-X+5/3,'edgecolor','green');
%}
%Untitled
%% Detector Modeling

N=10; % Æò¸éÀ» ³ª´­ È½¼ö, maximum 255

[X Y] = meshgrid(linspace(-5,5,N+1),linspace(-5,5,N+1)); % 
Z = 5*ones(length(X),length(Y));
Source = [1/2 0 sqrt(3)/2] * 1100;  %% (pi/6,0)
Angle=60; % phi ·Î È¸Àü
Q = fix(N/2);

% X=5
R1X = 5*ones(N,N);
%R1Y = ones(1,N^2);
%R1Z = ones(1,N^2);

% X=-5
R2X = -5*ones(N,N);
%R2Y = ones(1,N^2);
%R2Z = ones(1,N^2);

% Y=5
%R3X = ones(1,N^2);
R3Y = 5*ones(N,N);
%R3Z = ones(1,N^2);

% Y=-5
%R4X = ones(1,N^2);
R4Y = -5*ones(N,N);
%R4Z = ones(1,N^2);

% Z=10
%R5X = ones(1,N^2);
%R5Y = ones(1,N^2);
R5Z = 10*ones(N,N);

% ¼±¿øÀÇ Åõ°ú °Å¸® ¼³Á¤  
M_D_X5 = ones(N,N); M_D_Z = ones(N,N);
r=1;

for k=1:N^2        
    % 30µµ
    R1y = linspace(-5+5/(N),5-5/(N),N);
    R1z = linspace(0+5/(N),10-5/(N),N);
    [R1Y,R1Z] = meshgrid(R1y,R1z);
    
    % 60µµ
    R2y = linspace(-5+5/(N),5-5/(N),N);
    R2z = linspace(0+5/(N),10-5/(N),N);
    [R2Y,R2Z] = meshgrid(R2y,R2z);
    
    R3x = linspace(-5+5/(N),5-5/(N),N);
    R3z = linspace(0+5/(N),10-5/(N),N);
    [R3X,R3Z] = meshgrid(R3x,R3z);
    
    R4x = linspace(-5+5/(N),5-5/(N),N);
    R4z = linspace(0+5/(N),10-5/(N),N);
    [R4x,R4z] = meshgrid(R4x,R4z);
    
    R5x = linspace(-5+5/(N),5-5/(N),N);
    R5y = linspace(-5+5/(N),5-5/(N),N);
    [R5X,R5Y] = meshgrid(R5x,R5y);
        
    
    rayVec = (Source - [R1X(k) R1Y(k) R1Z(k)]) ./ sqrt((Source(1)-R1X(k))^2+(Source(2)-R1Y(k))^2+(Source(3)-R1Z(k))^2);
    rayVec2 = (Source - [R5X(k) R5Y(k) R5Z(k)]) ./ sqrt((Source(1)-R5X(k))^2+(Source(2)-R5Y(k))^2+(Source(3)-R5Z(k))^2);

    
    %v0 = R1Z(k) / rayVec(3);  
    
    %tx5 = (R1X(k)-5) / rayVec(1);
    %ux5 = (R1Y(k)-5) / rayVec(2);
    %vx5 = (R1Z(k)-5) / rayVec(3);
    
    R1tnx5 = (R1X(k)+5) / rayVec(1);
    R5tnx5 = (R5X(k)+5) / rayVec2(1);
    %unx5 = (R1Y(k)+0) / rayVec(2);
    %vnx5 = (R1Z(k)+0) / rayVec(3);

    %tny5 = (R1X(k)+0) / rayVec(1);
    R1uny5 = (R1Y(k)+5) / rayVec(2);
    R5uny5 = (R5Y(k)+5) / rayVec2(2);
    %vny5 = (R1Z(k)+0) / rayVec(3);
    
    %ty5 = (R1X(k)) / rayVec(1);
    R1uy5 = (R1Y(k)-5) / rayVec(2);
    R5uy5 = (R5Y(k)-5) / rayVec2(2);
    %vy5 = (R1Z(k)) / rayVec(3);

    %tz0 = (R1X(k)) / rayVec(1);
    %uz0 = (R1Y(k)) / rayVec(2);
    R1vz0 = (R1Z(k)) / rayVec(3);
    R5vz0 = (R5Z(k)) / rayVec2(3);

    
    
    K = [ abs(R1tnx5) abs(R1uny5) abs(R1uy5) abs(R1vz0) ];
    [m,T] = min(K);
    
    L = [ abs(R5tnx5) abs(R5uny5) abs(R5uy5) abs(R5vz0) ];
    [m2,T2] = min(L); 
    
    %if(m==4)

    
    if (T==1)
    D(1) = -5;
    D(2) = R1Y(k)-m*rayVec(2);
    D(3) = R1Z(k)-m*rayVec(3);
    D(4) = -5;
    D(5) = R5Y(k)-m2*rayVec2(2);
    D(6) = R5Z(k)-m2*rayVec2(3);
       
    elseif (T==2)    
    D(1) = R1X(k)-m*rayVec(1);
    D(2) = -5;
    D(3) = R1Z(k)-m*rayVec(3);        
    D(4) = R5X(k)-m2*rayVec2(1);
    D(5) = -5;
    D(6) = R5Z(k)-m2*rayVec2(3);      
    
    elseif (T==3)    
    D(1) = R1X(k)-m*rayVec(1);
    D(2) = 5;
    D(3) = R1Z(k)-m*rayVec(3);        
    D(4) = R5X(k)-m2*rayVec2(1);
    D(5) = 5;
    D(6) = R5Z(k)-m2*rayVec2(3); 
    
    else 
    D(1) = R1X(k)-m*rayVec(1);
    D(2) = R1Y(k)-m*rayVec(2);
    D(3) = 0;    
    D(4) = R5X(k)-m2*rayVec2(1);
    D(5) = R5Y(k)-m2*rayVec2(2);
    D(6) = 0;    
    
    end
    
    D_X5 = sqrt((R1X(k)-D(1))^2+(R1Y(k)-D(2))^2+(R1Z(k)-D(3))^2);
    D_Z = sqrt((R5X(k)-D(4))^2+(R5Y(k)-D(5))^2+(R5Z(k)-D(6))^2);
    
    fprintf('Pixel Center (%f %f %f) -> Minimum Value is %f Region %f\n',R1X(k),R1Y(k),R1Z(k),m,T);     
    fprintf('Distance = %f \n', D_X5); 
    fprintf('Pixel Center (%f %f %f) -> Minimum Value is %f %f\n',R5X(k),R5Y(k),R5Z(k),m2,T2);             
    fprintf('Distance = %f \n\n', D_Z);  
    
   
        %M_D_Z (N+1-r,rem(k,N)+1) = D_Z  ;
    
    if (rem(k,N)==0)
       M_D_X5(1,r) = D_X5 ;
       M_D_Z(1,r) = D_Z ;       
        r=r+1;
    else 
        M_D_X5(N+1-rem(k,N),r) = D_X5 ;
        M_D_Z(N+1-rem(k,N),r) = D_Z ;
    end 
       
end

for i = 1:N
    for k = 1:N
        W1 = M_D_X5(N+1-i,k);
        W2 = M_D_Z(i,k);
        if (0<=W1 && W1<5)              % ¹ÝÀÀµµ °ª
            W1=1;                      
        elseif (5<=W1 && W1<9)
            W1=3;                       % ¹ÝÀÀµµ °ª
        elseif (9<=W1)
            W1=5;
        end

        surf(Z(i:i+1,i:i+1),X(k:k+1,k:k+1),Y(i:i+1,i:i+1)+5,'edgecolor','black','FaceColor',[0.9 1-30*W1/255 1]); % ¶Ñ²± V(1,0,0) 51        
        hold on;    
        if (0<=W2 && W2<5)
            W2=1;                        % ¹ÝÀÀµµ °ª
        elseif (5<=W2 && W2<9)
            W2=3;
        elseif (9<=W2)
            W2=5;
        end        
        surf(X(i:i+1,i:i+1),Y(k:k+1,k:k+1),Z(i:i+1,i:i+1)+5,'edgecolor','black','FaceColor',[1 1-30*W2/255 1]); %¶Ñ²± V(0,0,1) 0        
    end
end
%%
hold on;

axis vis3d off

%%%%%%%%%%%% Mask Rotation ¼öÇà %%%%%%%%    
%%%% Viewpoint ¼³Á¤%%
theta_R = [10,20,30,40,50,60];       % theta,
phi_R = 0 ;      % phi
rootfolder = pwd;   % ÇöÀç ÀÛ¾÷Æú´õ À§Ä¡ ¼³Á¤
 for k=0
    for j=0:0
           cd(rootfolder)
           foldername = sprintf('center_Slice_view_slat13_DetC(%d,%d)',theta_R(k+1),phi_R(j+1));
           mkdir (foldername);
           branchfolder = strcat(rootfolder,'\',foldername);
           cd(branchfolder);
        for i=1:360 
           %f = figure ;
           rotate(aa,[0,0,1],1);
           rotate(bb,[0,0,1],1);
           rotate(cc,[0,0,1],1);
           rotate(dd,[0,0,1],1);
           view(90+phi_R(j+1),90-theta_R(k+1)); % - ¹Ù²Ù±â
           filename = sprintf('center_Slice_view(%d,%d)_%d deg',theta_R(k+1),phi_R(j+1),i);    
           print(filename,'-dpng','-r300');
        end
     end
  end    
        
 % slat 13 => slat 13°³ :  9.5/7*11/9.5 °£°ÝÀ¸·Î ¼³Á¤
        
        
