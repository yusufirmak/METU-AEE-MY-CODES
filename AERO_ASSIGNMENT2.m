%PLOTTING Airfoil AND CF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY IMPORTANT!! PLEASE CHANGE THE BELOW LINES ACCORDING TO THE LOCATION
% OF MY AIRFOIL AND CP DATA FILES IN YOUR COMPUTER HOCAM. TO DO THAT
% EASILY: WRITE "airfoil_2522597" in the search bar next to the windows logo, right click and click
% to the "copy path".Then paste it into the below line↓ Do the same for cpx 
fidAirfoil= fopen('C:\Users\yusuf\Downloads\airfoil_2522597.txt');

dataBuffer= textscan(fidAirfoil, '%f %f', 'CollectOutput',1,...
                                 'Delimiter','','HeaderLines',0);



fclose(fidAirfoil);
%delete(saveFlnmAF);

XB=dataBuffer{1}(:,1);
YB=dataBuffer{1}(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% VERY IMPORTAN!! CHANGE THE PATH for "cpx_2522597.txt" ↓ %%%%
fidCP= fopen('C:\Users\yusuf\Downloads\cpx_2522597.txt');

dataBuffer= textscan(fidCP, '%f %f %f','Headerlines',3,...
                              'CollectOutput',1,...
                              'Delimiter','');


fclose(fidCP);
%delete(saveFlnmCp);

X_0= dataBuffer{1,1}(:,1);
Y_0= dataBuffer{1,1}(:,2);
Cp_0= dataBuffer{1,1}(:,3);


XB_U= XB(YB >= 0);
XB_U= flip(XB_U);
XB_L= XB(YB < 0);
YB_U= YB(YB >= 0);
YB_U= flip(YB_U);
YB_L= YB(YB < 0);

Cp_U = Cp_0(YB >= 0);
Cp_U= flip(Cp_U);
Cp_L = Cp_0(YB < 0);
X_U= X_0(YB >= 0);
X_U= flip(X_U);
X_L= X_0(YB < 0);

figure(3);
cla; hold on; grid off;
set(gcf,'Color','White');
set(gca,'FontSize',12);
plot(XB_U,YB_U,'b.-');
plot(XB_L,YB_L, 'r.-');
xlabel('X Coordinate');
ylabel('Y Coordinate');
axis equal;

figure(11);
cla; hold on; grid off;
set(gcf,'Color','White');
set(gca,'FontSize',12);
plot(XB_U,Cp_U,'b.-');
plot(XB_L,Cp_L, 'r.-');
xlabel('X Coordinate');
ylabel('Cp Distribution');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
N = length(XB_L);
dY_L= diff(YB_L);
dX_L= diff(XB_L);
der_L= dY_L./dX_L;
der_L(N) = der_L(N-1);

M= length(XB_U);
dY_U= diff(YB_U);
dX_U= diff(XB_U);
der_U= dY_U./dX_U; 
der_U(M) = der_U(M-1);


cn= trapz(X_L,Cp_L)-trapz(X_U,Cp_U)

CP_U_der_U= Cp_U.*der_U;
CP_L_der_L= Cp_L.*der_L;
ca= trapz(X_U,CP_U_der_U)-trapz(X_L,CP_L_der_L) 

alfa=3
cl=cn*cosd(alfa)-ca*sind(alfa)

cd=cn*sind(alfa)+ca*cosd(alfa)

cm_LE= trapz(X_U,Cp_U.*X_U)-trapz(X_L,Cp_L.*X_L) + trapz(X_U,CP_U_der_U.*YB_U)-trapz(X_L,CP_L_der_L.*YB_L) 
cm_quarter= cm_LE + (cl/4)
 
