function results = adjust_main(yiparaN,yipara,unObsN,jdObsN,jdObs,numJ,numX,xdObs,xdObsN,app_RobusteParameter_Value,app_K1Parameter_Value,app_K2Parameter_Value,app_Adjustment_MethodButtonGroup_SelectedObject_Text,app_InstrumentParameter_Value,app_Datum_PointCheckBox_Value,app_Basic_PointCheckBox_Value,app_PointofreferenceCheckBox_Value,app_ShortBaselineCheckBox_Value,app_residualCheckBox_Value,app_instrumentCheckBox_Value)             
feature('DefaultCharacterSet','GBK'); 
slCharacterEncoding('GBK');
f = waitbar(0.1,'Computation is processing');
numyp = length(yiparaN);                 

stN = zeros(numyp,2);
for i = 1:numyp
    for j = 1:9
        if j < 3
           if yipara(i,j) == 1
              stN(i,1) = stN(i,1) +1;
           end
        else
           if yipara(i,j)  == 1
              stN(i,2) = stN(i,2) +2;
           end
        end
    end
end
waitbar(0.2,f)
%************************************************
gravN = length(unObsN);                 
unknowN = gravN + sum(sum(stN));        
%************************************************
ss = 0;
for i = 1:numJ
    for j = 1:length(unObsN)
        if strcmp(jdObsN(i,1),unObsN(j)) == 1
           ss = ss + 1;
           jdObsn(ss,1) = jdObsN(i,1);
           jdObss(ss,1:2) = jdObs(i,:);
         end
    end
end
numJ = ss;
%************************************************
%******************************                           
%************************************************
waitbar(0.3,f);
Ajd = zeros(numJ,unknowN);             
Pjd = diag(jdObss(:,2));               
for i =1:numJ
    pp = find(strcmp(jdObsn(i,1),unObsN) == 1);
    Ajd(i,pp) = 1;
end
init = 9.78*10^8;conS = 10^3;
Ljd = (jdObss(:,1)-init)/conS;         
%************************************************
%******************************                           
%************************************************
Axd = zeros(numX,unknowN);              
Lxd = zeros(numX,1);
Pxd =  diag(xdObs(:,8));                
T = [220 110 73.33 36.67 7.33 3.67 1];                                
              
%***********************************************************************
for i = 1:numX
    %*********************************************
    %*********************************************

    poinS = find(strcmp(xdObsN(i,2),unObsN) == 1);
    poinE = find(strcmp(xdObsN(i,3),unObsN) == 1);
    Axd(i,poinS) = -1;
    Axd(i,poinE) = 1;

    ge = find(strcmp(xdObsN(i,1),yiparaN(:,1)) == 1);
    for j = 1:2

        if yipara(ge,j) == 1
           if j == 1
              if ge >1
                 Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+j) = (xdObs(i,3) - xdObs(i,4))/conS;
              else
                 Axd(i,gravN+j) =  (xdObs(i,3) - xdObs(i,4))/conS; 
              end
           else
              if ge >1
                 Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+j) = (xdObs(i,3)^2 - xdObs(i,4)^2)/conS^2;  
              else
                 Axd(i,gravN+j)  =   (xdObs(i,3)^2- xdObs(i,4)^2)/conS^2; 
              end
           end
        end
   end
   st = 0;
   %********************************************

   for j = 1:7
       if yipara(ge,2+j) == 1
          st = st + 1;
          if st > 0
             if ge >1
                Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+stN(ge,1)+2*st-1) = (cos(xdObs(i,1)*2*pi/T(j))-cos(xdObs(i,2)*2*pi/T(j)));
                Axd(i,gravN+sum(sum(stN(1:ge-1,:)))+stN(ge,1)+2*st) = (sin(xdObs(i,1)*2*pi/T(j))-sin(xdObs(i,2)*2*pi/T(j)));
             else
                Axd(i,gravN+stN(ge,1)+2*st-1) = (cos(xdObs(i,1)*2*pi/T(j))-cos(xdObs(i,2)*2*pi/T(j)));
                Axd(i,gravN+stN(ge,1)+2*st) =(sin(xdObs(i,1)*2*pi/T(j))-sin(xdObs(i,2)*2*pi/T(j)));
             end
          end
       end
   end
end
waitbar(0.4,f);
A = [Ajd;Axd];
L = [Ljd' Lxd']';
P = diag([diag(Pjd)' diag(Pxd)'])';
save('data.mat','A','L','P');
if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'LS')
   NN = A'*P*A;NL = A'*P*L;
   gravX = NN\NL; 
   Vjd = -(Ljd - Ajd*gravX)*conS;                   
   Vxd = -(Lxd - Axd*gravX)*conS;                   
   V = [Vjd' Vxd']';
   %****************************************
   gravX(1:gravN) = gravX(1:gravN)*conS+ init;      
   rms = Vjd'*Pjd*Vjd + Vxd'*Pxd*Vxd;               
   freedom = numJ + numX - unknowN;                 
   sigma = sqrt(rms/freedom);                       
   QQ  = inv(NN);                                   
   varX = sigma*sqrt(QQ);                           
   varM = sigma*sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN); 
else strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'Robust')
    e = app_RobusteParameter_Value;
    k1 =app_K1Parameter_Value ;k2= app_K2Parameter_Value;
    [X_robust,Pend,iter] = robust(A,L,P,k1,k2,numJ,e);
    gravX = X_robust;
    Vxd = -(Lxd - Axd*gravX)*conS;    
    Vjd = -(Ljd - Ajd*gravX)*conS;    
    V = [Vjd' Vxd']';                 
    %**************************************************
    %**************************************************
    gravX(1:gravN)= X_robust(1:gravN)*conS+ init;
    num0 = 0;
    for i = numJ+1:numJ+numX 
        if Pend(i,i) == 10^-10
           num0 = num0 + 1;
        end
    end
    Pend_point = find(diag(Pend) == 10^-10) - numJ;
    sigma = sqrt((V'*Pend*V)/(numJ + numX - unknowN));  
    QQ  = inv(A'*Pend*A);                                      
    varX = sigma*sqrt(QQ);                                     
    varM = sigma*(sqrt(sum(diag(QQ(1:gravN,1:gravN)))/gravN)); 
    app_iter_number_Value = num2str(iter);
end
waitbar(0.5,f);
%**********************************************************************

if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'LS')
   PP = P;
else
   PP = Pend;
end
PV = diag(1./diag(PP)) - A*((A'*PP*A)\A');
PVV = sigma*sqrt(diag(PV)); 

Cpath = strsplit(app_InstrumentParameter_Value,'\'); 
if app_residualCheckBox_Value
   [stat_chao,numVV1,numVV2,stat_V,stat_VJ,stat_VX] = stat_res(V,PVV,Vjd,Vxd,numJ);
end
if app_instrumentCheckBox_Value
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'Instrument.txt'),'w');
   AS = xdObsN;
   num_AS = length(AS(:,1));
   for i = 1:num_AS
       AS(i,4) = num2cell(V(numJ+i));
   end
   AS = sortrows(AS,1);
   num_AS = length(AS(:,1));
   aa ={};
   aa{1} = AS{1,1}(1);
   %************************************************
   for i = 2:num_AS
       aa{i} = AS{i,1}(1);
   end
   leii = unique(aa);  
   %***********************************************
   for j = 1:length(leii)
       ASJ = [];
       posi_yiqi = find(strcmp(aa,leii{j}) == 1);
       ASJ = cell2mat(AS(posi_yiqi,4));
       num_lei = length(posi_yiqi);
       yiqi_pre = sqrt(sum((ASJ - mean(ASJ)).^2)/num_lei/(num_lei-1));
       if j == 1
          fprintf(fid,' instrument     number of measurement     precision\n');
       end
       fprintf(fid, '%5s     %5d    %7.3f\n',strcat(leii{j},'instrument'),num_lei, yiqi_pre);
       figure(j+6)
       histogram(ASJ(intersect(find(ASJ<150), find(ASJ>-150))));
       title(strcat(leii{j},'residual of instrument='));xlabel('residual(uGal)');ylabel('number');
   end
   msgbox('Operation Successfully');   
   fclose(fid);
end                 
waitbar(0.6,f);
%*******************************************************************
[yiE,QyiE,An,alf,QAn,Qalf]= stat_yiqi(yipara,gravX,varX,stN,conS,gravN,numyp);
%***************************************************

if strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'robust') ==1   
   fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'robust_outlier.txt'),'w');
   for i = 1:num0
        fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{Pend_point(i),1:3},xdObs(Pend_point(i),:));
   end
elseif strcmp(app_Adjustment_MethodButtonGroup_SelectedObject_Text,'outlier') ==1 
       fid = fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'quasi_detection.txt'),'w');
       for i = 1:num0
           fprintf(fid,' %5s  %6s %6s %11.3f %11.3f %11.3f %11.3f %11.3f %2d %2d %2d\n',xdObsN{posi_det(i),1:3},xdObs(posi_det(i),:));
       end
       fclose(fid)  
end
waitbar(0.7,f);
%*****************************************************
fid=fopen(strrep(app_InstrumentParameter_Value,Cpath{end},'Result of free network adjustment.txt'),'w','n','UTF-8');
fprintf(fid,'%7s   %7d      %7s   %4d\n','1.Number of measurement for adjustment computation:  ',numX+numJ,' 2.Number of gravity point and instrument:  ',unknowN);
fprintf(fid,'%7s  %7f    %7s   %7f\n\n','3.Sigma:  ',sigma,'4.Vaiance of gravity point:          ',varM);
fprintf(fid,'%6s %6s    %7s     %7s\n','5.The ending gravity value of observaiton point:  ','Seriers number of point   ','gravity value ','precision');
results = {};
for i = 1:gravN
    results(i,1:3) = {unObsN{i}  gravX(i) varX(i,i)};
end
results = sortrows(results,1);
for i =1:gravN
    if i == gravN
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n\n',results{i,:}); 
    else
       fprintf(fid,'                    %5s       %11.3f     %6.3f\n',results{i,:});     %输出各个重力点的精度
    end
end

C = firstalp(results);
if app_Datum_PointCheckBox_Value                         
   CC1 = find(strcmp(C,'A') == 1);
   num_Datumpoint = length(CC1);
   abs_pre = sqrt(sum(cell2mat(results(CC1,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.precision statistics of gravity point :');
   fprintf(fid,'The type of gravity point     number      precision\n');
   fprintf(fid,'                   datum point     %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_Basic_PointCheckBox_Value      
    CC2 = find(strcmp(C,'B') == 1);
    abs_pre = sqrt(sum(cell2mat(results(CC2,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.precision statistics of gravity point :');
    fprintf(fid,'The type of gravity point     number      precision\n');
    fprintf(fid,'                   basic point      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if app_PointofreferenceCheckBox_Value   
   CC3 = find(strcmp(C,'C') == 1);
   num_Datumpoint = length(CC3);
   abs_pre = sqrt(sum(cell2mat(results(CC3,3)).^2)/num_Datumpoint); 
   fprintf(fid,'6.precision statistics of gravity point :');
   fprintf(fid,'The type of gravity point     number      precision\n');
   fprintf(fid,'                   Leading point      %4d     %7.3f\n',num_Datumpoint,abs_pre);
end
if  app_ShortBaselineCheckBox_Value     
    CC4 = find(strcmp(C4,'S') == 1);
    num_Datumpoint = length(CC4);
    abs_pre = sqrt(sum(cell2mat(results(CC4,3)).^2)/num_Datumpoint); 
    fprintf(fid,'6.precision statistics of gravity point :');
    fprintf(fid,'The type of gravity point     number      precision\n');
    fprintf(fid,'                   Short baseline      %4d     %7.3f\n',num_Datumpoint,abs_pre);   
end
waitbar(0.8,f);
%***************************************************************

fprintf(fid,'\n');
fprintf(fid,'%12s  %4d  %4s %4s   %9s\n','7.grid value and parameter of periodic error parameter:',unknowN-gravN,'C1','C2','Xi,Yi(i=1:7)');
for i= 1:numyp
    if i == numyp
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n\n',yiparaN{i,1},yiE(i,:));
    else
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n',yiparaN{i,1},yiE(i,:));
    end
end
fprintf(fid,'%12s\n','8.The sigma of grid value and parameter of periodic error parameter:');
for i= 1:numyp
    if i == numyp
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n\n',yiparaN{i,1},QyiE(i,:));
    else
       fprintf(fid,'    %5s   %13.10f   %13.10f  %13.10f  %13.10f  %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f %13.10f   %13.10f  %13.10f  %13.10f\n',yiparaN{i,1},QyiE(i,:));
    end     
end
waitbar(0.9,f);
%*****************************************************

waitbar(1,f,{'The processing is ending' 'Operation Successfully'});
results = 0;