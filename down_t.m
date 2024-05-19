%-------------------------Downweight function-----------------------------
function tao = down_t(abVV,k1,k2,numJ)
%�����Ƕ�����Ȩ�����µĿ������
%P:weight matrix
%k1,k2:ΪIGGIII�����Ĳ���
%A:��ʾϵ������
%L:�����Ҷ���
%down_p��Ϊ���յ�Ȩ
%----------------------ȷ����Ȩ����-------------------------
n = length(abVV);
tao = ones(n,1);
for i = numJ+1:n
       if abVV(i) <= k1
          tao(i,1) = 1;
       elseif abVV(i) >k1 && abVV(i) <=k2
          tao(i,1) = k1/abVV(i)*((k2-abVV(i))/(k2-k1))^2;
       else
          tao(i,1) = 10^-10;
       end
end
%***********************************************
% [n,m] = size(A);
% N = A'*P*A; 
% S = A*(N\A');
% R = eye(n) - S*P;         %ƽ��������
% V = R*L;                  %�в�
% PV = (diag(1./diag(P))-S);%�в��Ȩ����
% num = 0;
%**********************************************
%�����׼���в��ͳ��Ȩ��Ϊ��ĸ���
% for i = 1:n
%       VV(i) = V(i)/sqrt(PV(i,i));   %�����׼���в����û�д�����λȨ����
% %     if P(i,i) <= 10^-8
% %         num = num +1;             %ͳ��ȨΪ��ĸ���
% %     end
% end
%**********************************************
%��׼���в����sigma�ķ�����������
% mad_sigma = sqrt((V'*P*V)/(n-m-num));   %�������λȨ����ȷ����λȨ�����޳�ȨΪ��Ĺ۲��������ɶ����ȥ��ȨΪ��ĸ���
% % mad_sigma = mad(abs(VV),1)*1.4826;     %������λ����ȷ����λȨ�����
% % mad_sigma =sqrt(median(VV.^2))*1.4826; %���ÿ���LMS(�������¶�����λ������);
% %**********************************************
% % VV = VV/mad_sigma;                       %�õ����յı�׼���в�

%**********************************************
% ���ݽ�ȨIGGIII���������㽵Ȩ����
% tao = zeros(n,1);
% for i = 1:n
%        abVV = abs(VV(i));
%        if abVV  <= k1
%           tao(i,1) = 1;
%        elseif abVV >k1 && abVV <=k2
%           tao(i,1) = k1/abVV*((k2-abVV)/(k2-k1))^2;
%        else
%           tao(i,1) = 10^-8;
%        end
% end
%***************************************
%���㷽����ȷ����397��ʼ��ʱ��ű�Ȩ�����ڷɻ��̶�Ȩʱ
%***************************************
%����1������֪������׼ֵ��Ȩȷ��Ϊ���޴󣬹��ڿ�չ�������ʱ����Ȩ����ı�
%������֪������׼�㹲8��
% tao(2011:end) = 1;                  %ȷ����Թ۲���ǰ396���۲��Ȩ����:����2ѡ��
%***************************************
%����2������֪������׼ֵ������֪����н����Ҳ�����ƽ����㣬���ڿ�չ�������ʱ�����迼�Ƕ�Ӧ��Ȩ��
% tao(numJ+1:numJ+396) = 1;           %ȷ����Թ۲���ǰ396���۲��Ȩ����
% tao(ts) = 1;                        %396�Ժ�Ϊ�ɻ���Ȩ����
% % tao(1:85,1) = ones(1,85);
% tao(1:numJ,1) = ones(1,numJ);
% for i = 86:n
%     if P(i,i) <= 10^-5
%         tao(i) = 1;
%     end
%     if tao(i) == 0
%         P(i,i) = 10^-8;
%         tao(i) = 1;
%     end
% end
% for i = 145:n
%     if P(i,i) <= 10^-5
%         tao(i) = 1;
%     end
%     if tao(i) == 0
%         P(i,i) = 10^-8;
%         tao(i) = 1;
%     end
% end
% tao(2023+85:end,1) = ones(1,1114);
% down_p = diag(tao.*diag(P));        %�õ����յ�Ȩ����
% num0 = 0;
% for i = 1:n 
%     if down_p(i,i) == 10^-8
%        num0 = num0 + 1;
%     end
% end
% num0
%***********************************************
%        for i = 1:n
%            tao(i) = V(i)/sqrt(PV(i,i));    
%        end
% %      mad_sigma = median(abs(tao-ones(n,1)*median(tao)))*1.4826;
%        mad_sigma = mad(abs(tao),1)*1.4826;
%        tao = tao/mad_sigma;
%        w = min(1,tune./abs(tao));           %�������HuberȨ���� 
%        ts = numJ + weiN;
%        w(ts) = 1;
%        down_p = diag(w'.*diag(P)); 
%mad_sigma = mad(abs(tao),1)*1.4826;
%tao = tao/mad_sigma;
%w = min(1,tune./abs(tao));                   %�������HuberȨ���� 
%w = 1./w;
%down_Q = diag([w])*Q*diag([w]);
%***********************************************