function display_rates(nsp,species,number_of_reactions)

h=(get(0,'children'));
ind=find(h~=1);
close(ind);

data=guidata(1);  
            
%Density evolution
figure
scs=get(0,'ScreenSize');
width=.7*scs(3);
height=.25*scs(4);
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
h=semilogy(data.hist_t,data.hist_nall,'linewidth',1);
hold on
legend(h,data.species,'location','eastoutside','numcolumns',2)
xlabel('Time');
ylabel('Density');

% %Density of all species
% E_levels=[0,11.548,11.624,11.723,11.828,12.907,13.116,13.273,13.295,13.328,13.48,13.884,13.994,14.09,14.229,14.252,14.304,14.509,14.69,14.792,14.906,14.976,15.028,15.083,15.153,15.205,15.215,15.282,15.324,15.347];
% figure
%  %h=semilogy(1:length(data.species),data.hist_n(end,:));
%  h=semilogy(E_levels,data.hist_n(end,1:30),'r');
%     xlabel('Energy of Level (eV)');
%    
% hold on

% Each single species
for i=1:nsp
    figure(i+2)
    set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
    h=semilogy(data.hist_t,data.hist_nall(:,i),'linewidth',1);
    hold on
    xlabel('Time');
    ylabel('Density');
    title(species{i});
end

%Te evolution
figure
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
h=plot(data.hist_t,data.hist_Te,'linewidth',1);
hold on
xlabel('Time');
ylabel('Electron temperature (eV)');
title('Te');

figure
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
plot(data.hist_t,data.hist_dS_ape,'r',data.hist_t,data.hist_dS_el,'g',data.hist_t,data.hist_dS_epro,'b',data.hist_t,data.hist_dS_rea,'m',data.hist_t,data.hist_dS_wall,'k');
legend('��λ������������','������ײ����ɢʧ','�����ܶȱ仯����ĵ����¶ȱ仯','��Ӧ��ѧ�ܴ����������仯','�����ڱ��ϵ���ʧ�����������仯');
xlabel('Time');
ylabel('��λʱ���ڵ����������eV��');

%Tv evolution
figure
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
h=plot(data.hist_t,data.hist_Tv*11600,'linewidth',1);
hold on
xlabel('Time');
ylabel('Temperature of vibrations (K)');
title('Tv');

%Particle balance for each species
figure
width=.75*scs(3);
height=.75*scs(4);
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
for isp=1:data.nsp
    subplot(ceil(sqrt(data.nsp)),ceil(data.nsp/ceil(sqrt(data.nsp))),isp);
    bar(data.X(isp,:));
    title(data.species(isp));
    YL=max(abs(data.X(isp,:)));
    xlim([0 number_of_reactions+1]);
    xlabel('Reaction number')
    set(gca,'xtick',0:5:number_of_reactions+1);
    ylabel('Frequency (s^-^1)');
    %Relative reaction contribution for each species
    gen_reactions=find(data.X(isp,1:number_of_reactions)>0);
    los_reactions=find(data.X(isp,1:number_of_reactions)<0);
    reaction_contribution(isp,gen_reactions)=100*data.X(isp,gen_reactions)/sum(data.X(isp,gen_reactions));
    reaction_contribution(isp,los_reactions)=100*data.X(isp,los_reactions)/sum(data.X(isp,los_reactions));
end
figure
width=.5*scs(3);
height=.65*scs(4);
set(gcf,'position',[(scs(3)-width)/2 (scs(4)-height)/2 width height]);
subplot(2,1,1)
bar(sum(reaction_contribution)/max(sum(reaction_contribution))*100)
title('Reaction relevance');
xlabel('Reaction number')
ylabel('Relative Relevance (a.u.)');%���ȼ���ĳ����Ӧ��ĳ���Ӳ�����Ӱ�������з�Ӧ�Ը����Ӳ�����Ӱ��֮�͵ı�ֵ��Ȼ��Ը÷�Ӧ���������ӵĸý��������õ���ֵ����ͨ���ٷֱ���Ӽ��㣬Ϊ��ԡ�
subplot(2,1,2)%��ĳ��Ӧ�Ե�һ���ӵ�Ӱ�����ռ��90%�������ֵ��ܴ���Ϊ90%�ܴ󡣶��������ͨ����ֵ���㣬����������ڷ�Ӧ���ܱ仯����С����÷�Ӧ��Ӧ��90%��ʾ����ֵ΢����΢
bar(sum(abs(data.X(1:nsp,1:number_of_reactions)))/max(sum(abs(data.X(1:nsp,1:number_of_reactions))))*100)
xlabel('Reaction number')
ylabel('Absolute Relevance (a.u.)');%��ֱֻ�ۿ������Ӧ�Է�Ӧ��������������Ũ��Ӱ�����ֱ��ͨ����ֵ����

% figure
% aux=(sum(reaction_contribution)/max(sum(reaction_contribution))*100);
% aux=[(1:length(aux))' aux'];
% aux=sortrows(aux,-2);
% subplot(2,1,1)
% bar(aux(:,2))
% set(gca,'XTickLabel',[]);
% for i=1:length(aux)
%     text(i-.4,aux(i,2)+5,sprintf('R%d',aux(i,1)))
% end
% title('Reaction relevance');
% ylabel('Relative Relevance (a.u.)');
% 
% aux=(sum(abs(data.X))/max(sum(abs(data.X)))*100)
% aux=[(1:length(aux))' aux'];
% aux=sortrows(aux,-2);
% subplot(2,1,2)
% bar(aux(:,2))
% set(gca,'XTickLabel',[])
% for i=1:length(aux)
%     text(i-.4,aux(i,2)+5,sprintf('R%d',aux(i,1)))
% end
% ylabel('Absolute Relevance (a.u.)');

%Spectrum
% figure
% index=find(data.lambda>0);
% stem(data.lambda(index),data.R(index),'r.');
% xlabel('Wavelength')
% ylabel('Intensity')
% axis([300 1000 0 max(data.R(index))]);
