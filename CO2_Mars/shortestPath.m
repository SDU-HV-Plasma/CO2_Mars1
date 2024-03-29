function minWeightMatrix=shortestPath(G,nodeNum)
 
minWeightMatrix=zeros(nodeNum,nodeNum);
for i=1:nodeNum %��ȡÿ�����㵽�����ڵ����·��Ȩֵ��Dijkstra Դ�����
    D=G(i,:);   %��ʼ���ڵ㵽�����ھӽڵ��Ȩֵ
    D(D==0)=inf;
    tag=zeros(1,nodeNum); %0 Ϊ��Ӧλ�ýڵ�δ������·��
    tag(i)=1;
    while sum(tag,2)~=nodeNum
        noTag=find(tag==0); %��δ�����̾���Ľڵ�
        minValue=min(D(noTag));
        if minValue~=inf
            tempIndex=find(D==minValue); %�ҳ���ѡ�ڵ㣬����ֵ
            tag(tempIndex)=1;
            % ����ѡ�ڵ�ȥ��
            for candidateNode=1:length(tempIndex)
                noTag(noTag==tempIndex(candidateNode))=[];
            end
            % �޸�Դ�㵽����δ��ǵĶ����ľ���
            for j=1:length(noTag) %������δ�����̾���Ľڵ���и���
                minTemp=D(noTag(j));
                for k=1:length(tempIndex)
                    if minTemp>(minValue+G(tempIndex(k),noTag(j)))&&G(tempIndex(k),noTag(j))~=0
                        minTemp=minValue+G(tempIndex(k),noTag(j));
                    end
                end
                D(noTag(j))=minTemp;
            end            
        else
            break;
        end
    end
    minWeightMatrix(i,:)=D;
end
minWeightMatrix(isinf(minWeightMatrix))=0;
minWeightMatrix=minWeightMatrix-diag(diag(minWeightMatrix));


