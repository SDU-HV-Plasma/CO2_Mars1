function findreaction(B,species)

global B species

m1=['�������һ����Ӧ�'];
m1=input(m1);
M1=strcmp(species,m1);
N1=unique([B{:,M1}]);
m2=['������ڶ�����Ӧ�'];
m2=input(m2);
M2=strcmp(species,m2);
N2=unique([B{:,M2}]);
m3=['�������������Ӧ�'];
m3=input(m3);
Thirdspecies=strcmp(m3,'null');
if(~Thirdspecies)
    M3=strcmp(species,m3);
    N3=unique([B{:,M3}]);
    N12=intersect(N1,N2);
    Reaction=intersect(N12,N3);
else
    Reaction=intersect(N1,N2);
end
Reaction=num2str(Reaction);

if(~isempty(Reaction))
    disp(strcat('����Ҫ��ķ�ӦΪ��',Reaction));
else
    disp('�޷��ҵ�����Ҫ��ķ�Ӧ��');
end

end

