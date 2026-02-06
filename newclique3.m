function sol = newclique3(D,kappa)
%INPUT:
%Clique size limit approach using HCGP3
%kappa :: clique size limit
%D : dataset

%OUTPUT: characteristic imset
tic
s=size(D);

nvar=s(2);
n=nvar;

last=-nvar;
for i=1:kappa+1 
    old=last;
    last=last+nchoosek(nvar,i);
end

%Calculate combinations
indK=nchoosek(nvar,2);

for i=2:kappa
    indK(i)=indK(i-1)+nchoosek(n,i+1);
end 

setS=gerarsetS2(nvar,kappa,indK(kappa));
setS=setS';

res=computeSubsetFrequenciesclique(D,kappa);
res = sortSubsetResultsclique(res,nvar);

[~, n2]= size(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation and writing of the objective function to a .txt file%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation of multiinformation function
m=zeros(length(res),1);

for i=n2:length(res)
    
    v=res{i}.subset;
        
    for j=1:length(res{i}.outcomes)
        dim=length(res{i}.subset);
        prod=1;
        for k=1:dim
            prod=prod*res{v(k)+1}.relativeFrequencies(res{i}.outcomes(j,k)+1);
        end 

        m(i)=m(i)+length(D)*res{i}.relativeFrequencies(j)*log2(res{i}.relativeFrequencies(j)/prod);
    end 
    
end 

%Calculate objective function
n = numel(res);  
result = zeros(n, 1);  

% Run set A
for i = n2+2:length(res)
    
    A = res{i}.subset; 
    card_A = numel(A);
    soma = 0;

    % Iterate over the sets B_i that precede A and have at least 2 elements
    for j = n2+2:i
        B = res{j}.subset;
        if numel(B) >= 2 && all(ismember(B, A)) % Verifica se B está contido em A
           soma = soma + ((-1)^(card_A - numel(B)))*(m(j))-log(s(2))*(2^numel(B)-numel(B)-1)/2;
        end
    end
    
    result(i) = soma;
end

%open lp file
fid = fopen('newclique3.lp', 'w');

%%%%%%%%%%%%%%%%%%%%
%Objective function%
%%%%%%%%%%%%%%%%%%%%

%Coeficientes função objetivo nao calculados 
% objind=zeros(length(res),1);

% Escrever a palavra "maximize"
fprintf(fid, 'Maximize \n ');
alfa=1/2;
objind=result;

%MDL
objind=objind(2+nvar:end)-alfa*log10(length(D));

fprintf(fid,'  obj: ');
for i = 1:length(objind)
    coef = round(objind(i), 4); 
    if coef < 0
        fprintf(fid, '- %.4f x%d ', abs(coef), i);
    else
        if i == 1
            fprintf(fid, '%.4f x%d ', coef, i);
        else
            fprintf(fid, '+ %.4f x%d ', coef, i);
        end
    end
    if mod(i,3)==0
        fprintf(fid,'\n');
        fprintf(fid,'        ');
    end
end

fprintf(fid, '\n'); fprintf(fid,'\n');
fprintf(fid, 'Subject To \n');

r=0; 

for i=indK(1)+1:indK(kappa)
    positions=calcsubsets(setS{i},nvar,indK);

    if i <=indK(kappa-1)
        
        restr1=[' r',num2str(r+3),': ', num2str(-length(setS{i})+2),' x',num2str(i)];
    
        %Restriction 14
        for j=1:length(positions)
            restr1=[restr1,' + x',num2str(positions(j))];
    
            %Restriction 15
            restr2=[' r',num2str(r),': x',num2str(i),' - x',num2str(positions(j)),' <= 0'];
            
            fprintf(fid, ' %s\n', restr2);
            r=r+1;
        end 
    
        restr1=[restr1,' <= 2'];
        fprintf(fid, ' %s\n', restr1);
        r=r+1;
    else

        restr1=[' r',num2str(r),': '];

        for j=1:length(positions)
            restr1=[restr1,' + x',num2str(positions(j))];
        end

        fprintf(fid, ' %s <= 2\n', restr1);
        r=r+1;

    end 
    
end 

%Create Ind2
n=nvar;
Ind2=zeros(n,n);
for i=1:nchoosek(n,2)
    aux=setS{i};
    Ind2(aux(1),aux(2))=i;
    Ind2(aux(2),aux(1))=i;
end 

%Create Ind3
Ind3=zeros(n,n,n);
for i=nchoosek(n,2)+1:nchoosek(n,2)+nchoosek(n,3)
    aux=setS{i};
    Ind3(aux(1),aux(2),aux(3))=i;
    Ind3(aux(2),aux(1),aux(3))=i;
    Ind3(aux(3),aux(1),aux(2))=i;
    Ind3(aux(1),aux(3),aux(2))=i;
    Ind3(aux(2),aux(3),aux(1))=i;
    Ind3(aux(3),aux(2),aux(1))=i;
end

%Permutation of size 4
P=[1 2 3 4;
    1 2 4 3;
    1 3 2 4;
    1 3 4 2;
    1 4 2 3;
    1 4 3 2;
    2 1 3 4;
    2 1 4 3;
    2 3 1 4;
    2 4 1 3;
    3 1 2 4;
    3 2 1 4
    ];

for i=indK(2)+1:indK(3)
    v=setS{i}; 

    if kappa == 3
        for j=1:12
                restr=[' r',num2str(r),': x',num2str(Ind2(v(P(j,1)),v(P(j,2)))),' + x',num2str(Ind2(v(P(j,1)),v(P(j,4)))),' + x',num2str(Ind2(v(P(j,2)),v(P(j,3)))),' + x',num2str(Ind2(v(P(j,3)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,2)),v(P(j,3)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,2)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,3)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,2)),v(P(j,3)),v(P(j,4)))),' <= 2'];
                fprintf(fid, '%s\n',restr);
            r=r+1;
        end

    else

        for j=1:12
            restr=[' r',num2str(r),': x',num2str(Ind2(v(P(j,1)),v(P(j,2)))),' + x',num2str(Ind2(v(P(j,1)),v(P(j,4)))),' + x',num2str(Ind2(v(P(j,2)),v(P(j,3)))),' + x',num2str(Ind2(v(P(j,3)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,2)),v(P(j,3)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,2)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,1)),v(P(j,3)),v(P(j,4)))),' - x',num2str(Ind3(v(P(j,2)),v(P(j,3)),v(P(j,4)))),' + 2 x',num2str(i),' <= 2'];
            fprintf(fid, '%s\n',restr);
            r=r+1;
        end 
    end
end 

fprintf(fid,'\n');

%Binary
fprintf(fid, 'Binaries \n ');
for i=1:old
    fprintf(fid, ' x%d', i);
    if mod(i,10)==0
        fprintf(fid,'\n ');
    end 
end

fprintf(fid,'\n\nEnd');

%Close file
fclose(fid);
toc

model=gurobi_read('newclique3.lp');

%Solve ILP
sol=gurobi(model);

end 

