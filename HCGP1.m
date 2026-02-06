function sol = HCGP1(D)
%INPUT: Dataset D
%OUTPUT: Homogeneous decomposable model learned by HCGP^rel(1)

s=size(D);

%nvar:: number of variables 
nvar=s(2);

%Compute frequencies
setS=gerarsetS(nvar);
setS=setS';
setS=setS(nvar+2:end,1);

res=computeSubsetFrequencies(D);
res = sortSubsetResults(res);

n2=nvar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation and writing of the objective function to a .txt file%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation of the multiinformation function
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
n = numel(res);  % number of sets
result = zeros(n, 1);  % save the results

% Run set A
for i = n2+2:length(res)
    
    A = res{i}.subset; 
    card_A = numel(A);
    soma = 0;
    
    % Iterate over the sets B_i that precede A and have at least 2 elements
    for j = n2+2:i
        B = res{j}.subset;
        if numel(B) >= 2 && all(ismember(B, A)) % Check if B is a subset of A
           soma = soma + ((-1)^(card_A - numel(B)))*(m(j))-log(s(2))*(2^numel(B)-numel(B)-1)/2;
        end
    end
    
    result(i) = soma;
end

%open lpfile
fid = fopen('HCGP1.lp', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write objective function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "maximize"
fprintf(fid, 'Maximize \n ');

objind=result;

%BIC
objind=objind(2+nvar:end); 

fprintf(fid,'  obj: ');
for i = 1:length(objind)
    coef = round(objind(i), 4); 
    if coef < 0
        % for negative values write "-"
        fprintf(fid, '- %.4f x%d ', abs(coef), i);
    else
        % For positive values write "+", except in the first term
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

fprintf(fid, '\n');
fprintf(fid,'\n');

fprintf(fid, 'Subject To \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% WRITE RESTRICTIONS IN THE .txt file %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setind=zeros(2^nvar,1);
for i=1:length(setS)
    Si=setS{i};
    aux=0;
    for k=1:length(Si)
        aux=aux+2^(Si(k)-1);
    end
    setind(aux)=i;
end

value=2+nchoosek(nvar,2)+nvar;
newvalue=value-nvar-1;
newvalue2=value-nvar-1+nchoosek(nvar,3);

%%%%%%%%%%%%%%%%%%
%RESTRICTIONS (6)%
%%%%%%%%%%%%%%%%%%
r=0; 

for j=newvalue:length(setS)
    
    Sj=setS{j};
    aux=0;

    for i=1:length(Sj)
        aux =aux+ 2^(Sj(i)-1);
    end 
    
    restr=[' r',num2str(r),': ',num2str(-length(Sj)+2),' x',num2str(setind(aux))];

    for k=1:length(Sj)
        aux2=aux-2^(Sj(k)-1); 
        restr=[restr,' + x',num2str(setind(aux2))];               
    end

    restr=[restr,' <= 2'];

    fprintf(fid, ' %s\n', restr);
    r=r+1;
    
end

%%%%%%%%%%%%%%%%%%%
%RESTRICTIONS (7) %
%%%%%%%%%%%%%%%%%%%

for j=newvalue:length(setS)   
    Sj=setS{j};
    
    for v1=1:length(Sj)-1
        for v2=v1+1:length(Sj)             
              fprintf(fid, '  r%d: -x%d + x%d <= 0\n',r, setind(2^(Sj(v1)-1) + 2^(Sj(v2)-1)), j);           
              r=r+1;
        end
    end       
end

%%%%%%%%%%%%%%%%%%%
%RESTRICTIONS (11)%
%%%%%%%%%%%%%%%%%%%
newvalue3= nchoosek(nvar,2)+nchoosek(nvar,3)+nchoosek(nvar,4);

Perms=calcPerms(setS(newvalue2:newvalue3));

for j = 1:length(Perms)
    Sj=Perms{j};
    
    fprintf(fid, '  r%d: x%d + x%d + x%d - x%d - x%d + x%d <= 2\n',r,setind(2^(Sj(1)-1) + 2^(Sj(2)-1) ), setind(2^(Sj(2)-1) + 2^(Sj(3)-1) ),setind(2^(Sj(3)-1) + 2^(Sj(4)-1) ), setind(2^(Sj(1)-1) + 2^(Sj(2)-1)+2^(Sj(3)-1) ), setind(2^(Sj(2)-1) + 2^(Sj(3)-1)+2^(Sj(4)-1) ), setind(2^(Sj(1)-1) + 2^(Sj(2)-1) + 2^(Sj(3)-1)+2^(Sj(4)-1) ) );
    r=r+1;    

end 

fprintf(fid,'\n');

%Binary variables
fprintf(fid, 'Binaries \n ');
for i=1:length(setS)
    fprintf(fid, ' x%d', i);
    if mod(i,10)==0
        fprintf(fid,'\n ');
    end 
end

fprintf(fid,'\n\nEnd');

fclose(fid);

model=gurobi_read('HCGP1.lp');

%Solve ILP using Gurobi
sol=gurobi(model);

end