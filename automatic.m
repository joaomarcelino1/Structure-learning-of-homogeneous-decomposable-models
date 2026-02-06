%Compute all results in the article
nvar=[5,7,9,11];

ndim=6;
nmodels=35;
        
SHDHCGP = cell(1,length(nvar));
      
SHDCSZL = cell(1,length(nvar));


for i = 1:length(nvar)
          
    SHDHCGP{i} = cell(1,ndim);
    SHDCSZL{i} = cell(1,ndim);

    for j=1:ndim
        SHDHCGP{i}{j}=zeros(nmodels,3);
        SHDCSZL{i}{j}=zeros(nmodels,5);
    end
end


TimeHCGP=SHDHCGP;
TimeCSZL=SHDCSZL;

%ii in nvar
%j number of reps
%k in dim 

dim=[100,500,1000,2000,4000,6000];

for ii=1:length(nvar)
    for j= 1:35
            
            %read edges
            filename = sprintf("edges%d_%d.xlsx", nvar(ii),j);
            T = readtable(filename);
        
            subsets2 = table2array(T(:,1:2));
            
            subsets2=subsets2+1;
            subsets2 = sort(subsets2, 2);
            N=nvar(ii);
            n = max(subsets2(:));

            %Combinatorial indexes
            indK=nchoosek(N,2);
            for i=2:N-1
                indK(i)=indK(i-1)+nchoosek(n,i+1);
            end 
            
            % chimset of truth-graph
            chtrue=zeros(nchoosek(N,2),1);
            for i=1:length(subsets2)
                S=subsets2(i,:);
                soma = setind2new(S,N,indK);
                chtrue(soma)=1;
            end 

        for k=1:length(dim)

            %Adjacency matrix
            n = max(subsets2(:));
                        
            %Read sample
            filename = sprintf("amostras%d_%d_%d.xlsx", nvar(ii),j,dim(k));
            D = readtable(filename);
            
            D= table2array(D);
            D(1,:) = [];
                        
            %HCGP1
            tic
            sol1 = HCGP2(D);
            time1=toc;
            
            ch1=sol1.x(1:nchoosek(N,2));
            
            SHDHCGP{ii}{k}(j,1) = SHD(ch1,chtrue);
            TimeHCGP{ii}{k}(j,1) =time1;

            %HCGP2
            tic
            sol2 = HCGP2(D);
            time2=toc;
            
            ch2=sol2.x(1:nchoosek(N,2));
            
            SHDHCGP{ii}{k}(j,2) = SHD(ch2,chtrue);
            TimeHCGP{ii}{k}(j,2) =time2;

            %HCGP3
            tic
            sol3 = HCGP3(D);
            time3=toc;
            
            ch3=sol3.x(1:nchoosek(N,2));
            
            SHDHCGP{ii}{k}(j,3) = SHD(ch3,chtrue);
            TimeHCGP{ii}{k}(j,3) =time1;
            if k==5
                CSZL
                for kappa=3:4
                    if kappa < nvar(ii)
                        tic
                        solC = newclique3(D,kappa);
                        timeCSZ = toc;
                        
                        TimeCSZL{ii}{k}(j,kappa-2) = timeCSZ;
                        chlearned = solC.x(1:nchoosek(N,2));
                        SHDCSZL{ii}{k}(j,kappa-2) = SHD(chlearned,chtrue);
                    end 
                end   
            end
        end    
    end     
end