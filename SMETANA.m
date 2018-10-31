function alignment=SMETANA(Net_id_list,input_folder,id_flag,out_file)
%  SMETANA multiple network alignment algorithm
%   alignment=SMETANA(Net_id_list,input_folder,out_file)
%
%  SMETANA finds the multiple alignment of a set of (biological) networks
%
%   input   Net_id_list     the input networks ids list
%           input_folder    the foder, where input files are placed
%           id_flag         This flag is used to expedite the reading
%                           process of inout files. The flag that indicates
%                           whether the nodes are in numeric format
%                           (id_flag=1) or not (id_flag=0). If nodes are in
%                           the formatof "species id+number", the id_flag
%                           should be 1, unless it is zero. For instance if
%                           nodes are as a0,a1,b1,b3,... the id_flag is 1.
%           out_file        output file address where the alignment result
%                           will be written there
%   output  alignment       Alignment result including the obtained aligned
%                           nodes
%
%   Example:
%        alignment=SMETANA({'dm','hs','sc'},'test',1,'output.txt')
%        alignment=SMETANA({'a','b','c'},'test',1,'output.txt')
%
%   Input files format:
%       The Net_id_list provides the list of networks (spicies) names. For example
%       here we have three species named as 'a', 'b', and 'c'. For these
%       specises we need the following files:
%
%       - Network files: a.net, b.net, c.net
%               These are tab-separated files that list the (undirected)
%               interactions in each network. For instanse, a.net may be as
%               follows:
%
%               a1	a2
%               a3	a1
%               a4	a2
%               a2	a3
%
%               Network files can also include the edge weights. In tha
%               case we may have a third column:
%
%               a1	a2  0.5
%               a3	a1  0.2
%               a4	a2  0.8
%               a2	a3  0.9
%
%       - Similarity scores files: a-b.sim, a-c.sim, b-c.sim
%
%               These files consist of the list of similarity scores
%               between nodes of different species. We have used BLAST
%               bit-score for our tests, but other scores such as
%               log(E-value) can also be used.
%
%               IMPORTANT note 1:
%               The similarity filename should have the species names in
%               lexicographic order, i.e., a-b.sim is expected, not
%               b-a.sim.
%
%               IMPORTANT note 2:
%               For a-b.sim file, the first column should contain node from
%               speices a and the second column should contain IDs from
%               species b.
%
%               As an example, a-b.sim may be as follows:
%
%               a1	b1	153
%               a1	b3	55
%               a1	b7	49
%               a2	b3	444
%               a3	b3	211
%               a3	b4	122
%               a4	b5	251
%               a4	b8	71
%
%  Output files format:
%       The Output files includes the list of aligned nodes. For instance:
%
%               a4 b5 c4
%               a1 b1 b7 c1
%               a2 b3 c2
%               a3 b4 c3
%
% For more information on the algorithms, please see:
%
% Sahraeian, S.M.E., Yoon, B.J. (2013), SMETANA: Accurate and Scalable
% Algorithm for Probabilistic Alignment of Large-Scale Biological Networks,
% submitted.
%
% By Sayed Mohammad Ebrahim Sahraeian and Byung-Jun Yoon
% Contact: bjyoon@ece.tamu.edu
%-------------------------------------------------------------------------


M=length(Net_id_list);
nets=cell(1,M);
nodes=cell(1,M);
edges=cell(1,M);
Sims=cell(M,M);
Net_id_list=lower(Net_id_list);
for i=1:M
    fid = fopen(strcat(input_folder,'/',Net_id_list{i},'.net'));
    if id_flag==0
        nets{i} = textscan(fid,'%s\t%s\t%s', 1000000);
    else
        nets{i} = textscan(fid,strcat(Net_id_list{i},'%d\t',Net_id_list{i},'%d\t%s'), 1000000);
    end;
    fclose(fid);
    if isempty(nets{i}{3}{1})
        nets{i}(3)=[];
    end;
end;

for i=1:M
    for j=i+1:M
        fid = fopen(strcat(input_folder,'/',Net_id_list{i},'-',Net_id_list{j},'.sim'));
        if id_flag==0
            Sims{i,j} = textscan(fid,'%s\t%s\t%f', 10000000);
        else
            Sims{i,j}= textscan(fid,strcat(Net_id_list{i},'%d\t',Net_id_list{j},'%d\t%f'), 1000000);
        end;
        fclose(fid);
    end;
end;

if (id_flag==0)
    for i = 1: M
        nodes{i} = unique([nets{i}{1}; nets{i}{2}]);
        for j = i + 1: M
            nodes{i} = unique([nodes{i}; Sims{i, j}{1}]);
        end
        
        for j = 1: i - 1
            nodes{i} = unique([nodes{i}; Sims{j, i}{2}]);
        end
        ids = cell(1, M);
        for j = 1: M
            ids{j} = [];
        end
        
        for j = i + 1: M
            [~, tempList] = ismember(Sims{i, j}{1}, nodes{i});
            Sims{i, j}{1} = {tempList};
        end
        
        for j = 1: i - 1
            [~, tempList] = ismember(Sims{j, i}{2}, nodes{i});
            Sims{j, i}{2} = {tempList};
        end
    end
    
    for i = 1: M
        edges{i} = zeros(length(nets{i}{1}), 3);
        [~, edges{i}(:, 1)] = ismember(nets{i}{1}, nodes{i});
        [~, edges{i}(:, 2)] = ismember(nets{i}{2}, nodes{i});
        
        if length(nets{i}) == 2
            edges{i}(:, 3) = 1;
        else
            edges{i}(:, 3) = str2double(nets{i}{3});
        end
    end
    
    for i = 1: M
        for j = i + 1: M
            Sims{i, j}=[cell2mat(Sims{i, j}{1}), cell2mat(Sims{i, j}{2}), (Sims{i, j}{3})];
        end
    end
else
    nnodes=cell(1,M);
    for i=1:M
        nnodes{i}=max([nets{i}{1};nets{i}{2}]);
        for j=i+1:M
            nnodes{i}=max([nnodes{i};Sims{i,j}{1}]);
        end;
        for j=1:i-1
            nnodes{i}=max([nnodes{i};Sims{j,i}{2}]);
        end;
        
    end;
    
    for i=1:M
        nodes{i}=1:nnodes{i};
    end;
    for i=1:M
        if length(nets{i})==2
            edges{i}=[nets{i}{1},nets{i}{2},ones(length(nets{i}{1}),1)];
        else
            edges{i}=[nets{i}{1},nets{i}{2},str2double(nets{i}{3})];
        end;
    end;
    
    
    for i=1:M
        for j=i+1:M
            Sims{i,j}=[Sims{i,j}{1},Sims{i,j}{2},Sims{i,j}{3}];
        end;
    end;
    
    
end;

%---------Perform Alignment using SMETANA-------------------------------
tic;


%---initializtion
M=length(edges);
L=zeros(1,M);
for i=1:M
    L(i)=length(nodes{i});
end;
G=cell(1,M);
for i=1:M
    ee = double(edges{i});
    Q = sparse(ee(:, 1), ee(:, 2), ee(:, 3), L(i), L(i));
    indices  = sub2ind(size(Q), ee(:, 2), ee(:, 1));
    Q(indices) = ee(:, 3);
    G{i} = Q;
end;
S=cell(M,M);
for i=1:M
    for j=i+1:M
        S{i,j}=sparse(double(Sims{i,j}(:,1)),double(Sims{i,j}(:,2)),double(Sims{i,j}(:,3)),L(i),L(j));
    end;
end;

%---Removing non-homologue nodes
J=cell(1,M);
for i=1:M
    J{i}=1:L(i);
end;
Rj=J;
for i=1:M
    for j=i+1:M
        Oj=find(sum(S{i,j})==0);
        Oi=find(sum(S{i,j}')==0);
        Rj{i}=intersect(Rj{i},Oi);
        Rj{j}=intersect(Rj{j},Oj);
    end;
end;
[GG,SS,J,LL]=update(L,G,S,Rj);


%Estimation of probabilistic node correspondence scores through semi-Markov random walk
C=ComputeCScore(GG,SS);

%Perform Intra-network probabilistic consistency transformation
C2=PCT_i(GG,C,LL);

%Perform Cross-network probabilistic consistency transformation
C2=PCT(C2,LL);

%Alignment construction
AL=AlignmentConstruction(C2,J,LL);
t2=toc;

sprintf('Elapsed time for Alignment: %f seconds',t2)
%Print the alignment result
fod = fopen(out_file,'w');
alignment=cell(1,length(AL));
for i=1:length(AL)
    alignment{i}=cell(1,size(AL{i},1));
    for j=1:size(AL{i},1)
        if id_flag==0
            fprintf(fod,'%s ',nodes{AL{i}(j,1)}{AL{i}(j,2)});
            alignment{i}{j}=nodes{AL{i}(j,1)}{AL{i}(j,2)};
        else
            fprintf(fod,'%s%d ',Net_id_list{AL{i}(j,1)},AL{i}(j,2));
            alignment{i}{j}=strcat(Net_id_list{AL{i}(j,1)},num2str(AL{i}(j,2)));
        end;
        
    end;
    fprintf(fod,'\n');
end;
fclose(fod);

end
%----------------------------------------------------------\



%----------------------------------------------------------
%Perform Cross-network probabilistic consistency transformation
%----------------------------------------------------------
function C2=PCT(C,L)
M=size(C,2)-1;
C2=[];
H=eye(M);
for i=1:M
    H(i,i)=1;
end;
for i=1:M
    for j=i+1:M
        J=cell(1,2);
        J{1}=1:L(i);
        J{2}=1:L(j);
        xy=find(C(:,i)>0 & C(:,j)>0);
        
        MMM=AlignmentConstruction(C(xy,[i,j,M+1]),J,[L(i),L(j)]);
        CC=sparse(C(xy,i),C(xy,j),C(xy,M+1),L(i),L(j));
        p=0;
        cnt=0;
        for ii=1:length(MMM)
            f1=find(MMM{ii}(:,1)==1);
            f2=find(MMM{ii}(:,1)==2);
            p=p+sum(sum(CC(MMM{ii}(f1,2),MMM{ii}(f2,2))));
            cnt=cnt+length(f1)*length(f2);
        end;
        H(i,j)= p/cnt;
        H(j,i)= H(i,j);
    end;
end;

for i=1:M
    for j=i+1:M
        P=sparse(L(i),L(j));
        for k=1:M
            if k==i
                Pxz=speye(L(k));
            else
                xz=find(C(:,i)>0 & C(:,k)>0);
                Pxz=sparse(C(xz,i),C(xz,k),C(xz,M+1),L(i),L(k));
            end;
            if k==j
                Pzy=speye(L(k));
            else
                zy=find(C(:,k)>0 & C(:,j)>0);
                Pzy=sparse(C(zy,k),C(zy,j),C(zy,M+1),L(k),L(j));
            end;
            
            P=P+Pxz*Pzy*H(i,k)*H(k,j)/sum(H(i,:).*H(:,j)');
        end;
        xy=find(C(:,i)>0 & C(:,j)>0);
        Pxy=sparse(C(xy,i),C(xy,j),C(xy,M+1),L(i),L(j));
        ll=length(xy);
        [cc,ii]=sort(nonzeros(P(Pxy==0)),'descend');
        mmm=min(round(ll*.01),length(cc));
        if mmm>0
            th=cc(min(round(ll*.01),length(cc)));
        else
            th=0;
        end;
        P=P.*(Pxy>0)+P.*(Pxy==0).*(P>th);
        [ind3,ind4,dd2]=find(P/M);
        z=-ones(length(dd2),M+1);
        z(:,i)=ind3;
        z(:,j)=ind4;
        z(:,M+1)=dd2;
        C2=[C2;z];
    end;
end;
end


%----------------------------------------------------------
%Perform Intra-network probabilistic consistency transformation
%----------------------------------------------------------
function C2=PCT_i(G,C,L)
M=size(C,2)-1;
C2=[];
for i=1:M
    for j=i+1:M
        P=sparse(L(i),L(j));
        xy=find(C(:,i)>0 & C(:,j)>0);
        Pxy=sparse(C(xy,i),C(xy,j),C(xy,M+1),L(i),L(j));
        
        aaa=.9;
        P=aaa*Pxy+(1-aaa)*G{i}*Pxy*G{j}';
        ll=length(xy);
        [cc,ii]=sort(nonzeros(P(Pxy==0)),'descend');
        mmm=min(round(ll*.01),length(cc));
        if mmm>0
            th=cc(min(round(ll*.01),length(cc)));
        else
            th=0;
        end;
        P=P.*(Pxy>0)+P.*(Pxy==0).*(P>th);
        [ind3,ind4,dd2]=find(P);
        z=-ones(length(dd2),M+1);
        z(:,i)=ind3;
        z(:,j)=ind4;
        z(:,M+1)=dd2;
        C2=[C2;z];
    end;
end;
end






%----------------------------------------------------------
%Estimation of probabilistic node correspondence scores
%through semi-Markov random walk
%----------------------------------------------------------
function C=ComputeCScore(G,S)

M=length(G);
L=zeros(1,M);
for i=1:M
    L(i)=size(G{i},1);
end;
Pi=cell(1,M);
for i=1:M
    f=find(sum(G{i})==0);
    for ii=1:length(f)
        G{i}(f(ii),f(ii))=1;
    end;
    qq=sum(G{i})';
    [h,j,k]=find(G{i});
    for ii=1:length(h)
        G{i}(h(ii),j(ii))=k(ii)/qq(h(ii));
    end;
    
    Pi{i}=sparse(Steady(G{i},repmat(1/L(i),1,L(i))));
    
    C=[];
    cnnt=0;
end;

for i=1:M
    for j=i+1:M
        [ind1,ind2,dd]=find(S{i,j});
        F=sparse(L(i),L(j));
        
        for k=1:length(dd)
            F(ind1(k),ind2(k))=Pi{i}(ind1(k))*Pi{j}(ind2(k));
        end;
        
        
        FF=F.*S{i,j};
        FF=FF/sum(sum(FF));
        for ttt=1:1
            s1=spdiags(1./sum(FF)',0,size(FF,2),size(FF,2));
            s2=spdiags(1./sum(FF')',0,size(FF,1),size(FF,1));
            FF=(FF*s1+s2*FF)/2;
        end;
        
        [ind3,ind4,dd2]=find(FF);
        z=-ones(length(dd2),M+1);
        z(:,i)=ind3;
        z(:,j)=ind4;
        z(:,M+1)=dd2;
        C=[C;z];
        cnnt=size(C,1);
    end;
end;
end


%----------------------------------------------------------
%Computing steady state distribution
%----------------------------------------------------------
function Pi=Steady(A,p)
err=1;
J=p';
cnt=0;
while ((err>.00001)&(cnt<100))
    J0=J;
    JJ=(J'*A)';
    J=JJ/norm(JJ,1);
    err=norm(J0-J,2);
    cnt=cnt+1;
    if (cnt>1000)
        aaa=1;
    end;
end;
Pi=J;

end


%----------------------------------------------------------
%Update the networks after Removing non-homologue nodes
%----------------------------------------------------------
function [GG,SS,J,LL]=update(L,G,S,Rj)
M=length(L);
J=cell(1,M);
GG=cell(1,M);
for i=1:M
    J{i}=1:L(i);
    J{i}(Rj{i})=[];
    GG{i}=sparse(length(J{i}),length(J{i}));
    GG{i}=G{i}(J{i},J{i});
end;
SS=cell(M,M);
for i=1:M
    for j=i+1:M
        SS{i,j}=sparse(length(J{i}),length(J{j}));
        SS{i,j}=S{i,j}(J{i},J{j});
    end;
end;

LL=zeros(1,M);
for i=1:M
    LL(i)=size(GG{i},1);
end;

end



%----------------------------------------------------------
%Computing Semi-Markov correspondence scores
%----------------------------------------------------------
function MMM=AlignmentConstruction(C,J,L)
M=length(J);

[s,y]=sort(C(:,M+1),'descend');

CC=cell(M,M);
for i=1:M
    CC{i,i}=sparse(L(i),L(i));
end;

for i=1:M
    for j=i+1:M
        f=find(C(:,i)>0 & C(:,j)>0);
        CC{i,j}=sparse(C(f,i),C(f,j),C(f,M+1),L(i),L(j));
        CC{j,i}=CC{i,j}';
    end;
end;


MM=cell(1,M);
for i=1:M
    MM{i}=sparse(1,L(i));
end;
kkkk=.8;
k=1;
g_max=1;
max_l=ones(1,M)*10;
EC={};
missed=[];
ssss=[];
while (k<=length(y))
    h=C(y(k),1:M);
    if C(y(k),M+1)<0
        k=k+1;
        continue;
    end;
    if (prod(h)~=0)
        F=find(h~=-1);
        f1=F(1);
        f2=F(2);
        h1=h(f1);
        h2=h(f2);
        c1=MM{f1}(h1);
        c2=MM{f2}(h2);
        if ((c1==0)&(c2==0))
            MM{f1}(h1)=g_max;
            MM{f2}(h2)=g_max;
            g_max=g_max+1;
            EC=[EC,[f1,h1;f2,h2]];
        elseif ((c1>0)&(c2==0))
            if length(find(EC{c1}(:,1)==f2))>=max_l(f2)
                k=k+1;
                continue;
            end;
            if ~isempty((find(EC{c1}(:,1)==f2)))
                others=find((EC{c1}(:,1)~=f2));
                own=find((EC{c1}(:,1)==f2));
                m_1=zeros(1,length(own));
                m_2=0;
                for j=1:length(others)
                    jj=others(j);
                    for i=1:length(own)
                        ii=own(i);
                        m_1(i)=m_1(i)+CC{f2,EC{c1}(jj,1)}(EC{c1}(ii,2),EC{c1}(jj,2));
                    end;
                    m_2=m_2+CC{f2,EC{c1}(jj,1)}(h2,EC{c1}(jj,2));
                end;
                
                if m_2<(kkkk*mean(m_1))
                    k=k+1;
                    continue;
                end;
            end;
            EC{c1}=[EC{c1};f2,h2];
            MM{f2}(h2)=MM{f1}(h1);
        elseif ((c1==0)&(c2>0))
            if length(find(EC{c2}(:,1)==f1))>=max_l(f1)
                k=k+1;
                continue;
            end;
            if ~isempty((find(EC{c2}(:,1)==f1)))
                others=find((EC{c2}(:,1)~=f1));
                own=find((EC{c2}(:,1)==f1));
                m_1=zeros(1,length(own));
                m_2=0;
                for j=1:length(others)
                    jj=others(j);
                    for i=1:length(own)
                        ii=own(i);
                        m_1(i)=m_1(i)+CC{f1,EC{c2}(jj,1)}(EC{c2}(ii,2),EC{c2}(jj,2));
                    end;
                    m_2=m_2+CC{f1,EC{c2}(jj,1)}(h1,EC{c2}(jj,2));
                end;
                
                if m_2<(kkkk*mean(m_1))
                    k=k+1;
                    continue;
                end;
            end;
            EC{c2}=[EC{c2};f1,h1];
            MM{f1}(h1)=MM{f2}(h2);
        elseif (c1~=c2)
            ec_temp=[EC{c1};EC{c2}];
            flag=0;
            [tt1,htt1]=hist(ec_temp(:,1),1:M);
            if max(tt1)>1
                k=k+1;
                continue;
            end;
            
            for i=1:M
                MM{i}(find(MM{i}==c1 | MM{i}==c2))=g_max;
            end;
            EC=[EC,[EC{c1};EC{c2}]];
            EC{c1}=[];
            EC{c2}=[];
            g_max=g_max+1;
        end;
    end;
    k=k+1;
end;



MMM=cell(g_max,1);
for j=1:g_max
    for i=1:M
        ff=find(MM{i}==j);
        if ~isempty(ff)
            MMM{j}=[MMM{j};ones(length(ff),1)*i,J{i}(ff)'];
        end;
    end;
end;
cnt=1;
while cnt<=length(MMM)
    if isempty(MMM{cnt})
        MMM(cnt)=[];
    else
        cnt=cnt+1;
    end;
end;

end






