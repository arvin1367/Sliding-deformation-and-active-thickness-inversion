function thickness_cal(ve,vn,vu,dx,f,lamda,up)
%%%%%%%%%%%%%%%%%%%%
% This code is used to invert the thickness of a landslide using the 
% three-dimensional flow velocity of the landslide. Gh=vz
% The material balance method is used. The code is written by Huiyuan Luo.
% Note: The CVX library (https://cvxr.com/cvx/) should be installed before use.
%ve:East-west deformation rate .tif
%vn:North-south deformation rate .tif 
%vu:Vertical deformation rate .tif
%dx:Grid resolution 
%f:rheological parameter 
%lamda:Regularization coefficient     recommend 0.1
%up:Maximum constraint. 
%%%%%%%%%%%%%%%%%%%

%load data
[vx,R] = readgeoraster(ve);
[vy,~] = readgeoraster(vn);
[vz,~] = readgeoraster(vu);
info=geotiffinfo(vu);

%%Construct coefficient matrix based on central finite difference method
%Find the solvable index of vz, as the input grid may contain a large number of empty values or zeros.
nonz=0;cc=vx;cc(:,:)=0;
for i=2:size(vx,1)-1
    for j=2:size(vx,2)-1
        if vx(i,j)~=0 && vx(i-1,j)~=0 && vx(i+1,j)~=0 && vx(i,j-1)~=0 && vx(i,j+1)~=0
            nonz=nonz+1;
            cc(i,j)=1; 
        end
    end
end
clear i j

vz(cc==0)=0;
vzc=reshape(vz',[],1);
vzc(vzc==0)=[];
z=double(-2*dx*vzc);
clear vz vzc

%Find edge value index

for i=2:size(vx,1)-1
    for j=2:size(vx,2)-1
        if cc(i,j)==0
            if cc(i-1,j)==1 || cc(i+1,j)==1 || cc(i,j-1)==1 || cc(i,j+1)==1
                cc(i,j)=2;
            end
        end
    end
end
for i=2:size(vx,1)-1
    for j=2:size(vx,2)-1
        if cc(i,j)==1
            if cc(i-1,j)==0
                cc(i-1,j)=2;
            elseif cc(i+1,j)==0
                cc(i+1,j)=2;
            elseif cc(i,j-1)==0
                cc(i,j-1)=2;
            elseif cc(i,j+1)==0
                cc(i,j+1)=2;
            end
        end
    end
end
nonzb=length(find(cc==2));
clear i j

%Obtain all solvable thickness indices
t=1;
hij=ones(nonz+nonzb,2);
for i=1:size(vx,1)
    for j=1:size(vx,2)
        if cc(i,j)~=0
            hij(t,1)=i;
            hij(t,2)=j;
            t=t+1;
        end
    end
end
clear i j t

%Build G-matrix
k=1;
G=zeros(nonz,nonz+nonzb);

for i=1:size(vx,1)
    for j=1:size(vx,2)
        if cc(i,j)==1
            for t=1:size(hij,1)
                if hij(t,1)==i && hij(t,2)==j
                    ind=t;
                elseif hij(t,1)==i-1 && hij(t,2)==j
                    inds=t;
                elseif hij(t,1)==i+1 && hij(t,2)==j
                    indx=t;
                elseif hij(t,1)==i && hij(t,2)==j-1
                    indz=t;
                elseif hij(t,1)==i && hij(t,2)==j+1
                    indy=t;
                end
            end
            clear t
            G(k,ind)=vx(i+1,j)-vx(i-1,j)+vy(i,j+1)-vy(i,j-1);%(vxi+1j-vxi-1j+vyij+1-vyij-1)hij own
            G(k,inds)=-vx(i,j);%vxij,hi-1j up
            G(k,indx)=vx(i,j);%vxij,hi+1j down
            G(k,indz)=-vy(i,j);%vyij,hij-1 left
            G(k,indy)=vy(i,j);%vyij,hij+1 right
            k=k+1;
            clear ind inds indx indz indy
        end
    end
end
clear i j k 
G=G*f;

%Building a regularization matrix
l=zeros(length(hij),length(hij));
for i=1:length(hij)
    k=0;
    for t=1:size(hij,1)
        if hij(t,1)==hij(i,1) && hij(t,2)==hij(i,2)%own
            ind=t;
        elseif hij(t,1)==hij(i,1)-1 && hij(t,2)==hij(i,2)%up
            l(i,t)=1;
            k=k+1;
        elseif hij(t,1)==hij(i,1)+1 && hij(t,2)==hij(i,2)%down
            l(i,t)=1;
            k=k+1;
        elseif hij(t,1)==hij(i,1) && hij(t,2)==hij(i,2)-1%left
            l(i,t)=1;
            k=k+1;
        elseif hij(t,1)==hij(i,1) && hij(t,2)==hij(i,2)+1%right
            l(i,t)=1;
            k=k+1;
        end
    end
    l(i,ind)=-k;
    clear k ind
end
clear i t

%Solve using the least squares norm and add positive constraints.
cvx_begin
    variable h(size(hij,1),1);

    minimize(square_pos(norms(G*h-z))+square_pos(norms(lamda*l*h)));

    subject to
      0<= h <=up;
cvx_end

%Display and save results
hh=double(vx);
for i=1:size(hij,1)
        hh(hij(i,1),hij(i,2))=h(i);
end
clear i
imshow(hh,[0 120])
colormap summer
geotiffwrite('h.tif',hh,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
clear hh h
end