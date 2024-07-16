classdef myProcess_sliding
    properties (SetAccess=public)
        aspect
        m2
        lamda
        normv
        flag
    end

    methods
        function obj=myProcess_sliding(preobj,lamda)
            %Compute
            disp('Bulid storage output:...')
            obj.lamda=lamda;
            
            obj.aspect.length=size(preobj.nan_mask,1);
            obj.aspect.width=size(preobj.nan_mask,2);
            obj.aspect.phase=size(preobj.date.delta_t,1);
            obj.aspect.a=size(preobj.date.ascending,1);
            obj.aspect.d=size(preobj.date.descending,1);

            obj.aspect.tikhonov=(length(preobj.date.delta_t)-1)*2;
            obj.aspect.pairs=obj.aspect.a+obj.aspect.d+obj.aspect.tikhonov;

            vr=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vi=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));

            a.a=preobj.m.ascending_los;a.d=preobj.m.descending_los;

            tstart=tic;
            disp('Starting calculation may take a long time, please be patient and wait:...');
            count=1;all_count=obj.aspect.length;lable=0;str='Completed  0%';
            h=waitbar(0,'Please be patient and wait.....');
            waitbar(0,h,str);
            norm_x=0;countt=0;norm_axy=0;

            for i=1:obj.aspect.length
                for j=1:obj.aspect.width
                    if ~isnan(preobj.nan_mask(i,j)) %Determine whether it is a null value
                        b=make_b(obj,a,i,j); %b is deformation data pairs
                        A=make_A(obj,preobj,i,j); %A is time coefficient matrix
                        temp=lsqminnorm(A,b);
                        %temp=pinv(A)*b; %SVD

                        for k=1:obj.aspect.phase
                            vr(i,j,k)=temp((k-1)*2+1);
                            vi(i,j,k)=temp((k-1)*2+2);
                        end
                        norm_t=squeeze(cat(3,vr(i,j,:),vi(i,j,:)));
                        if  norm(norm_t)<10000
                            norm_x=norm_x+norm(norm_t);
                            norm_axy=norm_axy+norm(A*temp);
                            countt=countt+1;
                        end
                        clear temp A b k norm_t
                    end
                end
                count_a=100*count/all_count;
                count_b=floor(count_a);
                count_c=count_a-count_b;
                if lable>count_c
                    str=['Completed  ',num2str(count_b),'%'];
                    disp(str);
                    waitbar(i/all_count,h,str);
                end
                lable=count_c;
                count=count+1;
            end
            
            obj.normv.x=norm_x/countt;
            obj.normv.axy=norm_axy/countt;
            disp(['Computed ||x|| and ||Ax-Y|| norms: ',num2str(obj.normv.x),'    ',num2str(obj.normv.axy)]);
            close(h)
            clear h lable count_c count_a all_count count_b str count i j q a

            disp('Storing result as mat........myFile_out-----velocity');
            for k=1:obj.aspect.phase
                vr(:,:,k)=del_std(obj,vr(:,:,k));
                vi(:,:,k)=del_std(obj,vi(:,:,k));
            end
            save('myFile_out.mat',"vr",'-v7.3');
            obj.m2 = matfile('myFile_out.mat','Writable',true);
            obj.m2.vi=vi;

            tnow=toc(tstart);
            disp(['Completion time:',num2str(tnow/3600),'h']);
            clear tstart tnow

            %Create cumulative time series
            disp('Start building cumulative time series output.....');
            tstart=tic;
            vr_cm=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vi_cm=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            for k=1:obj.aspect.phase
                if k==1
                    vr_cm(:,:,k)=vr(:,:,k).*preobj.date.delta_t(k);
                    vr_cm(:,:,k)=del_std(obj,vr_cm(:,:,k));
                    vi_cm(:,:,k)=vi(:,:,k).*preobj.date.delta_t(k);
                    vi_cm(:,:,k)=del_std(obj,vi_cm(:,:,k));
                else
                    vr_cm(:,:,k)=vr(:,:,k).*preobj.date.delta_t(k)+vr_cm(:,:,k-1);
                    vr_cm(:,:,k)=del_std(obj,vr_cm(:,:,k));
                    vi_cm(:,:,k)=vi(:,:,k).*preobj.date.delta_t(k)+vi_cm(:,:,k-1);
                    vi_cm(:,:,k)=del_std(obj,vi_cm(:,:,k));
                end
            end
            obj.m2.vr_cm=vr_cm;
            obj.m2.vi_cm=vi_cm;

            vr_cm_v=(vr_cm(:,:,end)./days(preobj.date.time(end)-preobj.date.time(1))).*365;
            vr_cm_v=del_std(obj,vr_cm_v);
            vi_cm_v=(vi_cm(:,:,end)./days(preobj.date.time(end)-preobj.date.time(1))).*365;
            vi_cm_v=del_std(obj,vi_cm_v);
            obj.m2.vr_cm_v=vr_cm_v;
            obj.m2.vi_cm_v=vi_cm_v;
            tnow=toc(tstart);
            disp(['completion time: ',num2str(tnow),'s']);
            clear tstart tnow k

            %Export as tif
            disp('Export as tif write to disk.....');
            tstart=tic;
            if isunix==1
                mkdir results_velocity/vr
                mkdir results_velocity/vi
                mkdir results_cumulative/vr
                mkdir results_cumulative/vi
            else
                mkdir results_velocity\vr
                mkdir results_velocity\vi
                mkdir results_cumulative\vr
                mkdir results_cumulative\vi
            end
            for k=1:obj.aspect.phase
                str_vr_v=['results_velocity',preobj.sy,'vr',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vr_v,vr(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vi_v=['results_velocity',preobj.sy,'vi',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vi_v,vi(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vr_cm=['results_cumulative',preobj.sy,'vr',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vr_cm,vr_cm(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vi_cm=['results_cumulative',preobj.sy,'vi',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vi_cm,vi_cm(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

            end
            geotiffwrite('vr_ave_year.tif',vr_cm_v,preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite('vi_ave_year.tif',vi_cm_v,preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            tnow=toc(tstart);
            disp(['Completion time:',num2str(tnow),'s']);
            clear vr vi vr_cm_v vi_cm_v vr_cm vi_cm str tstart tnow str_vr_v str_vi_v str_vr_cm str_vi_cm
        end

        function A=make_A(obj,preobj,i,j)
            %Cteate time coefficient matrix
            A=zeros(obj.aspect.pairs,length(preobj.date.delta_t)*2);
            time=preobj.date.time;
            delta_t=preobj.date.delta_t;

            alpha=preobj.sup.dem_slope(i,j);
            belta=preobj.sup.dem_aspect(i,j);
            thetaa=preobj.angle.a.range;
            thetad=preobj.angle.d.range;
            phia=preobj.angle.a.azimuth;
            phid=preobj.angle.d.azimuth;

            p1a=sind(thetaa)*sind(phia)*cosd(alpha)*cosd(belta)-sind(thetaa)*cosd(phia)*cosd(alpha)*sind(belta)-cosd(thetaa)*sind(alpha);
            p3a=cosd(thetaa)*cosd(alpha)+sind(thetaa)*sind(phia)*sind(alpha)*cosd(belta)-sind(thetaa)*cosd(phia)*sind(alpha)*sind(belta);

            p1d=sind(thetad)*sind(phid)*cosd(alpha)*cosd(belta)-sind(thetad)*cosd(phid)*cosd(alpha)*sind(belta)-cosd(thetad)*sind(alpha);
            p3d=cosd(thetad)*cosd(alpha)+sind(thetad)*sind(phid)*sind(alpha)*cosd(belta)-sind(thetad)*cosd(phid)*sind(alpha)*sind(belta);

            %1:a
            for k=1:obj.aspect.a
                index1=find(time==preobj.date.ascending(k,1));
                index2=find(time==preobj.date.ascending(k,2))-1;
                A(k,2*(index1-1)+1)=delta_t(index1)*p1a;
                A(k,2*(index1-1)+2)=delta_t(index1)*p3a;

                while (index2-index1)~=0
                    A(k,2*(index2-1)+1)=delta_t(index1)*p1a;
                    A(k,2*(index2-1)+2)=delta_t(index1)*p3a;
                    index2=index2-1;
                end
                clear index1 index2
            end

            %a:d
            for k=obj.aspect.a+1:obj.aspect.d+obj.aspect.a
                index1=find(time==preobj.date.descending(k-obj.aspect.a,1));
                index2=find(time==preobj.date.descending(k-obj.aspect.a,2))-1;
                A(k,2*(index1-1)+1)=delta_t(index1)*p1d;
                A(k,2*(index1-1)+2)=delta_t(index1)*p3d;

                while (index2-index1)~=0
                    A(k,2*(index2-1)+1)=delta_t(index1)*p1d;
                    A(k,2*(index2-1)+2)=delta_t(index1)*p3d;
                    index2=index2-1;
                end

                clear index1 index2
            end

            %d:Tikhonov
            for k=obj.aspect.a+1+obj.aspect.d:obj.aspect.pairs
                q=k-obj.aspect.a-obj.aspect.d;
                A(k,q)=-obj.lamda;
                A(k,q+2)=obj.lamda;
                clear q
            end

            clear k time delta_t phia phid thetaa thetad alpha belta
        end

        function b=make_b(obj,a,i,j)
            b=zeros(obj.aspect.pairs,1);
            for k=1:obj.aspect.a
                b(k)=a.a(i,j,k);
            end
            for k=1:obj.aspect.d
                b(k+obj.aspect.a)=a.d(i,j,k);
            end
            clear k
        end

        function e=del_std(obj,a)
            b=std(a,1,'all','omitnan' );
            c=mean(a,'all','omitnan' );
            d=[c-10000*b c+10000*b];
            a(a>d(2))=nan;
            a(a<d(1))=nan;
            e=a;
            clear b c d a
        end


    end
end