classdef myPreprocess_sliding
    properties (SetAccess=public)
        files
        geo
        m
        sy
        index
        date
        nan_mask
        sup
        angle
        strsup
    end
    methods
        function obj=myPreprocess_sliding(dirc,ar,af,br,bf)
            obj.angle.a.range=ar;
            obj.angle.a.azimuth=af;
            obj.angle.d.range=br;
            obj.angle.d.azimuth=bf;

            if isunix==1
                obj.sy='/';
            else
                obj.sy='\';
            end
            obj=creat_folder(obj,dirc);
            %Define file path
            imds(1).path=imageDatastore(obj.files(1).folder,"FileExtensions",".tif","IncludeSubfolders",true);
            imds(2).path=imageDatastore(obj.files(2).folder,"FileExtensions",".tif","IncludeSubfolders",true);
            %Check data
            ishasdata(obj,imds(1).path);
            ishasdata(obj,imds(2).path);

            disp("Continue");
            disp('Start reading data:...');

            %load data
            tstart=tic;
            [~,obj.geo.R]=readgeoraster(imds(1).path.Files{1});
            obj.geo.info=geotiffinfo(imds(1).path.Files{1});
            ascending_los=load_data(obj,imds(1).path,dirc);
            save('myFile.mat',"ascending_los",'-v7.3');
            clear ascending_los
            obj.m = matfile('myFile.mat','Writable',true);

            descending_los=load_data(obj,imds(2).path,dirc);
            obj.m.descending_los=descending_los;
            clear descending_los

            tnow=toc(tstart);
            express=['Data read completed !!!','Completion time:',num2str(tnow),'s'];
            disp(express);
            clear tstart tnow

            %mask index
            tstart=tic;

            indexa=extract_mask_zeros(obj,obj.m.ascending_los,imds(1).path,dirc);
            indexd=extract_mask_zeros(obj,obj.m.descending_los,imds(2).path,dirc);
            obj.index=unique(cat(1,indexa,indexd));

            %Storage mask tif
            temp=obj.m.descending_los(:,:,1);
            temp(obj.index)=nan;
            obj.nan_mask=temp;
            temp(~isnan(temp))=1;
            temp(isnan(temp))=0;
            geotiffwrite('nan_mask.tif',temp,obj.geo.R,'GeoKeyDirectoryTag',obj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            clear temp
            tnow=toc(tstart);
            express=['Extract mask completed !!!','Completion time:',num2str(tnow),'s'];
            disp(express);

            %load date
            disp('Calculate the time:...');
            d1=[dirc,obj.sy,'date_ascending.xls'];d2=[dirc,obj.sy,'date_descending.xls'];
            obj.date.ascending=datetime(readmatrix(d1),"ConvertFrom",'yyyyMMdd');
            obj.date.descending=datetime(readmatrix(d2),"ConvertFrom",'yyyyMMdd');
            obj.date.time=unique([obj.date.ascending(:,1);obj.date.ascending(:,2);obj.date.descending(:,1);obj.date.descending(:,2)]);
            obj.date.delta_t=days(diff(obj.date.time));
            clear d1 d2

            %Import supporting files
            obj.strsup=[dirc,obj.sy,'supporting_documents'];
            if isfolder(obj.strsup)==1
                disp('Import supporting files:...');
                sup1=[dirc,obj.sy,'supporting_documents'];
                imds(3).path=imageDatastore(sup1,"FileExtensions",".tif","IncludeSubfolders",true);
                obj=load_sup_data(obj,imds(3).path);
            end
            disp('Preprocessing completed !!!!!!!!!!!!!!!');

        end

        function obj=creat_folder(obj,dirc)
            files1(1).name='ascending_los';
            files1(1).folder=[dirc,obj.sy,'ascending_los'];
            files1(2).name='descending_los';
            files1(2).folder=[dirc,obj.sy,'descending_los'];
            obj.files=files1;
        end

        function ishasdata(obj,data)
            if hasdata(data)==0
                disp("error: No data in folder");
            end
        end

        function data=load_data(obj,imds,dirc)
            dirc1=[dirc,obj.sy];
            h=waitbar(0,'Start reading data.....');
            [eg,~]=readgeoraster(imds.Files{1});
            data=single(zeros(size(eg,1),size(eg,2),length(imds.Files)));
            k=length(imds.Files);
            for i=1:k
                [tt,~]=readgeoraster(imds.Files{i});
                tt(tt>10000)=0;tt(tt<-10000)=0;
                data(:,:,i)=tt;
                str=['Exporting  ',extractAfter(imds.Folders{1,1},dirc1),',  current progress is ',num2str(i),'/',num2str(k)];
                waitbar(i/k,h,str);
                clear tt
            end
            disp(['Exporting  ',extractAfter(imds.Folders{1,1},dirc1)]);
            close(h)
            clear h i k str eg
        end

        function indexout=extract_mask_zeros(obj,a,imds,dirc)
            h=waitbar(0,'Start extract masking data.....');
            index1=find(a(:,:,1)==0);
            dirc1=[dirc,obj.sy];
            for i=1:size(a,3)
                index2=find(a(:,:,i)==0);
                index3=cat(1,index1,index2);
                index1=unique(index3);
                str=['Extract masking ',extractAfter(imds.Folders{1,1},dirc1),',  current progress is ',num2str(i),'/',num2str(size(a,3))];
                waitbar(i/size(a,3),h,str);
                clear index2 index3
            end
            indexout=index1;
            close(h)
            disp(['Exporting  ',extractAfter(imds.Folders{1,1},dirc1)]);
            clear i index1 a
        end

        function obj=load_sup_data(obj,imds)
            str1=['supporting_documents',obj.sy];
            str2='.tif';
            for i=1:length(imds.Files)
                str_name=extractAfter(imds.Files{i},str1);
                str_name=extractBefore(str_name,str2);
                if strcmp(str_name,'dem_slope')==1
                    obj.sup.dem_slope=readgeoraster(imds.Files{i});
                elseif strcmp(str_name,'dem_aspect')==1
                    obj.sup.dem_aspect=readgeoraster(imds.Files{i});
                end
                clear str_name
            end
        end
    end

end