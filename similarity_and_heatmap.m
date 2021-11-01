clc;
clear;

[file,path] = uigetfile('C:\Users','MultiSelect','on');

experi_num=length(file);



load('data_compiled_optrode_test91.mat');
load('CS_responsive28.mat');

do_only_CS=1;


if do_only_CS == 1
    
    data_compiled=data_compiled_test(1,CS_responsive(:));
    
end
    
%%
    
before_cue=000;
after_cue=15000;
total_length=20000+before_cue+after_cue;

invivo_bin_size=100;

cueNdelay=bin_cueNdelay(data_compiled, before_cue,after_cue,invivo_bin_size)';
cueNdelay_mean_10sec=sum(cueNdelay,1)./10/22;

for j=1:experi_num
    
    filename=[path file{1,j}];  
    r=load(filename);
    r=r.r;
    r=r(:,1:9000);
    
    
    
    r_mean=mean(r,1);
    r_sem=std(r)/sqrt(length(r));
    
    time=(1:length(r))*0.005;
    
    bin_time=invivo_bin_size/1000;
    
    bin_size=bin_time/0.005;
    
    bin_r=reshape(sum(reshape(r,1000,bin_size,[]),2),1000,[]);
    bin_r_mean=mean(bin_r,1);
    bin_r_sem=std(bin_r)/sqrt(length(bin_r));
    
    
    %heatmap
    figure()
    bin_r_sum=sum(bin_r,2);
    [B,Index_sort]=sort(bin_r_sum);
    
    bin_r=bin_r(Index_sort,:);
    
    timepoint=bin_time:bin_time:65;
    
    N=ones(size(bin_r_mean,2),1)*size(r,1);
    
    xlabel=NaN(size(bin_r,2),1);
    xlabel_num=ceil(size(bin_r,2)/10)-1;
    xlabel(1:100:(xlabel_num+1)*10)=0:10:xlabel_num;
    
    ylabel=NaN(size(bin_r,1),1);
    ylabel(1:100:901)=0:100:900;
    
    
    [pks,locs] = max(bin_r,[],2);
    temp = [locs bin_r];
    temp = temp(:,2:end);
    h=heatmap(temp);
    h.XDisplayLabels=xlabel;
    h.YDisplayLabels=ylabel;
    h.MissingDataLabel='';
    h.GridVisible='off';
    h.Colormap=parula;
    h.ColorbarVisible='off';
    
    %similarity score
    corr_r=corrcoef(bin_r(:,10/bin_time+1:size(bin_r,2)));
    corr_invivo=corrcoef(cueNdelay);
    analysis_length=length(corr_invivo);
    
    for i=1:analysis_length-1
        for k=1:analysis_length-i
            corr_temp(k)=corr_invivo(k+i,k);
        end
        corr_invivo_mean(i)=nanmean(corr_temp);
        clear corr_temp;
    end
    
    for i=1:analysis_length-1
        for k=1:analysis_length-i
            corr_temp(k)=corr_r(k+i,k);
        end
        corr_r_mean(i)=nanmean(corr_temp);
        clear corr_temp;
    end
    
    corr_sum=nansum(abs(corr_invivo_mean(1:analysis_length-1)...
        -corr_r_mean(1:analysis_length-1)));
    similarity(j)=1-corr_sum/(analysis_length-1);
    %similarity(j)=1-corr_sum/abs(nansum(corr_invivo_mean(1:analysis_length-1)));
    
    period_activity(j,:)=[mean(bin_r_mean(10/bin_time+1:30/bin_time))...
        mean(bin_r_mean(30/bin_time+1:35/bin_time)) mean(bin_r_mean(35/bin_time+1:45/bin_time)) ];
    norm_period_activity=normalize(period_activity,2,'norm');
%     
%     invivo_mean=mean(cueNdelay(CS_responsive_neurons,:),1);
%     period_activity_invivo=[mean(invivo_mean(1:20)) mean(invivo_mean(21:25))...
%         mean(invivo_mean(26:50)) ];
%     norm_period_activity_invivo=normalize(period_activity_invivo,2,'norm');
%     
%     period_activity_noinhi(j,:)=[mean(bin_r_mean(10/bin_time+1:30/bin_time))...
%         mean(bin_r_mean(35/bin_time+1:65/bin_time))];
%     norm_period_activity_noinhi=normalize(period_activity_noinhi,2,'norm');
    
    
end

%% period activity normlaized

for i=1:size(period_activity,1)
    cue_norm_activity(i,:)=period_activity(i,:)./period_activity(i,1);
end

