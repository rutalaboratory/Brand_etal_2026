% code accompanying Brand et al. (2026). Ecological and social cues orchestrate concerted courtship in Drosophila. Nature 
%
% loads post processed flytracker output .mat file 
% run through Brand_et_al_26_flytracker_postProcessing.m
%
% script loads all replicate experiments of a given experimental condition in a folder 
% (as indicated by the name) and groups them
% correct data loading requires the following naming scheme of files
% Imaging data: species_individuals_age_housing_chemicalEnvironment.mat


%%
clear;
fps = 42;
assay_length = fps*60*40; %40 minute assay length cutoff

start_path = fullfile('/Users/pbrand/science/behavior/free_behavior/assay_tracks/');
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end

filePattern = sprintf('%s/*.mat', topLevelFolder);
baseFileNames = dir(filePattern);
numberOfMatFiles = length(baseFileNames);
    
spec = {'SPECIES'};
inds = {'INDIVIDUALS'};
age ={'AGE'};
housing = {'HOUSING'};
chem = {'CHEMeNIRONMENT1','CHEMeNVIRONMENT2'};

groups = strcat({'group'},spec,{'_'},inds, {'_'}, age, {'_'}, housing, {'_'},chem);
numGroups=length(groups);

patterns = strcat(spec,{'.*'},inds, {'.*'},age, {'.*'}, housing, {'.*'},chem);

%prepare data structures to store parameters saved in .mat files in
att = struct();
ii=struct();
ai=struct();
ai_maps=struct();
focal_maps=struct();
sex = struct();
cop =struct();
ai_blocks=struct();
euc_dist=struct();

for pat = 1:length(patterns)
    thisPattern=patterns{pat};
    %fprintf('checking for %s\n',thisPattern)
    
    hits={};
    %get matches to pattern - these are assays within the same group
    cop.(groups{pat}).cpm = [];
    cop.(groups{pat}).cpf = [];
    for f = 1 : numberOfMatFiles
        if regexp(baseFileNames(f).name, thisPattern)
            %fprintf('\tfound %s\n',baseFileNames(f).name)
            hits{length(hits)+1} = baseFileNames(f).name;
        end
    end
        
    numHits=length(hits);
        
    for hit = 1:numHits
        fullFile = fullfile(topLevelFolder, hits{hit});
        fprintf('loading file %s\n',hits{hit})
        load(fullFile);
        mobile.(groups{pat})(hit).name=hits{hit};
        
        hitref=sprintf('map%d',hit);

        %save and calculate parameters from .mat file data

        att.(groups{pat}).(hitref) = struct('mm',[],'mf',[],'fm',[],'ff',[]);

        sex.(groups{pat}).(hitref)=results.sex;
        
        wing_map.(groups{pat}).(hitref)=results.wing_map;
        
        approach_map.(groups{pat}).(hitref)=results.approach_map;
        
        numflies=length(results.flies);
        
        for individual1=1:length(results.flies)-1
           
            for individual2=individual1+1:length(results.flies)
                       

                euc_dist.(groups{pat}).(hitref).eucCell{individual1,individual2}=sqrt( (results.x(:,individual2)-results.x(:,individual1)).^2 + (results.y(:,individual2)-results.y(:,individual1)).^2  );
                
            end
            
        end
        
        
        
        [ai_maps.(groups{pat}).(hitref),focal_maps.(groups{pat}).(hitref)]=mapLoader(results.interaction_map,results.approach_map,assay_length,numflies);
        
        [ii.(groups{pat}).(hitref).total,ii.(groups{pat}).(hitref).m,ii.(groups{pat}).(hitref).f]=interactionCounter(results.mating_blocks,focal_maps.(groups{pat}).(hitref),assay_length,numflies,results.sex,fps);
      
        [ai.(groups{pat}).(hitref).total,ai.(groups{pat}).(hitref).m,ai.(groups{pat}).(hitref).f]=attractionIndex(ai_maps.(groups{pat}).(hitref),assay_length,numflies,results.sex,fps);
        
        [att.(groups{pat}).(hitref).mm,att.(groups{pat}).(hitref).mf, att.(groups{pat}).(hitref).fm,att.(groups{pat}).(hitref).ff,att.(groups{pat}).(hitref).fa] = actionThroughTime(focal_maps.(groups{pat}).(hitref),fps,assay_length,numflies,results.sex);
        
        cop.(groups{pat}).(hitref).ind=cellfun('size',results.mating_blocks,1);
        
        cop.(groups{pat}).(hitref).blocks=results.mating_blocks;
        
        cop.(groups{pat}).(hitref).cpm = sum(cop.(groups{pat}).(hitref).ind(sex.(groups{pat}).(hitref)==1))/sum(sex.(groups{pat}).(hitref)==1);
        cop.(groups{pat}).(hitref).cpf = sum(cop.(groups{pat}).(hitref).ind(sex.(groups{pat}).(hitref)==0))/sum(sex.(groups{pat}).(hitref)==0);
        cop.(groups{pat}).cpm = [cop.(groups{pat}).cpm cop.(groups{pat}).(hitref).cpm];
        cop.(groups{pat}).cpf = [cop.(groups{pat}).cpf cop.(groups{pat}).(hitref).cpf];

        %x and y
        positions.(groups{pat})(hit).x=results.x(1:assay_length,:);
        positions.(groups{pat})(hit).y=results.y(1:assay_length,:);

        iblocks.(groups{pat}).(hitref)=results.interaction_blocks;
        velocity.(groups{pat}).(hitref)=results.vel;
        ablocks.(groups{pat}).(hitref)=results.approach_blocks;
        
        ai_blocks.(groups{pat}).(hitref)=blocker(ai_maps.(groups{pat}).(hitref),10);
        wing_blocks.(groups{pat}).(hitref)=blocker(wing_map.(groups{pat}).(hitref),1);
        
        ai_block_flip.(groups{pat}).(hitref) = cell(numflies,1);
        %flip blocks from 'received' to 'administered' (from target focus to executing focus)
        
        for tar = 1:length(ai_blocks.(groups{pat}).(hitref))
           
            for bl = 1:size(ai_blocks.(groups{pat}).(hitref){tar},1)
                
                check= size(ai_block_flip.(groups{pat}).(hitref),1);
                
                if check<1
                    ai_block_flip.(groups{pat}).(hitref){ ai_blocks.(groups{pat}).(hitref){tar}(bl,3),1 } =  ai_blocks.(groups{pat}).(hitref){tar}(bl,:);
                    ai_block_flip.(groups{pat}).(hitref){ ai_blocks.(groups{pat}).(hitref){tar}(bl,3),1 }(end,3)=tar;
                else
                    tmper=ai_blocks.(groups{pat}).(hitref){tar}(bl,:);
                    tmper(3) = tar;
                    ai_block_flip.(groups{pat}).(hitref){ ai_blocks.(groups{pat}).(hitref){tar}(bl,3),1 } = [ai_block_flip.(groups{pat}).(hitref){ ai_blocks.(groups{pat}).(hitref){tar}(bl,3),1 }; tmper];

                end
            end
            
        end
        
    end
end

%% plot interaction and mating stats
figure(61)
assay_vals_by_group= struct();
cmap=cool();
subplot(2,2,[1 2])
ii_males_groups={};
ii_females_groups={};
modl=1;
for i = modl:numGroups
    
   %totals first 
   numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));

   colors=cmap(round(linspace(1,length(cmap),numAssays)),:);
   
   vals_m=[];
   vals_f=[];
   cop_m=[];
   cop_f=[];
   iscop_m=[];
   iscop_f=[];
   ci(i).ci=[];
   for j = 1:numAssays
       
       assay_sex=sex.(groups{i}).(assays{j});
       assay_total=ii.(groups{i}).(assays{j}).total;
       males=find(assay_sex==1);
       females=find(assay_sex==0);
       
       for sm = 1:length(males)
          
           if length(cop.(groups{i}).(assays{j}).blocks{males(sm)})
               
               cop_m=[cop_m cop.(groups{i}).(assays{j}).blocks{males(sm)}(1)];
               iscop_m=[iscop_m 2];
           else
               cop_m=[cop_m nan()];
               iscop_m=[iscop_m 1];
           end
       end
       
       for sf = 1:length(females)
          
           if length(cop.(groups{i}).(assays{j}).blocks{females(sf)})
               
               cop_f=[cop_f cop.(groups{i}).(assays{j}).blocks{females(sf)}(1)];
               iscop_f=[iscop_f 2];
           else
               cop_f=[cop_f nan()];
               iscop_f=[iscop_f 1];
           end
       end
       
       
       
%        assay_vals_by_group.m_means{i,1}(j)=nanmean(assay_total(5,males));
%        assay_vals_by_group.m_vals{i,j}=assay_total(5,males);
%        assay_vals_by_group.f_means{i,1}(j)=nanmean(assay_total(5,females));
%        assay_vals_by_group.f_vals{i,j}=assay_total(5,females);
%        assay_vals_by_group.m_stds{i,1}(j)=nanstd(assay_total(5,males));
%        assay_vals_by_group.f_stds{i,1}(j)=nanstd(assay_total(5,females));


 ci(i).ci=[ci(i).ci; assay_total(5,males)];
       vals_m=[vals_m mean(assay_total(5,males),'omitnan')];
       vals_f=[vals_f mean(assay_total(5,females),'omitnan')];
       
       
       scatter(i*ones(length(males),1)',assay_total(5,males),[],colors(j,:),'filled','jitter', 'on', 'jitterAmount', 0.1)
       xlim([0.5 (numGroups*2)+0.5])
       ylim([0 1])
       hold on
       scatter((i+numGroups)*ones(length(females),1)',assay_total(5,females),[],colors(j,:),'filled','jitter', 'on', 'jitterAmount', 0.1)

   end
   
   ii_males_groups{i,1}=vals_m;
   ii_females_groups{i,1}=vals_f;

    plot([i-0.2 i+0.2],[nanmean(vals_m) nanmean(vals_m)],'Linewidth',3,'Color','k')
    plot([i+numGroups-0.2 i+numGroups+0.2],[nanmean(vals_f) nanmean(vals_f)],'Linewidth',3,'Color','k')
    xlabel('Interaction index')
    xticks([1:numGroups*2])
    
    fprintf('%s: mean ii males %s\n', groups{i}, nanmean(vals_m));
    fprintf('%s: mean ii females %s\n', groups{i}, nanmean(vals_f));
    
    name_specs={};
    if size(chem,2)>1
        name_specs=chem;
    else 
        if size(housing,2)>1
           name_specs=housing;
        else
            
            if size(inds,2)>1
                name_specs=inds;
                
            else
                    
                  if size(spec,2)>1
                        name_specs = spec;
                  else
                        name_specs=age;
                  end
            end
        end
    end
    
    xlabs=strcat ( 'males','-',name_specs);
    tmp = strcat ( 'females','-',name_specs);
    xlabs=[xlabs tmp];
    xticklabels(xlabs);
    ylabel('Interaction index')
    
    cop.(groups{i}).cop_f=cop_f;
    cop.(groups{i}).court_f=vals_f;
    cop.(groups{i}).iscop_f=iscop_f;
    
    cop.(groups{i}).cop_m=cop_m;
    cop.(groups{i}).court_m=vals_m;
    cop.(groups{i}).iscop_m=iscop_m;
end

subplot(2,2,3)
nums=zeros(3,numGroups);
for i = modl:numGroups
   numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));

   colors=cmap(round(linspace(1,length(cmap),numAssays)),:);
  x= i*ones(length(cop.(groups{i}).cpm),1);
  scatter(x,cop.(groups{i}).cpm,'filled','jitter', 'on', 'jitterAmount', 0.1)
  xlim([0.75 numGroups+0.25])
  ylim([0 1])
  mean_cop = nanmean(cop.(groups{i}).cpm);
  hold on

  plot([i-0.15 i+0.15], [mean_cop mean_cop],'LineWidth',2,'Color','k')
  ylabel('Copulations per male')
  xticks([1:numGroups]);
  
  
   
       name_specs={};
    if size(chem,2)>1
        name_specs=chem;
    else 
        if size(housing,2)>1
           name_specs=housing;
        else
            
            if size(inds,2)>1
                name_specs=inds;
                
            else
                    
                  if size(spec,2)>1
                        name_specs = spec;
                  else
                        name_specs=age;
                  end
            end
        end
    end
    
  xticklabels(name_specs);
  
  
  for j = 1:numAssays
  
   male_cops= cop.(groups{i}).(assays{j}).ind(sex.(groups{i}).(assays{j})==1);
   
  
   for k = 1:length(male_cops)
      
       nums(male_cops(k)+1,i)=nums(male_cops(k)+1,i)+1;
       
   end
  
  
  end
  
end

subplot(2,2,4)

for i = 1: size(nums,2)
   
    nums(:,i)=nums(:,i)/sum(nums(:,i));
    
end

for i =modl:numGroups
    if i==1
        nums=[nums;  sum(cop.(groups{i}).iscop_f==2)/length(cop.(groups{i}).iscop_f) zeros(1,numGroups-1)];
    else
        nums(size(nums,1),i)=sum(cop.(groups{i}).iscop_f==2)/length(cop.(groups{i}).iscop_f);
    end
end

bar(nums)
ylabel({'Fraction of males'})
xticks([1 2 3 4])
xticklabels({'none', 'one', 'two','female rate'});
xlabel('Matings')
ylim([0 1])


%%%%plot activity through time
numGroups=length(groups);
cols={'#00CCFF','#FF9933'};
colsRGB=[0 0.6 1; 1 0.6 0];

cases={'fm','ff','fa','mm','mf'};
for int = 1:5
figure(int)
    
    bins = length(att.(groups{1}).map1.(cases{1}));
    x=1:bins;

    means = zeros(numGroups,bins);
    sds = means;
    for i = modl: numGroups

        numassays=length(fieldnames(att.(groups{i})));
        assays = fieldnames(att.(groups{i}));


        merged=zeros(numassays,bins);

        merged(1,:)=att.(groups{i}).(assays{1}).(cases{int});
        
        subplot(numGroups+1,1,i)

        if mod(i,2)>0
            colType=1;
        else
            colType=2;
        end
        plot(x,att.(groups{i}).(assays{1}).(cases{int}),'Color',cols{colType})
        ylabel('Fraction interacting')
        
        sub={};
        if size(chem,2)>1
            sub=chem;
        else 
            if size(housing,2)>1
               sub=housing;
            else
                if size(inds,2)>1
                sub=inds;
                
                else
                    
                    if size(spec,2)>1
                        sub = spec;
                    else
                        sub=age;
                    end
                end
            end
        end

        
        title( strcat ( sub(i),'-',cases{int}) );
        
        ylim([-0.05 1.05])
        xlim([0 bins])
        xticks([1 60 120 180 240 300 360 420 478])
        xticklabels({'0','5','10','15','20','25','30','35','40'})

        hold on
        for j = 2:numassays
            merged(j,:)=att.(groups{i}).(assays{j}).(cases{int});

            plot(x,att.(groups{i}).(assays{j}).(cases{int}),'Color',cols{colType})

        end
            plot(x,nanmean(merged),'Color','k')
            means(i,:) = nanmean(merged);
            sds(i,:) = nanstd(merged)/sqrt(numassays);

    end
        subplot(numGroups+1,1,i+1)
        
        
        for k = 1:numGroups
            curve1 = means(k,:) + sds(k,:);
            curve2 = means(k,:) - sds(k,:);
            x2 = [x, fliplr(x)];
            inBetween = [curve1, fliplr(curve2)];

            if mod(k,2)
                colType=1;
            else
                colType=2;
            end

            h=fill(x2, inBetween, colsRGB(colType,:));
            h.FaceAlpha = 0.5;
            set(h,'EdgeColor','none');
            ylim([-0.05 1.05])
            xlim([0 bins])
            
            xticks([1 60 120 180 240 300 360 420 478])
            xticklabels({'0','5','10','15','20','25','30','35','40'})
            
%             xlim([0 60]) %use for first 5 min window
%              xticks([1 12 24 36 48 60])
%              xticklabels({'0','1','2','3','4','5'})
            
            ylabel('Fraction interacting')
            xlabel('time [min]')
        
            hold on
            plot(x,means(k,:),'Color','k')
        end
            hold off

end   

for i = modl:numGroups
   
    iscop=sum(cop.(groups{i}).cpm);
    xcumul=zeros(1,bins);
    if iscop>0
        copSort=sort(cop.(groups{i}).cop_f);
        copSort=round(copSort/210.8787); %map on bins
        totalCops=sum(~isnan(copSort));
        for k=1:totalCops
            xcumul(copSort(k):end)=xcumul(copSort(k):end)+1;
        end
        xcumul=xcumul/totalCops;
    end
    hold on
    
    if mod(i,2)>0
        colType=1;
    else
        colType=2;
    end

    plot(xcumul,'--','Color',cols{colType},'LineWidth',2)
    
end


%% plot mean velocity and velocity pre and post interaction start
%what are females doing prior to and after start of a courtship bout?
t_test_vel={};
for i = 1:numGroups
   figure(6)
   numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
  
   vel_snippet_means=[];
   mean_delta_vel=[];
   for j = 1:numAssays

       %counts only pre-copulation velocity!
       
      numFliegen=size(sex.(groups{i}).(assays{j}),2);
      %vel_mean.(groups{i}).(assays{j})=nanmean(velocity.(groups{i}).(assays{j}),1);
      vel_mean.(groups{i}).(assays{j})=nan(1,numFliegen);
      for fliege = 1:numFliegen
            vel_end=assay_length;
            if cop.(groups{i}).(assays{j}).ind(fliege)==1
                
                cbl=cop.(groups{i}).(assays{j}).blocks{fliege};
                
                vel_end=cbl(1)-1;
            end
            vel_pre_cop=velocity.(groups{i}).(assays{j})(1:vel_end,fliege);
            vel_mean.(groups{i}).(assays{j})(fliege)=nanmean(vel_pre_cop);
      end


      vel_snippets=[];
      if sum(sex.(groups{i}).(assays{j})==1)>0
        mint_blocks=iblocks.(groups{i}).(assays{j}){sex.(groups{i}).(assays{j})==1,1}; %male interaction blocks
      
        for k = 1:size(mint_blocks,1)
         
          if sex.(groups{i}).(assays{j})(mint_blocks(k,3))==0 && mint_blocks(k,4)>=fps*5 && mint_blocks(k,1)+fps*10<=assay_length%interaction w/ female for longer or equal 5 seconds?
          
              midpoint=mint_blocks(k,1);
              pre=midpoint-fps*10;
              post=midpoint+fps*10;
              target=mint_blocks(k,3);
              
              add_pre=0;
              if pre<1
                  add_pre=pre*(-1);
                 pre=1; 
                 target_vel=[NaN(add_pre+1,1); velocity.(groups{i}).(assays{j})(pre:post,target)];
              else
                  target_vel=velocity.(groups{i}).(assays{j})(pre:post,target);
              end
              
              
              
              if max(target_vel)<=35 %remove tracks containing speeds >35 mm/sec
                  vel_snippets=[vel_snippets target_vel];

                  x=1:length(target_vel);

                  velCol=[0 0.5 1];%blue

                  %plot(x,target_vel,'Color',velCol)
                  %ylim([0 35])
                  %hold on
              end
          end
        end
      
      end
      if ~isempty(vel_snippets) && size(vel_snippets,2)>2
          mean_cur=nanmean(vel_snippets');
        vel_snippet_means=[vel_snippet_means; mean_cur];
        
        delta_vel=[];
        
        for dv=fps+1:fps:fps*19+1
            sec_mean1=nanmean(mean_cur(dv-fps:dv-1));
            sec_mean2=nanmean(mean_cur(dv:dv+fps));
            delta_vel=[delta_vel sec_mean2-sec_mean1];
        
        end
        mean_delta_vel=[mean_delta_vel; delta_vel];
      end
   end
   
   subplot(2,1,1)

    col=[0 0.6 1; 1 0.6 0];
    hold on
    mean_vel_snip_mean=nanmean(vel_snippet_means,1);
    sd_vel_snip_means=nanstd(vel_snippet_means,1)/sqrt(size(vel_snippet_means,1));
   
        if mod(i,2)>0
            colType=1;
        else
            colType=2;
        end

    if ~isnan(mean_vel_snip_mean)
        curve1 = mean_vel_snip_mean + sd_vel_snip_means;
                curve2 = mean_vel_snip_mean - sd_vel_snip_means;
                x2 = [x, fliplr(x)];
                inBetween = [curve1, fliplr(curve2)];
                h=fill(x2, inBetween, col(colType,:));
                h.FaceAlpha = 0.5;
                set(h,'EdgeColor','none');

        plot(x,mean_vel_snip_mean,'k','LineWidth', 1.5)

        xline(fps*10+1,'r--','LineWidth',3)
        xticks([1    84   168   252   336   421   504   588   672   756   840])
        xticklabels({'-10','-8','-6','-4','-2','0','+2','+4','+6','+8','+10'})
        xlabel('time from courtship bout start [sec]')
        ylabel('mean female velocity [mm/s]')
        xlim([1 840])
    end
       subplot(2,1,2)
        hold on
       x=1:size(mean_delta_vel,2);
       mean_dv=nanmean(mean_delta_vel,1);
       sds_dv=nanstd(mean_delta_vel,1)/sqrt(size(mean_delta_vel,2));
       
       if ~isnan(mean_dv)
            curve1 = mean_dv + sds_dv;
                    curve2 = mean_dv - sds_dv;
                    x2 = [x, fliplr(x)];
                    inBetween = [curve1, fliplr(curve2)];
                    h=fill(x2, inBetween, col(colType,:));
                    h.FaceAlpha = 0.5;
                    set(h,'EdgeColor','none');


            plot(x,mean_dv,'k','LineWidth', 1.5)
            xlim([1 19])
            xticks([2:2:19])
            xticklabels({'-8','-6','-4','-2','0','+2','+4','+6','+8'})
            ylabel('∆velocity [mm/s]')
            xlabel('time from courtship bout start [sec]')
    
       end
%     if sum(sex.(groups{i}).(assays{j})==1)>0
%         %plot(x,nanmean(snip1'),'b','LineWidth',2)
%         %plot(x,nanmean(snip2'),'g','LineWidth',2)
%         %plot(x,nanmean(snip3'),'r','LineWidth',2)
%         %plot(x,nanmean(snip4'),'y','LineWidth',2)
%         plot(x,nanmean(vel_snippets'),'k','LineWidth',2)
%         title(name_specs(i));
%     end
    figure(7)
    velm.m=[];
    velm.f=[];
%     per_bout_mean=nanmean(vel_snippets);
%     t_test_vel{i,1}=per_bout_mean;
%     scatter(i*ones(length(per_bout_mean),1),per_bout_mean,'filled')
%     line([i-0.1 i+0.1],[nanmean(per_bout_mean) nanmean(per_bout_mean)],'Color','k','LineWidth',2)
    
    if sum(sex.(groups{i}).(assays{j})==1)>0
        vel_merge.m=vel_mean.(groups{i}).(assays{1})(sex.(groups{i}).(assays{1}) ==1);
        velm.m(1)=mean(vel_mean.(groups{i}).(assays{1})(sex.(groups{i}).(assays{1}) ==1));
    end
    
    if sum(sex.(groups{i}).(assays{j})==0)>0
        vel_merge.f=vel_mean.(groups{i}).(assays{1})(sex.(groups{i}).(assays{1}) ==0);
        velm.f(1)=mean(vel_mean.(groups{i}).(assays{1})(sex.(groups{i}).(assays{1}) ==0));
    end
    
    for v = 2:numAssays
        
       if sum(sex.(groups{i}).(assays{j})==1)>0 
        vel_merge.m= [vel_merge.m vel_mean.(groups{i}).(assays{v})(sex.(groups{i}).(assays{v}) ==1)];
        velm.m(v)=mean(vel_mean.(groups{i}).(assays{v})(sex.(groups{i}).(assays{v}) ==1));
       end
       
       if sum(sex.(groups{i}).(assays{j})==0)>0
        vel_merge.f= [vel_merge.f vel_mean.(groups{i}).(assays{v})(sex.(groups{i}).(assays{v}) ==0)];
        velm.f(v)=mean(vel_mean.(groups{i}).(assays{v})(sex.(groups{i}).(assays{v}) ==0));
       end
    end
    
    if sum(sex.(groups{i}).(assays{j})==0)>0
        hold on
        scatter(i*ones(length( vel_merge.f),1)', vel_merge.f, 'filled')
        plot( [i-0.1 i+0.1],[nanmean(vel_merge.f) nanmean(vel_merge.f)],'Color','k','LineWidth',2 ) 
    end
    
    if sum(sex.(groups{i}).(assays{j})==1)>0
        hold on
        scatter(i+numGroups*ones(length( vel_merge.m),1)', vel_merge.m, 'filled')
        plot( [i+numGroups-0.1 i+numGroups+0.1],[nanmean(vel_merge.m) nanmean(vel_merge.m)],'Color','k','LineWidth',2 ) 
    end
    
    ylim([0 15])
    xlim([0.5 numGroups*2+0.5])
    xticks([1:numGroups*2])
    
    xlabs=strcat ( 'females','-',name_specs,'-',inds);
    tmp = strcat ( 'males','-',name_specs,'-',inds);
    xlabs=[xlabs tmp];
    xticklabels(xlabs);
  
    ylabel('Mean velocity [mm/s]')
    vel(i).m=vel_merge.m;
    vel(i).f=vel_merge.f;
    vel(i).mm=velm.m;
    vel(i).mf=velm.f;
end

cumul_mint_blocks=struct(); %cumulative interation behavior of all males in a group
cumul_mapp_blocks=struct(); %cumulative approach behavior of all males in a group
cumul_mai_blocks=struct(); %cumulative interaction + approach behavior of all males in a group
cumul_mwing_blocks=struct();



%% ineraction plots
for i = 1:numGroups
     l=1;
   numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
   mcol=[0 0.8 0; 0 0.5 1];
   
  subplot(numGroups,1,i);
  m_mean_mean_lens=[];
  f_mean_mean_lens=[];
    
  cumul_mint_blocks.(groups{i})=[];
  cumul_mapp_blocks.(groups{i}) = [];
  cumul_mai_blocks.(groups{i})=[];
  cumul_mwing_blocks.(groups{i})=[];

   for j = 1:numAssays
      
       males=find(sex.(groups{i}).(assays{j})==1);
       females=find(sex.(groups{i}).(assays{j})==0);
       
        for in = 1:length(females)
            if ~isempty(females)
            fint_blocks=iblocks.(groups{i}).(assays{j}){females(in),1};
             fint_blocks=fint_blocks(fint_blocks(:,4)>=10,:);
              f_mean_lens=nanmean(fint_blocks(:,4));
              f_mean_mean_lens=[f_mean_mean_lens f_mean_lens]; 
            end
            
        end
        
        for in = 1:length(males)
            if ~isempty(males)
                mint_blocks=iblocks.(groups{i}).(assays{j}){males(in),1};

                mint_blocks=mint_blocks(mint_blocks(:,4)>=10,:);%remove interactions below 10 frames
                if ~isempty(ai_block_flip.(groups{i}).(assays{j}))
                    if ~isempty (ai_block_flip.(groups{i}).(assays{j}){males(in),1})
                        ai_tmp=ai_block_flip.(groups{i}).(assays{j}){males(in),1}(ai_block_flip.(groups{i}).(assays{j}){males(in),1}(:,4)>=10,:);
                        ai_tmp=ai_tmp(ai_tmp(:,1)<=assay_length,:);

                        cumul_mapp_blocks.(groups{i}) = [cumul_mapp_blocks.(groups{i}); ablocks.(groups{i}).(assays{j}){males(in),1}];
                        cumul_mai_blocks.(groups{i})=[cumul_mai_blocks.(groups{i}); size(ai_tmp,1)];
                        cumul_mwing_blocks.(groups{i})=[cumul_mwing_blocks.(groups{i}); wing_blocks.(groups{i}).(assays{j}){males(in),1}];

                        cumul_mint_blocks.(groups{i})=[cumul_mint_blocks.(groups{i}); mint_blocks];
                    else
                        cumul_mai_blocks.(groups{i})=[cumul_mai_blocks.(groups{i}); 0];

                    end
                else
                    cumul_mai_blocks.(groups{i})=[cumul_mai_blocks.(groups{i}); 0];
                end
                
                m_mean_lens=nanmean(mint_blocks(:,4));
          
                m_mean_mean_lens=[m_mean_mean_lens m_mean_lens];

            

            end

           
           
           
           

           l=l+1;
        end
   end


   plot([0 60], [nanmean(m_mean_mean_lens) nanmean(m_mean_mean_lens)],'Color','r','LineStyle','--','LineWidth',2)
    xlim([0 31])
    
    figure(9)
    scatter(i*ones(length(m_mean_mean_lens),1),m_mean_mean_lens,'filled')
    hold on
    plot([i-0.2 i+0.2],[nanmean(m_mean_mean_lens) nanmean(m_mean_mean_lens)],'Color','k','LineWidth',2)
    scatter((i+numGroups)*ones(length(f_mean_mean_lens),1),f_mean_mean_lens,'filled')
    plot([i+numGroups-0.2 i+numGroups+0.2],[nanmean(f_mean_mean_lens) nanmean(f_mean_mean_lens)],'Color','k','LineWidth',2)
    ylim([0 500])
    xlim([0 numGroups*2+1])
    ylabel('Mean interaction length')
    xticks([1:numGroups*2])
    
        name_specs={};
    if size(chem,2)>1
        name_specs=chem;
    else 
        if size(housing,2)>1
           name_specs=housing;
        else
            
            if size(inds,2)>1
                name_specs=inds;
                
            else
                    
                  if size(spec,2)>1
                        name_specs = spec;
                  else
                        name_specs=age;
                  end
            end
        end
    end
    xlabs=strcat ( 'males','-',name_specs);
    tmp = strcat ( 'females','-',name_specs);
    xlabs=[xlabs tmp];
    xticklabels(xlabs);
  
    cumul_mapp_blocks.(groups{i})=cumul_mapp_blocks.(groups{i})(cumul_mapp_blocks.(groups{i})(:,4)>=10,:);
    
        cumul_mapp_blocks.(groups{i})=cumul_mapp_blocks.(groups{i})(cumul_mapp_blocks.(groups{i})(:,1)<assay_length,:);
        cumul_mint_blocks.(groups{i})=cumul_mint_blocks.(groups{i})(cumul_mint_blocks.(groups{i})(:,1)<assay_length,:);
        cumul_mwing_blocks.(groups{i})=cumul_mwing_blocks.(groups{i})(cumul_mwing_blocks.(groups{i})(:,1)<assay_length,:);

end

%% time to courtship initiation
figure()
ttci=struct();
ttii=struct(); %time to interaction initiation -- female courtship initiation pendant
wing_per_ai_bout=struct();
vel_per_ai_bout=struct();
vel_per_ai_boutF=struct(); %only works for 1m1f assays. Disregard for any other

for i = 1:numGroups
    
       numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
    
   wing_per_ai_bout.(groups{i}).means=[];
   vel_per_ai_bout.(groups{i}).means=[];
   vel_per_ai_boutF.(groups{i}).means=[];
   ttci.(groups{i}).all=[];
   for j = 1:numAssays
      
       
       
       for flieses = 1:size(ai_maps.(groups{i}).(assays{j}),2)
           if ~isempty(ai_block_flip.(groups{i}).(assays{j})) && size(ai_block_flip.(groups{i}).(assays{j}),1)>=flieses
               
               if size(ai_block_flip.(groups{i}).(assays{j}){flieses},1)>1

                    ai_block_flip.(groups{i}).(assays{j}){flieses} = sortrows(ai_block_flip.(groups{i}).(assays{j}){flieses},1) ;
               end
           end
       end
       
       %k = find(sex.(groups{i}).(assays{j})==1);
       %wing_per_ai_bout.(groups{i}).(assays{j}).means=nan(k,1);
       
       %males to females
       for k = find(sex.(groups{i}).(assays{j})==1) %counts all males if 1
           if isempty(ai_block_flip.(groups{i}).(assays{j}))
                int_occur=0;
           else
                int_occur=size(ai_block_flip.(groups{i}).(assays{j}){k},1);
           end
           
           
           
           if int_occur>0
                long_enough = ai_block_flip.(groups{i}).(assays{j}){k}(ai_block_flip.(groups{i}).(assays{j}){k}(:,4)>fps,:);
                to_females=long_enough(sex.(groups{i}).(assays{j})(long_enough(:,3))==0,:); %interactions with females if 0

                if ~isempty(to_females)

                    initial_length_required=fps*5; %first interaction length minimum to call initiation
                    if sum(to_females(:,4)>=initial_length_required)>=1
                        first=find(to_females(:,4)>=initial_length_required);

                        ttci.(groups{i}).(assays{j})(k)= to_females(first(1),1);
                        
                        %calculate wing extension frames per bout frames for
                        %each frame for all males
                        % k is male id
                        tmp_block=[];
                        tmp_vel_block=[];
                        for ai =1:size(long_enough,1)
                            wing_frames=nnz(wing_map.(groups{i}).(assays{j})(long_enough(ai,1):long_enough(ai,2),k)>0);
                            mean_male_vel=nanmean(velocity.(groups{i}).(assays{j})(long_enough(ai,1):long_enough(ai,2),k));
                            if k==1
                                fem=2;
                            else
                                fem=1;
                            end
                            mean_female_vel=nanmean(velocity.(groups{i}).(assays{j})(long_enough(ai,1):long_enough(ai,2),fem));
                            
                            wing_per_block=0;
                            if wing_frames
                               wing_per_block= wing_frames/long_enough(ai,4);
                            end
                            tmp_block=[tmp_block; long_enough(ai,:), wing_frames, wing_per_block];
                            tmp_vel_block=[tmp_vel_block; long_enough(ai,:), mean_male_vel, mean_female_vel];
                        end
                        wing_per_ai_bout.(groups{i}).(assays{j}){k,1}=tmp_block;
                        wing_per_ai_bout.(groups{i}).means=[wing_per_ai_bout.(groups{i}).means; nanmean(tmp_block(:,6))];
                        
                        vel_per_ai_bout.(groups{i}).(assays{j}){k,1}=tmp_vel_block;
                        vel_per_ai_bout.(groups{i}).means=[vel_per_ai_bout.(groups{i}).means; nanmean(tmp_vel_block(:,5))];
                        vel_per_ai_boutF.(groups{i}).means=[vel_per_ai_boutF.(groups{i}).means; nanmean(tmp_vel_block(:,6))];
                    else
                        ttci.(groups{i}).(assays{j})(k)= assay_length;
                    end
                else
                    ttci.(groups{i}).(assays{j})(k)= assay_length;
                end
           else
               ttci.(groups{i}).(assays{j})(k)= assay_length;
               
               wing_per_ai_bout.(groups{i}).(assays{j}){k,1}=[];
               wing_per_ai_bout.(groups{i}).means=[wing_per_ai_bout.(groups{i}).means; nan()];
               vel_per_ai_bout.(groups{i}).(assays{j}){k,1}=[];
               
              vel_per_ai_bout.(groups{i}).means=[vel_per_ai_bout.(groups{i}).means; nan()];
               vel_per_ai_boutF.(groups{i}).means=[vel_per_ai_bout.(groups{i}).means; nan()];

           end
       end
       
       if sum(sex.(groups{i}).(assays{j})==1)>0 && sum(sex.(groups{i}).(assays{j})==0)>0
        ttci.(groups{i}).(assays{j})(ttci.(groups{i}).(assays{j})==0)=[];
       end
       %females to males or anyone - switch exists
       
       for k = find(sex.(groups{i}).(assays{j})==0) %counts all males if 1
           if isempty(ai_block_flip.(groups{i}).(assays{j})) 
                int_occur=0;
           else
               if size(ai_block_flip.(groups{i}).(assays{j}),1)<k
                int_occur=0;
               else
                int_occur=size(ai_block_flip.(groups{i}).(assays{j}){k},1);
               end
           end
           
           
           if int_occur>0
                long_enough = ai_block_flip.(groups{i}).(assays{j}){k}(ai_block_flip.(groups{i}).(assays{j}){k}(:,4)>fps,:);
                %to_males=long_enough(sex.(groups{i}).(assays{j})(long_enough(:,3))==1,:); %interactions with females if 0
                to_males=long_enough; %this is set to start counting as soon as a female initiates interaction with another individual independent of sex
                
                if ~isempty(to_males)
                    ttii.(groups{i}).(assays{j})(k)= to_males(1,1);
                    
                     
                    
                else
                    ttii.(groups{i}).(assays{j})(k)= assay_length;
                end
           else
               ttii.(groups{i}).(assays{j})(k)= assay_length;
               
              
           end
       end
       
       if sum(sex.(groups{i}).(assays{j})==0)>0
        ttii.(groups{i}).(assays{j})(ttii.(groups{i}).(assays{j})==0)=[];
       end
       
       if sum(sex.(groups{i}).(assays{j})==1)>0 && sum(sex.(groups{i}).(assays{j})==0)>0
           col='b';
           hold on
           scatter(i*ones(sum(sex.(groups{i}).(assays{j})),1),ttci.(groups{i}).(assays{j}),[],col,'filled','jitter', 'on', 'jitterAmount', 0.1)
           xlim([0.5 numGroups+0.5])
           ylim ([0 assay_length])
           
           mttci=mean(ttci.(groups{i}).(assays{j}));
           ttci.(groups{i}).means(j)=mttci;
           ttci.(groups{i}).all=[ttci.(groups{i}).all; ttci.(groups{i}).(assays{j})];
           
       end
   end
   allvals=[];
   
   if sum(sex.(groups{i}).(assays{j})==1)>0 && sum(sex.(groups{i}).(assays{j})==0)>0
       for ij = 1:numAssays
        allvals=[allvals ttci.(groups{i}).(assays{ij})] ;
       end
   end
   hold on
   plot([i-0.2 i+0.2],[nanmean(allvals) nanmean(allvals)],'Color','k','LineWidth',2)
   
   xticks(1:numGroups)
    xlabs=strcat ( 'males','-',name_specs);
    xticklabels(xlabs);
    ylabel({'courtship initiation [min]'})
    ytick_labs=0:fps*60*10:assay_length;
    yticks(ytick_labs)
    yticklabels({'0','10', '20','30', '40'})
    ylim([0 assay_length])
    
%     ytick_labs=0:fps*60:fps*60*5;
%     yticks(ytick_labs)
%     yticklabels({'0','1','2','3','4','5'})
%     ylim([0 5*fps*60])
    
    wing_per_ai_bout.(groups{i}).means(wing_per_ai_bout.(groups{i}).means>100)=[];
end



%% plot Courtship index and courtship vigor
figure(11)

for i = 1:numGroups
    
       numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
   
   vigor.(groups{i}).all=[];
   f_interactionVigor.(groups{i}).all=[];
   
   cuin.(groups{i}).all=[];
   for j = 1:numAssays
      
       
       malenum=1;
       femalenum=1;
       for flieses = 1:size(ai_block_flip.(groups{i}).(assays{j}),1)
       %for flieses = 1:numflies
           vigor_length=fps*60*10;
           
           if sex.(groups{i}).(assays{j})(flieses)==1
               %male?
              if ttci.(groups{i}).(assays{j})(malenum)<assay_length
               %did court?
                                                
                    accum_interaction_time=0;
                    
                    for row = 1:size(ai_block_flip.(groups{i}).(assays{j}){flieses},1)
                       
                            foc_row=ai_block_flip.(groups{i}).(assays{j}){flieses}(row,:);
                        if foc_row(1)>=ttci.(groups{i}).(assays{j})(malenum)
                           %row after determined CI?
                            if foc_row(1)<= ttci.(groups{i}).(assays{j})(malenum)+vigor_length
                               %row within 10 minutes of courtship init?
                               
                               if foc_row(4)>10
                                  %interaction longer than 10 frames?
                                  
                                  accum_interaction_time=accum_interaction_time+foc_row(4);
                                  %then add length of interaction to
                                  %overall interaction within 10 minute
                                  %time window post initiation.
                                   
                               end
                            end
                            
                        end
                        
                    end
                    
                    
                    if ~isempty(cop.(groups{i}).(assays{j}).blocks{flieses})
                        %if individual copulates
                        if cop.(groups{i}).(assays{j}).blocks{flieses}(1,1)<ttci.(groups{i}).(assays{j})(malenum)+vigor_length
                            %start copulating earlier than 10 minutes after
                            %courtship initiation?
                           vigor_length =  ttci.(groups{i}).(assays{j})(malenum)+ cop.(groups{i}).(assays{j}).blocks{flieses}(1,1);
                           %adjust stretch of time used to calculate vigor
                        end
                    end
                  vigor.(groups{i}).(assays{j})(flieses)=accum_interaction_time/vigor_length;
                  
                    
              else
                    vigor.(groups{i}).(assays{j})(flieses)=0;
              end
              malenum=malenum+1;
              f_interactionVigor.(groups{i}).(assays{j})(flieses)=nan();
              
           else
               %female
               
               %%%%%%%%
               %%%%%%%%
               %%%%%%%%
               
               if ttii.(groups{i}).(assays{j})(femalenum)<assay_length
               %did court?
                                                
                    accum_interaction_time=0;
                    
                    for row = 1:size(ai_block_flip.(groups{i}).(assays{j}){flieses},1)
                       
                            foc_row=ai_block_flip.(groups{i}).(assays{j}){flieses}(row,:);
                        if foc_row(1)>=ttii.(groups{i}).(assays{j})(femalenum)
                           %row after determined CI?
                            if foc_row(1)<= ttii.(groups{i}).(assays{j})(femalenum)+vigor_length
                               %row within 10 minutes of courtship init?
                               
                               if foc_row(4)>10
                                  %interaction longer than 10 frames?
                                  
                                  accum_interaction_time=accum_interaction_time+foc_row(4);
                                  %then add length of interaction to
                                  %overall interaction within 10 minute
                                  %time window post initiation.
                                   
                               end
                            end
                            
                        end
                        
                    end
                    
                    
                    if ~isempty(cop.(groups{i}).(assays{j}).blocks{flieses})
                        %if individual copulates
                        if cop.(groups{i}).(assays{j}).blocks{flieses}(1,1)<ttii.(groups{i}).(assays{j})(femalenum)+vigor_length
                            %start copulating earlier than 10 minutes after
                            %courtship initiation?
                           vigor_length =  ttii.(groups{i}).(assays{j})(femalenum)+ cop.(groups{i}).(assays{j}).blocks{flieses}(1,1);
                           %adjust stretch of time used to calculate vigor
                        end
                    end
                  f_interactionVigor.(groups{i}).(assays{j})(flieses)=accum_interaction_time/vigor_length;
                  
                    
              else
                    f_interactionVigor.(groups{i}).(assays{j})(flieses)=0;
              end
              femalenum=femalenum+1;
               
               %%%%%%%
               %%%%%%%
               %%%%%%%
               
               
               vigor.(groups{i}).(assays{j})(flieses)=nan();
           end
           
       end
       hold on
       scatter((i+numGroups)*ones(sum(sex.(groups{i}).(assays{j})==1),1), vigor.(groups{i}).(assays{j})(sex.(groups{i}).(assays{j})==1), 'filled','jitter', 'on', 'jitterAmount', 0.1)
       xlim([0.5 numGroups*2+0.5])
       ylim([0 1])
       
       vigor.(groups{i}).all=[vigor.(groups{i}).all vigor.(groups{i}).(assays{j})(sex.(groups{i}).(assays{j})==1)];
       if length(f_interactionVigor.(groups{i}).(assays{j})) < length(sex.(groups{i}).(assays{j})==0)
        f_interactionVigor.(groups{i}).all=[f_interactionVigor.(groups{i}).all 0];
       else
        
        f_interactionVigor.(groups{i}).all=[f_interactionVigor.(groups{i}).all f_interactionVigor.(groups{i}).(assays{j})(sex.(groups{i}).(assays{j})==0)];   
       end
       
       assay_sex=sex.(groups{i}).(assays{j});
       assay_total=ii.(groups{i}).(assays{j}).total;
       males=find(assay_sex==1); %should be 1 not 0
       
       cuin.(groups{i}).all=[cuin.(groups{i}).all assay_total(5,males)];%should be males not females

       
       scatter(i*ones(length(males),1)',assay_total(5,males),'filled','jitter', 'on', 'jitterAmount', 0.1)

   end
   
   hold on
   
   plot([i-0.2 i+0.2],[nanmean(cuin.(groups{i}).all) nanmean(cuin.(groups{i}).all)],'Linewidth',3,'Color','k')
    plot([i+numGroups-0.2 i+numGroups+0.2],[nanmean(vigor.(groups{i}).all) nanmean(vigor.(groups{i}).all)],'Linewidth',3,'Color','k')
    ylabel('CI / Vigor')
    xticks([1:numGroups*2])
    
        xlabs=strcat ( 'CI','-',name_specs);
    tmp = strcat ( 'vigor','-',name_specs);
    xlabs=[xlabs tmp];
    xticklabels(xlabs);
    
    
    %xticklabels({'CI-agar','CI-acv','vigor-agar','vigor-acv'});
    
    cop.(groups{i}).vig_m=vigor.(groups{i}).all;
    cop.(groups{i}).vig_f=f_interactionVigor.(groups{i}).all;
end

%% plot probability of courtship over time

figure()
colsRGB=[0 0.6 1; 1 0.6 0];
meanSaver=[];
for i = 1:numGroups
   

   numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));

    if mod(i,2)
        sp=1;
    else
        sp=2;
    end

   %subplot(2,1,sp)
   hold on
   concat=[];
   for j = 1:numAssays
        
       numMales=sum(sex.(groups{i}).(assays{j}));
       maleIds=find(sex.(groups{i}).(assays{j})==1);
       femaleIds=find(sex.(groups{i}).(assays{j})==0);
       for k= 1:numMales

           courtInit=ttci.(groups{i}).(assays{j})(k);
           courtEnd=assay_length;
           if cop.(groups{i}).(assays{j}).ind(k)==1
                cb=cop.(groups{i}).(assays{j}).blocks{k,1};
                courtEnd=cb(1)-1;                           
           end

           
           courtBeh=focal_maps.(groups{i}).(assays{j})(courtInit:courtEnd,maleIds(k));

           for whoto = maleIds
                %courtBeh(courtBeh==whoto)=0;
           end
           courtBeh(courtBeh>1)=1;
           if length(courtBeh)<assay_length
               difference=assay_length-length(courtBeh);
                filler=[courtBeh; nan(difference,1)];
           end
           concat=[concat filler];
           
%            hold on
%            
%            plot(smoothdata(filler,'movmean',fps*60),'Color',[0.2 0.2 0.2])
            
           
           
          
       end

   end
   concat(1,find(concat(1,:)==0))=nan();
   res(sp).cb=concat;
   
   means=nanmean(concat,2);
   res(sp).means=means;

   sds=nanstd(concat,[],2)/sqrt(numAssays*numMales);
   
   
   min_inc=1; %how long?
   res(sp).sds=sds;
           x=1:fps*60*min_inc; %row vector
            curve1 = smoothdata(means(1:fps*60*min_inc) + sds(1:fps*60*min_inc),'movmean',fps*2,'omitnan'); %column vector
            curve2 = smoothdata(means(1:fps*60*min_inc) - sds(1:fps*60*min_inc),'movmean',fps*2,'omitnan'); %column
            x2 = [x, fliplr(x)];
            inBetween = [curve1; flipud(curve2)];

            if mod(k,2)
                colType=1;
            else
                colType=2;
            end

         
                if mod(i,2)>0
                    colType=1;
                else
                    colType=2;
                    
                end

                colPL='k';

                
            %ha=fill(x2, inBetween', colsRGB(colType,:));
            %ha.FaceAlpha = 0.5;
            %set(ha,'EdgeColor','none');
            fill(x2, inBetween', colsRGB(colType,:),'FaceAlpha',0.5, 'EdgeColor', 'none')
            xticks(0:fps*60:fps*60*10)
            xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
            xlabel('Time post courtship initiation (min)')
            hold on
            tmpp=smoothdata(means,'movmean',fps*2,'includenan');
            plot(tmpp,colPL)

            %yline(0.5,'--r')
            ylim([0 1])
            ylabel('p(Courtship)')
            xlim([0 fps*60*min_inc])
            %AUC???
            meanSaver=[meanSaver tmpp];


end



%% calculate pairwise distances
numGroups=length(groups);
steppo=1350/38;
exclude=steppo*3; %nan values below 3mm (~1 fly length)

normal=0 + 1350*rand(size(f));

colors={'b','g'};
edges=0:10:1350;

for i =1:numGroups
       numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
   histEuc.(groups{i})=[];
   
   for j = 1:numAssays
        euc_dist.(groups{i}).(assays{j}).eucmat=[];
        %for ind1=1:numflies-1
        for ind1=find(sex.(groups{i}).(assays{j})==1) %males

            %for ind2=ind1+1:numflies
            for ind2=find(sex.(groups{i}).(assays{j})==0) %females

                euc_dist.(groups{i}).(assays{j}).eucmat=[euc_dist.(groups{i}).(assays{j}).eucmat euc_dist.(groups{i}).(assays{j}).eucCell{ind1,ind2}];

            end
        end
        f=euc_dist.(groups{i}).(assays{j}).eucmat;
        f2=f;
        f2(f2<exclude)=nan();






       
       h1=histcounts(f2,edges,'Normalization','probability');
       histEuc.(groups{i})=[histEuc.(groups{i}); h1];
   end

end


%% plot pairwise distances when courting (fig 1) and for all data (fig 2)
edges=0:10:1350;
colsRGB=[0 0.6 1; 1 0.6 0];
figure
for i =1:numGroups
    
    numAssays = length(fieldnames(ii.(groups{i})));
    assays = fieldnames(ii.(groups{i}));
    
    for j=1:numAssays
        court_dist.(groups{i})(j).accum=[];
        court_dist.(groups{i})(j).allnocop=[];

        femaleIDcheck=find(sex.(groups{i}).(assays{j})==0);
        for k=find(sex.(groups{i}).(assays{j})==1)
            tmp=struct();
            %filter only those times where males actually court a female
            blocks=ai_blocks.(groups{i}).(assays{j}){k,1};
            blocks_f=blocks(blocks(:,4)>fps/2,:);
            blocks_fi=blocks_f(ismember(blocks_f(:,3),femaleIDcheck),:);
            

            for l=1:size(blocks_fi,1)

                start=blocks_fi(l,1);
                fin=blocks_fi(l,2);
                target=blocks_fi(l,3);
                extract=nan();
                if k<target %when male id is smaller than female id
                    extract=euc_dist.(groups{i}).(assays{j}).eucCell{k,target}(start:fin,:);
                else
                    extract=euc_dist.(groups{i}).(assays{j}).eucCell{target,k}(start:fin,:);
                end
                
                court_dist.(groups{i})(j).accum=[court_dist.(groups{i})(j).accum; extract];
                

            end

            %filer all pairwise distances between males and females but exclude
            %matings
            
            for fem=femaleIDcheck
                    
                if k<fem
                    tmp(fem).pdist=euc_dist.(groups{i}).(assays{j}).eucCell{k,fem};
                else
                    tmp(fem).pdist=euc_dist.(groups{i}).(assays{j}).eucCell{fem,k};

                end
            end

            if size(cop.(groups{i}).(assays{j}).blocks{k,1},1)>0
                    
                mblock=cop.(groups{i}).(assays{j}).blocks{k,1};

                for l=1:-1:size(mblock,1) %removes all pdists to all females when a male copulates
                    
                    start=mblock(l,1);
                    fin=mblock(l,2);
                    %target=mblock(l,3);
                    for f=femaleIDcheck
                        tmp(f).pdist(start:fin)=nan();
                    end
                end

            end

            allPsmale=vertcat(tmp.pdist);
            nonanall= allPsmale(~isnan(allPsmale));
            court_dist.(groups{i})(j).allnocop=[court_dist.(groups{i})(j).allnocop; nonanall]; %not the best naming but this saves all pdist of all males to all females while not copulating
        end

        h1=histcounts(court_dist.(groups{i})(j).accum,edges,'Normalization','probability');
        court_dist.(groups{i})(j).hist=h1;

        h2=histcounts(court_dist.(groups{i})(j).allnocop,edges,'Normalization','probability');
        court_dist.(groups{i})(j).allHist=h2;
    end
    all=vertcat(court_dist.(groups{i}).hist);
    m=mean(all);
    maxCourt=find(m==max(m));
    sds=std(all)/sqrt(numAssays);
   saver(i).courtOnly=all;
    x=1:size(all,2);
    curve1 = m + sds; %column vector
    curve2 = m - sds; %column
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    hold on
    
    fill(x2, inBetween, colsRGB(i,:),'FaceAlpha',0.5, 'EdgeColor', 'none')
    plot(m,'Color',colsRGB(i,:),'LineWidth',2)
end
step=135/38; % 1350/38/10 for real but same
ticks=0:step*3:135;

xlabel('distance [mm]')
ylabel('probability')
ylim([0 0.3])

xticks(ticks)
xticklabels(num2cell(0:3:38))  
%xlim([0 ticks(3)])

%pairwise distances total
figure

for i =1:numGroups
        all=vertcat(court_dist.(groups{i}).allHist);
        saver(i).dist=all;
    m=mean(all);
    sds=std(all)/sqrt(numAssays);
   
    x=1:size(all,2);
    curve1 = m + sds; %column vector
    curve2 = m - sds; %column
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    hold on
    
    fill(x2, inBetween, colsRGB(i,:),'FaceAlpha',0.5, 'EdgeColor', 'none')
    plot(m,'Color',colsRGB(i,:),'LineWidth',2)
    xline(maxCourt,'--k','LineWidth',2)

end
step=135/38; % 1350/38/10 for real but same
ticks=0:step*3:135;

xlabel('distance [mm]')
ylabel('probability')
ylim([0 0.04])

xticks(ticks)
xticklabels(num2cell(0:3:38))  
%xlim([0 ticks(3)])


%% plot 2D distribution heatmap
% Define the number of colors for each transition segment
nColors = 256*2;
nSegment = round(nColors / 2);  % Half for white-to-green, half for green-to-purple

% Transition from white to dark green
whiteToGreen = [linspace(1, 0, nSegment)', linspace(1, 0.5, nSegment)', linspace(1, 0, nSegment)'];

% Transition from dark green to purple
greenToPurple = [linspace(0, 0.5, nSegment)', linspace(0.5, 0, nSegment)', linspace(0, 0.5, nSegment)'];


% Combine the segments to form the full colormap
cmap = [whiteToGreen; greenToPurple];

figure()
for i=1:numGroups

    numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
   num_bins=50;
   totalCounts=zeros(num_bins,num_bins);
   for j = 1:numAssays
        
       allX=positions.(groups{i})(j).x(:);
       allY=positions.(groups{i})(j).y(:);

       
       x_edges=linspace(0, 1500, num_bins +1); %from 0 to 1500; 201 bins
       y_edges=linspace(0, 1500, num_bins + 1); %from 0 to 1500; 201 bins

       [counts, ~, ~] = histcounts2(allX, allY, x_edges, y_edges);
       
       positions.(groups{i})(j).counts=counts;
       positions.(groups{i})(j).totalC=sum(counts(:));
       totalCounts=totalCounts+counts;
     

   end  
    numC=sum(totalCounts(:));
    normTotal=totalCounts/numC;
    heatM(i).normTotal=normTotal;
end

globalMax=max(max([heatM.normTotal]));
globalMin=min(min([heatM.normTotal]));

for i=1:numGroups
            subplot(2,1,i)
        normTotal=heatM(i).normTotal;
        saver(i).dist=normTotal;
        imagesc(x_edges,y_edges, normTotal)
        colormap(cmap);
        colorbar;
        caxis([globalMin globalMax]);

end


%% plot time spent immobile and velocity 
% 1. how much time spent sitting?
% 2. of the time spent moving, how fast? (mean)
% 3. exclusively pre-mating
cutoff=0.05; %sitting cutoff

figure();
mobility=struct();
for i =1:numGroups
       numAssays = length(fieldnames(ii.(groups{i})));
   assays = fieldnames(ii.(groups{i}));
   
   %values each individual
   mobility.(groups{i}).males=[];
   mobility.(groups{i}).females=[];
   mobility.(groups{i}).mmales=[];
   mobility.(groups{i}).mfemales=[];

    %individual values
   sittingMassay=[]; %assay means
   sittingFassay=[];
   speedMassay=[];
   speedFassay=[];
   
   for j = 1:numAssays
   sittingMtmp=[];
   sittingFtmp=[];
   sittingFind=[]; %individual values
   sittingMind=[];
   speedFind=[];
   speedMind=[];
   
        velPrecise.(groups{i})(j).vel={};
        mobile.(groups{i})(j).binary=nan(assay_length,10);  
       for k =1:size(velocity.(groups{i}).(assays{j}),2)
            mated=0;
            max=assay_length;%size(velocity.(groups{i}).(assays{j}),1);
           %check if mated
            if cop.(groups{i}).(assays{j}).ind(k)==1
                mated=1;
                b1=cop.(groups{i}).(assays{j}).blocks{k};
                max=b1(1)-1;
            end

            mobile.(groups{i})(j).binary(1:max,k)=velocity.(groups{i}).(assays{j})(1:max,k)>=cutoff;
            sitting=sum(velocity.(groups{i}).(assays{j})(1:max,k)<cutoff);
            total=max;
            rel=sitting/total;
            moving=velocity.(groups{i}).(assays{j})(velocity.(groups{i}).(assays{j})(1:max,k)>=cutoff,k);
            speed=mean(moving);
            if sex.(groups{i}).(assays{j})(k)==1
                sittingMind=[sittingMind rel];
                sittingMtmp=[sittingMtmp rel];
                speedMind=[speedMind speed];

            else
                sittingFind=[sittingFind rel];
                sittingFtmp=[sittingFtmp rel];
                speedFind=[speedFind speed];
            end

            velPrecise.(groups{i})(j).vel{k}=velocity.(groups{i}).(assays{j})(1:max,k);
            velPrecise.(groups{i})(j).ai{k}=ai_maps.(groups{i}).(assays{j})(1:max,k);
            velPrecise.(groups{i})(j).movingTT{k}=velocity.(groups{i}).(assays{j})(1:max,k)>=0.06;
       end
    sittingMassay=[sittingMassay mean(sittingMtmp)];
    sittingFassay=[sittingFassay mean(sittingFtmp)];

    speedMassay=[speedMassay mean(speedMind)];
    speedFassay=[speedFassay mean(speedFind)];
    
    %individuals plotting
    subplot(2,2,1)
    hold on
    scatter(ones(length(sittingMind),1)*i,sittingMind,[],'m', 'Filled','jitter','on','jitterAmount', 0.1)

    mobility.(groups{i}).males=[mobility.(groups{i}).males sittingMind];

    mobility.(groups{i}).mmales=[mobility.(groups{i}).mmales mean(sittingMind)];

    subplot(2,2,2)
    hold on
    scatter(ones(length(speedMind),1)*i,speedMind,[],'m', 'Filled','jitter','on','jitterAmount', 0.1)
    
    subplot(2,2,3)
    hold on
    scatter(ones(length(sittingFind),1)*i,sittingFind,[],'b', 'Filled','jitter','on','jitterAmount', 0.1)
    mobility.(groups{i}).females=[mobility.(groups{i}).females sittingFind];

    mobility.(groups{i}).mfemales=[mobility.(groups{i}).mfemales mean(sittingFind)];

    subplot(2,2,4)
    hold on
    scatter(ones(length(speedFind),1)*i,speedFind,[],'b', 'Filled','jitter','on','jitterAmount', 0.1)

   end
    subplot(2,2,1)
    hold on
    plot([i-0.25 i+0.25], [mean(sittingMassay) mean(sittingMassay)],'k','LineWidth', 2)
    mobility.(groups{i}).mmean=mean(sittingMassay);

    subplot(2,2,2)
    hold on

    plot([i-0.25 i+0.25], [mean(speedMassay) mean(speedMassay)],'k','LineWidth', 2)
    
    subplot(2,2,3)
    hold on
    plot([i-0.25 i+0.25], [mean(sittingFassay) mean(sittingFassay)],'k','LineWidth', 2)
    mobility.(groups{i}).fmean=mean(sittingFassay);

    subplot(2,2,4)
    hold on
    plot([i-0.25 i+0.25], [mean(speedFassay) mean(speedFassay)],'k','LineWidth', 2)

end
labs={'fraction of time immobile','mean velocity if moving [mm/sec]','fraction of time immobile','mean velocity if moving [mm/sec]'};
for pl=1:4
    subplot(2,2,pl)
    xlim([0.5 numGroups+0.5])
    xticks([1:numGroups])
    xticklabels({'agar','acv','agar','acv'}) 
    ylabel(labs{pl})
    if pl==1 || pl==3
        ylim([0 1])
    else
        ylim([0 15])
    end
end


%% calculate probability of fly moving given other flies are moving and random draw (as baseline)

win_width=fps*2; %2 seconds
for i=1:numGroups

    for j=1:numel(mobile.(groups{i}))
    
        for k=1:10
            others_moving=[];
            
            target=mobile.(groups{i})(j).binary(:,k);
            others=mobile.(groups{i})(j).binary;others(:,k)=[];

            %for l=win_width+1:assay_length
            for l=1:assay_length-(win_width+1)

                %win=others(l-win_width:l,:);
                win=others(l:l+win_width,:);

                if target(l)==1
                    %check if other flies move within 5s before socal fly moves
                    check=sum(sum(win,'omitnan'));

                  
                    others_moving=[others_moving check];


                end

            end
            mobile.(groups{i})(j).othersMoving{:,k}=others_moving;
            times_others_moving=sum(others_moving>=1);
            times_focal_moving=length(others_moving);
            mobile.(groups{i})(j).othersFreq=times_others_moving/times_focal_moving;

             reps=[];
             for rep=1:10
                others_movingP=[];
                %now random draw, not only when focal fly moves
                %numbers = randperm(assay_length-(win_width+1), times_focal_moving)+win_width;
                numbers = randperm(assay_length-(win_width+1), times_focal_moving);

                for l=numbers
                    
                    %win=others(l-win_width:l,:);
                    win=others(l:l+win_width,:);
                    check=sum(sum(win,'omitnan'));
                    others_movingP=[others_movingP check];
                end
                times_others_moving=sum(others_movingP>=1);
                reps=[reps times_others_moving/times_focal_moving];
                
             end
            mobile.(groups{i})(j).othersFreqP=mean(reps);
            mobile.(groups{i})(j).othersP_reps=reps;
        end
    end


end
%% plot
figure; hold on
agF=[mobile.groupere_5m5f___agar.othersFreq];
agFp=[mobile.groupere_5m5f___agar.othersFreqP];
acF=[mobile.groupere_5m5f___acv.othersFreq];
acFp=[mobile.groupere_5m5f___acv.othersFreqP];

scatter(ones(size(agF)),agF,[],'filled','jitter', 'on')
plot([0.75 1.25],[mean(agF) mean(agF)],'k','linewidth',2)

scatter(2*ones(size(agFp)),agFp,[],'filled','jitter', 'on')
plot([1.75 2.25],[mean(agFp) mean(agFp)],'k','linewidth',2)

scatter(3*ones(size(acF)),acF,[],'filled','jitter', 'on')
plot([2.75 3.25],[mean(acF) mean(acF)],'k','linewidth',2)

scatter(4*ones(size(acFp)),acFp,[],'filled','jitter', 'on')
plot([3.75 4.25],[mean(acFp) mean(acFp)],'k','linewidth',2)

xlim([.5 4.5])
xticks(1:4)
xticklabels({'p agar','p agar ran','p acv','p acv ran'})


%% velocity of others around approach


bin_edges=0:fps*60:fps*60*40; % to count interactions per minute over time

   figure()


%%
function [ai,foc]=mapLoader(int_map,app_map,assay_length,numflies)

    ai = zeros(assay_length,numflies);
    foc = ai;

    for i = 1:numflies
        to= assay_length;

            %merge interaction and approach map to 'holistic interaction map'
        foc(1:assay_length,i)=int_map(1:assay_length,i);

        for j = 1: to

            if(app_map(j,i)>0)
                if foc(j,i)==0
                    foc(j,i)=app_map(j,i); 
                end
                ai(j,app_map(j,i))=i;
            end

            if(int_map(j,i)>0)
                      ai(j,int_map(j,i))=i;
            end
        end        
    end

end

function [ii,ii_m,ii_f] = interactionCounter(mat_blocks,foc_map,assay_length,numflies,sex,fps)

ii=zeros(5,numflies);
ii_m=ii;
ii_f=ii;


    for i = 1:numflies
        fin=assay_length;

        %ii only up until first mating
        if ~isempty(mat_blocks{i,1})
            if mat_blocks{i,1}(1,1)-1 < fin
                fin = mat_blocks{i,1}(1,1)-1;
            end
        end
        
        %total ii, ii_m, and ii_f for all individuals are saved in 5th row
        total_interactions = nnz(foc_map(1:fin,i)>0);
        ii(5,i)=total_interactions/fin;
        ii_m(5,i)=nnz(sex(foc_map(foc_map(1:fin,i)>0,i))==1)/fin;
        ii_f(5,i)=nnz(sex(foc_map(foc_map(1:fin,i)>0,i))==0)/fin;
        
        %bins=[1 fps*60*10; (fps*60*10)+1 fps*60*20; (fps*60*20)+1 fps*60*30;(fps*60*30)+1 fps*60*40];
        quart=round(assay_length/4);
        bins=[1 quart; quart+1 quart*2; (quart*2)+1 quart*3; (quart*3)+1 assay_length];
        %bin_length=fps*60*10;
        bin_length=quart;
        
        
        for t = 1:4
            bstart=bins(t,1);
            bend=bins(t,2);
            
            foc_map_cur=foc_map(bstart:bend,i);
            
 %fprintf('init start: %d, end: %d, first cop: %d\n',bstart,bend,fin)
            
            
            bend2=bin_length; %bin length
            if fin<bend
               if fin>bstart
                  bend2 = bend2-((bend-fin)+1); %bin length (bend2) is shorter by the difference between fin and bend
                  bin_length=bend2;
               else
                  bend2=nan();
               end
            
            end
            
  %fprintf('mod start: %d, end: %d, first cop: %d\n',1,bend2,fin)

            if ~isnan(bend2)
                total_interactions = nnz(foc_map_cur(1:bend2)>0);

            %ii_firstMating(t,i)=nnz(foc_map(bstart:bend,i)>0)/bend;
                ii(t,i)=total_interactions/bin_length;
                ii_m(t,i)=nnz(sex(foc_map_cur(foc_map_cur(1:bend2)>0))==1)/bin_length;
                ii_f(t,i)=nnz(sex(foc_map_cur(foc_map_cur(1:bend2)>0))==0)/bin_length;
            else
                ii(t,i)=bend2;
                ii_m(t,i)=bend2;
                ii_f(t,i)=bend2;
            end
        
        end
    end
end

function [ai,ai_m,ai_f] = attractionIndex(ai_map,assay_length,numflies,sex,fps)
%get attraction indeces for all individuals in a provided map
ai=zeros(5,numflies);
ai_m=ai;
ai_f=ai;
    for j = 1:numflies
   
        ai(5,j)=nnz(ai_map(:,j)>0)/assay_length;
        ai_m(5,j)= nnz(sex(ai_map(ai_map(:,j)>0,j))==1)/assay_length;
        ai_f(5,j)= nnz(sex(ai_map(ai_map(:,j)>0,j))==0)/assay_length;
        
        %bins=[1 fps*60*10; (fps*60*10)+1 fps*60*20; (fps*60*20)+1 fps*60*30;(fps*60*30)+1 fps*60*40];
        %bin_length=fps*60*10;

        quart=round(assay_length/4);
        bins=[1 quart; quart+1 quart*2; (quart*2)+1 quart*3; (quart*3)+1 assay_length];
        %bin_length=fps*60*10;
        bin_length=quart;
        
        for t=1:4
            bstart=bins(t,1);
            bend=bins(t,2);
            
            ai_map_cur=ai_map(bstart:bend,j);            
            bin_length=quart; %bin length

            ai(t,j)=nnz(ai_map_cur>0)/bin_length;
            ai_m(t,j)= nnz(sex(ai_map_cur(ai_map_cur>0))==1)/bin_length;
            ai_f(t,j)= nnz(sex(ai_map_cur(ai_map_cur>0))==0)/bin_length;
            
        end
    end
end

function [mm_sw,mf_sw,fm_sw,ff_sw,fa_sw] = actionThroughTime(foc_map,win_size,assay_length,numflies,sex)
%this function uses a sliding window approach to identify interactions
%through time. A first pass identifies the number of individuals of each
%sex that show interaction towards the other sexes within a one second
%window. A second pass then smoothes the curve in 10s windows with a 5s
%jump.

%set win_size to fps to get 1 second windows
 
 steps = 1:win_size:assay_length;
 
 m_m=zeros(length(steps)-1,1);
 m_f=m_m;
 f_m=m_m;
 f_f=m_m;
 f_a=m_m;
 
 for i = 2:length(steps)
     n_males=sum(sex==1);
     n_females=sum(sex==0);
     
     males_to_females=0;
     females_to_males=0;
     males_to_males=0;
     females_to_females=0;
     females_to_all=0;
     
     for j = 1:numflies
         to_male=0;
         to_female=0;
         nans=0;
                 
         uni=unique(foc_map(steps(i)-win_size:steps(i),j));
         
         for each = 1:length(uni)
             %fprintf('yo!\n')
             if isnan(uni(each))
                 %fprintf('yo!\n')
 
                 nans=nans+1;
                 
             else
                 if uni(each)>0
                     if sex(uni(each))==1
                        to_male =1;
                     else 
                         to_female =1;
                     end
                 end
             end
             
         end
         %fprintf('num nans: %d\n',nans);
         if nans == win_size
            %fprintf('yo! ind:%d win: %d - %d\n',j,steps(i)-win_size,steps(i));
             if sex(j)==1
                n_males=n_males-1;
             else
                 n_females=n_females-1;
             end
         end
         
         if sex(j)==1
            males_to_females=males_to_females+to_female;
            males_to_males=males_to_males+to_male;
         else
             females_to_females=females_to_females+to_female;
             females_to_males=females_to_males+to_male;
             females_to_all=females_to_all+to_female+to_male;
         end
 
     end
     
     m_m(i-1)=males_to_males/n_males;
     m_f(i-1)=males_to_females/n_males;
     f_m(i-1)=females_to_males/n_females;
     f_f(i-1)=females_to_females/n_females;
     f_a(i-1)=females_to_all/n_females;
         
     
 end
 
 %sliding window over 1 sec bins
 w_size=10; %window size
 slide=5;
 step=10; %end of current window;
 
 mm_sw=[];
 mf_sw=[];
 fm_sw=[];
 ff_sw=[];
 fa_sw=[];
 i=1;
 while step <length(m_m)
     
     mm_sw(i)=nanmean(m_m((step-w_size)+1:step));
     mf_sw(i)=nanmean(m_f((step-w_size)+1:step));
     fm_sw(i)=nanmean(f_m((step-w_size)+1:step));
     ff_sw(i)=nanmean(f_f((step-w_size)+1:step));
     fa_sw(i)=nanmean(f_a((step-w_size)+1:step));
     
     step=step+slide;
     i=i+1;
     
 end


end

function [f,gap_ext] = blocker(n,gap_closer) %put in matrix to translate to blocks; second argument tells the function to close gaps or not value is frames

    f = {}; %make cell array

    gap_ext = n;

    for fly = 1:size(n,2)
         
       blocks = zeros(1,4);
       start = NaN();
       for row =   2: size(n,1)
  
           if n(row,fly) >0 %is there interaction?
                  if n(row-1,fly) ~= n(row,fly) %is interaction new?
                     
                      if ~isnan(start) %was interaction with other fly before? i.e. is start defined?
                            blocks = [blocks; start row-1 n(row-1,fly) (row-1)-start+1]; %save range of last interaction
                           
                      
                      end
                      
                     start = row; %save start of current interaction
                     
                     if row == size(n,1) %if last one is singleton
                         blocks = [blocks; start row n(row,fly) 1];
                     end
                  else
                      if n(row-1,fly) == n(row,fly) && row == 2
                         start = row-1; 
                      else %if last one is interaction, save it
                          if row == size(n,1) && n(row,fly)>0
                              
                              blocks = [blocks; start row n(row,fly) (row)-start+1]; 
                          end
                      end
                  end
                  
           else
               if ~isnan(start) %if n(row,fly) is 0 but start is defined block needs to be saved
                  blocks = [blocks; start row-1 n(row-1,fly) (row-1)-start+1]; %save range of last interaction
                  start = NaN();
               end
           end
           
           
       end
        
       
       if gap_closer
          bl = 2;
          array_dim = size(blocks,1);
          while bl< array_dim
              if blocks(bl+1,1) - blocks(bl,2)<gap_closer && blocks(bl+1,3) == blocks(bl,3)%check if end of last block and begin of next block is less than gap_closer frames away from each other.
                  gap_ext(blocks(bl,2):blocks(bl+1,2),fly)=blocks(bl,3); %closes gap in map array
                  blocks(bl,2) = blocks(bl+1,2); %extends block length
                  blocks(bl+1,:) = []; %deletes obsolete row that was merged into previous block
                  blocks(bl,4) = blocks(bl,2)-blocks(bl,1)+1; %recalculates length of merged block
                  
                  
                  array_dim = array_dim-1;
              else
                  bl=bl+1;
              end
             
              
          end
       end
       
       f{fly,1} = blocks(2:end,:);
    end

end


function [within_vals, between_vals] = selectOr(C,numAss)
    within_corrs = cell(1, numAss);  % preallocate
    for a = 1:numAss
        idx = (1:10) + (a-1)*10;       % indices for assay a
        block = C(idx, idx);          % 10x10 within-assay block
        upper_vals = block(triu(true(10), 1));  % upper triangle without diagonal
        within_corrs{a} = upper_vals;          % store the 45 values
    end
    
    between_corrs = cell(numAss, numAss);  % initialize
    
    for a1 = 1:numAss
        idx1 = (1:10) + (a1-1)*10;
        for a2 = (a1+1):numAss  % only upper diagonal
            idx2 = (1:10) + (a2-1)*10;
            block = C(idx1, idx2);           % 10x10 between-assay block
            between_corrs{a1, a2} = block(:);  % store as 100-element vector
        end
    end
    
    within_vals = cell2mat(within_corrs);
    between_vals = cell2mat(between_corrs(~cellfun(@isempty, between_corrs)));
end

function [P] = pmove(M)
    
    M=logical(M);
    win = -30*5:30*5;                 % window around time t
    [N, T] = size(M);
    P = nan(N,1);
    
    for i = 1:N
        others = setdiff(1:N, i);
    
        % Build a boolean mask of times when any other is moving within the window
        mask = false(1,T);
        anyOthers = any(M(others,:), 1);        % base indicator
    
        % Dilate the indicator by 'win' (logical OR over shifted versions)
        for w = win
            idx = (1:T) + w;
            idx = idx(idx>=1 & idx<=T);
            mask(idx) = mask(idx) | anyOthers(idx - w);
        end
    
        P(i) = mean(M(i, mask));                % conditional probability in window
    end

end