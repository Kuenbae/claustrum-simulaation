function [rstore] = Recurent_simulation_function(tau,tauE,tauS,tauI,halo_power,g,gInh,nreps,...
    activation_power,recurrent_connectivity,otherarea_recurrent_connectivity,...
    claustrum2otherarea_connectivity,other2claustrum_connectivity)

for rep = 1:nreps
    N = 1000; sp = recurrent_connectivity; thr = 0.1;
    
    sp2=claustrum2otherarea_connectivity;
    sp3=other2claustrum_connectivity;
    
    sp4=otherarea_recurrent_connectivity;
    
 
    dt      = 0.005;
    
    time    = dt:dt:65;
    tOn(:,1)     = 10+rand(N,1)/2;
    tOff(:,1)    = tOn(:,1)+2*randn(N,1);
    
    inh_tOn= 30;
    inh_tOff= 35;
    
    stimRaw = double(time>tOn(:,1) & time<tOff(:,1));
    stim1   = smoothts(stimRaw,'e',4/dt);
    % stim1=stim1/max(eps+stim1);
    stim1(time>tOn+30)=0;
    stim3 = halo_power*double(time>inh_tOn & time<inh_tOff);
    time    = time - 10;
    
    reroll = 1;
    
    doOnline =1;
    
    
    if(doOnline)
        figure(10);clf;
        h = plot(rand(1,10)');
        xlim(time([1 end]));
    end
    
    if(reroll)
        [w1,w2,w3] = deal(zeros(N,1));
        fr = round(N/2);
        w1(randperm(N,fr/10)) = 1/10 * rand(fr/10,1);
        w3(randperm(N,N*0.4)) =10 * rand(N*0.4,1);
    end
    
    rstore = {};
    
    if(reroll)
        J = rand(N).*(rand(N)<sp)/sqrt(N*sp);
        J2= rand(N).*(rand(N)<sp4)/sqrt(N*sp4);
        if sp2==0
            JJ=zeros(N);
        else
            JJ= rand(N).*(rand(N)<sp2)/sqrt(N*sp2);
        end
        if sp3==0
            JJ2=zeros(N);
        else
            JJ2= rand(N).*(rand(N)<sp3)/sqrt(N*sp3);
        end
    end
    [r,x,p,pS,p2,pS2,r2,x2] = deal(zeros(N,length(stim1)));
    I       = zeros(1,length(stim1));
    I2      = zeros(1,length(stim1));
    IS      = zeros(1,length(stim1));
    [pS(:,1),p2(:,1),pS2(:,1),p(:,1)]=deal(rand(N,1).*0.01);
    x(:,1)=rand(N,1)*0.01;
    x2(:,1)=rand(N,1)*0.01;
    Iin=0;
    Iin2=0;
    
    disp(rep);
    for t = 2:length(time)
        
        p(:,t)   = r(:,t-1) + (p(:,t-1) - r(:,t-1)) * exp(-dt/tauE);
        p2(:,t)   = r2(:,t-1) + (p2(:,t-1) - r2(:,t-1)) * exp(-dt/tauE);
        if(t>1/dt)
            smRates = r(:,t-1).*((sum(r(:,max(t-1/dt,1):t-1),2) - 20)>0);
            smRates2 = r2(:,t-1).*((sum(r2(:,max(t-1/dt,1):t-1),2) - 20)>0);
            pS(:,t)  = smRates + (pS(:,t-1) - smRates) * exp(-dt/tauS);
            pS2(:,t)  = smRates2 + (pS2(:,t-1) - smRates2) * exp(-dt/tauS);
        elseif(stim3(t)==0)
            pS(:,t) = pS(:,t-1);
            pS2(:,t) = pS2(:,t-1);
        else
            pS(:,t)  = r(:,t-1) + (pS(:,t-1) - r(:,t-1)) * exp(-dt/tauS);
            pS2(:,t)  = r2(:,t-1) + (pS2(:,t-1) - r2(:,t-1)) * exp(-dt/tauS);
        end
        
        
        
        % membrane potentials
        
        Iin2     = g/2*J2*pS2(:,t-1) + g/2*J2*p2(:,t-1) - (g*(g+gInh)/1000*I2(t-1))...
            + w1'*stim1(:,t) +  rand(N,1)*0.25 +0.5*activation_power*JJ*pS(:,t-1)+0.5*activation_power*JJ*p(:,t-1);
        x2(:,t)  = Iin2 + (x2(:,t-1) - Iin2) * exp(-dt/tau);
        
        
        
        Iin     = g/2*J*pS(:,t-1) + g/2*J*p(:,t-1) - (g*(g+gInh)/1000*I(t-1))...
            + w1'*stim1(:,t) + rand(N,1)*0.25 +w3* stim3(t)+0.5*activation_power*JJ2*pS2(:,t-1)+0.5*activation_power*JJ2*p2(:,t-1);
        x(:,t)  = Iin + (x(:,t-1) - Iin) * exp(-dt/tau);
        
        
        
        % spikes
        r(:,t)  = (x(:,t)>=thr)*1/dt/100;
        x(r(:,t)~=0,t) = 0;
        
        r2(:,t)  = (x2(:,t)>=thr)*1/dt/100;
        x2(r2(:,t)~=0,t) = 0;
        
        %inhibition
        if otherarea_recurrent_connectivity==0
            I(t) = sum(r(:,t-1)) + (I(t-1) - sum(r(:,t-1))) * exp(-dt/tauI);
        end
        I2(t) = sum(r2(:,t-1)) + (I2(t-1) - sum(r2(:,t-1))) * exp(-dt/tauI);
        
        if(mod(t,200)==0 && doOnline)
            disp([num2str(t) '/' num2str(length(time))]);
            %             set(h,'xdata',time(1:t),'ydata',smoothts(mean(r(:,1:t),1),'b',.1/dt));
            set(h,'xdata',time(1:t),'ydata',mean(r(:,1:t),1));
            drawnow;
        end
        if(time(t)>56 && time(t)<60)
            r(:,t)=0;x(:,t)=0;I(t)=0;Iin=0;pS(:,t)=0;p(:,t)=0;
            %             smRate=zeros(size(smRates));
        end
    end
    
    % compute R1 and R2 from model raster:
    rIter1 = r(:,time > 0  & time<55)~=0;
    rIter2 = r(:,time > 60 & time<(60+55))~=0;
    rstore{rep,1} = rIter1;
    rstore{rep,2} = rIter2;

    
    graph_title1=['recurrent connectivity=' num2str(recurrent_connectivity) ' '...
        ' other2cla connectivity= ' num2str(other2claustrum_connectivity)...
        'other2other= ' num2str(otherarea_recurrent_connectivity)...
        'g= ' num2str(g)];
    graph_title2=['cla2other connectivity=' num2str(claustrum2otherarea_connectivity) ' ' ...
        ' halo power=' num2str(halo_power)...
        'g inh= ' num2str(gInh) ' taus=' num2str(tauS) '  nreps= ' num2str(rep)];
    
    
    
    thr         = 0;
    sig         = smoothts(smoothts(double(rstore{rep,1}),'e',.1/dt),'e',2.5/dt);
    response    = max(sig,[],2);
    sig         = sig(response > thr,:);
    
    [pks,locs] = max(sig,[],2);
    temp = sortrows([locs sig],1);
    temp = temp(:,2:end);
    temp = bsxfun(@minus,temp,median(temp(:,1:110),2));
    temp = bsxfun(@times,temp,1./max(temp,[],2));
    plottime = (1:size(temp,2))*dt;
    temp_image=uint8((temp')'*128+64/5);
    
    figure();clf;
    image(plottime,[],temp_image);
    xlim([0 55]);
    title(graph_title1,graph_title2)
    
    save([graph_title1 graph_title2 '.mat'],'r');
    save('I','I');
end

disp('finish')









