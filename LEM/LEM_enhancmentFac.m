function  ETA = LEM_enhancmentFac(sDose)
   % if two SSBs are generated within 25 bp, we score them as an additional DSB AND
   % increase the nummer of DSB induced by two SSBs.
   % Any additional strand break within the 25 bp INTERVAL including two
   % SSB or a DSB is discarded 
   
   alpha_SSB  = 1250;  % SSBs/cell/Gy
   alpha_DSB  = 30 ;   % DSBs/cell/Gy

   h_range    = 25;      % two SSB within h_range base pairs are considered as DSB
   bpCell     = 6*10^9;  % base pair of the cell
   
   N_SSB      = ceil(alpha_SSB * sDose);
   N_DSB      = ceil(alpha_DSB * sDose);
   
   vSSB       = sort(randi([1 bpCell],N_SSB,1));

   %vDSB       = randi([1 bpCell],N_DSB,1);
   vBucket    = linspace(0,bpCell,bpCell/h_range);
   ixBucket   = discretize(vSSB,vBucket);
   bins       = unique(ixBucket);
   
   %% process all values at once
   vCountsRef = histcounts(ixBucket(:),bins)'; 
   N_SSB2     = sum(vCountsRef > 1);
   ETA        = (N_DSB + N_SSB2) / (N_DSB);
   
   
    %% new approach
   vSSB = randi([1 3*10^9],round(N_SSB/2),2);
   vDSB = randi([1 6*10^9],N_DSB,1);

   vSSB(:,1) = sort(vSSB(:,1));
   vSSB(:,2) = sort(vSSB(:,2));

   vDSB      = sort(vDSB);
 
   ixBucketSSB  = discretize(vSSB,vBucket);
   ixBucketDSB  = discretize(vDSB,vBucket);

   bins1SSB  = unique(ixBucketSSB(:,1));
   bins2SSB  = unique(ixBucketSSB(:,2));
   ix        = ismember(bins1SSB,bins2SSB);
   N_SSB     = sum(ix);

   bins1DSB  = unique(ixBucketDSB(:,1));
  
   ixx       = ismember(bins1DSB,bins1SSB(ix));
   N_DSB_2   = N_DSB - sum(ixx);

   ETA       = (N_DSB_2 + N_SSB) / (N_DSB_2);
%    
   
   
   %%
   
   N_SSBref = 0;
   Num      =  numel(vSSB)/2;
   i        = 1;

   waitbar(0) 
   
    while true
        
        window = [vSSB(i,1) - h_range vSSB(i,1) + h_range];

        j = 1;
        
        while true
            
           if   vSSB(j,2) >= window(1) && vSSB(j,2) <= window(2) && i~=j
                N_SSBref = N_SSBref + 1;
                break;
           elseif  vSSB(j,2) > window(2) ||  vSSB(j,2) < window(1)
               break
           end
           
           j = j + 1;
           
            if j >= Num
                 break;
            end
           
        end
        
        i = i + 1;
          
        if mod(i,1000) == 0
          waitbar(i/Num) 
        end
   
        if i >= Num
            break;
        end
        
    end
   
   ETA = (N_DSB  + N_SSBref) / (N_DSB);
   
   
   
   
   %% process values in batches
   NumBatch   = ceil(numel(bins)/100);
   Flag       = true;
   vCounts    = zeros(numel(bins)-1,1);
   Cnt        = 0;
   
   while Flag  
      
      if ~issorted(bins(Cnt+1:Cnt+NumBatch))
          bins(Cnt+1:Cnt+NumBatch) = sort(bins(Cnt+1:Cnt+NumBatch));
      end
      vCnt = histcounts(ixBucket(:),bins(Cnt+1:Cnt+NumBatch))'; 
      
      vCounts(Cnt+1 : Cnt+NumBatch-1,1) = vCnt;
      
      Cnt = Cnt + NumBatch;
      
      if Cnt >= numel(bins)
          Flag = false;
      elseif Cnt + NumBatch > numel(bins)
          NumBatch = numel(bins) - Cnt; 
      end
      
   end
   
  %% calculate enhancement factor
  N_SSB2 = sum(vCounts > 1);  
  ETA    = (N_DSB + N_SSB2) / (N_DSB);
     
  
  %% some old code
% 
%   % get unique values
%   vSSB_disc   = discretize(vSSB,vBucket);
%   [~,ia,~]        = unique(vSSB_disc);
%    
%   % now obtain duplicate values
%   duplicate_ind   = unique(setdiff(1:size(vSSB_disc,1), ia));
%   SSB2            = vSSB_disc(duplicate_ind);
%   N_SSB2ref       = numel(SSB2);
%     
%   vDSB_disc   = discretize(vDSB,vBucket);
%   [~,ib,~]    = unique(vDSB_disc);
%   vDSB_disc   = vDSB_disc(ib);
%   N_DSB2      = numel(setdiff(vDSB_disc,vSSB_disc));
% 
%   ETA2  = (N_DSB + N_SSB2ref) / (N_DSB);
    
   N_SSB3 = 0;
    
   vSSB = sort(vSSB);
   Num  =  numel(vSSB);
   RunLoop = true;
   i       = 1;
   waitbar(0) 
   
    while RunLoop
        
        window = vSSB(i) + h_range;

          if window - vSSB(i+1) > 0
                N_SSB3 = N_SSB3 + 1;

               while true
                    i = i + 1;
                    if i > Num || vSSB(i) > window
                        break;
                    end
                    
               end
               
          else
               i = i + 1;
          end


        if mod(i,1000) == 0
          waitbar(i/Num) 
        end
   
        if i >= Num
            break;
        end
        
    end
   
   ETA3 = (N_DSB  + N_SSB3)/ (N_DSB);
   
   
   N_SSB4  = 0;
   
   intervall_range = 0;
   boolIntervall   = false;
   waitbar(0) 
    for i = 1:numel(vSSB)
            if intervall_range > vSSB(i) && i ~= 1
                if ~boolIntervall 
                   boolIntervall = true;
                   N_SSB4        = N_SSB4 + 1;
                end 
            else
                intervall_range = vSSB(i) + h_range;
                boolIntervall   = false;
            end 
            if mod(i,1000) == 0
              waitbar(i/numel(vSSB)) 
            end
    end
   
 ETA4 = (N_DSB  + N_SSB4)/ (N_DSB);
   
  
   
end









