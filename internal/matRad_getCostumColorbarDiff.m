% generates a red-blue colormap for difference plots
function costumMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,plane)
    if plane == 1
         img = 100*(squeeze(cube1(slice,:,:))-squeeze(cube2(slice,:,:)))./max(max(cube1(slice,:,:)));
    elseif plane == 2
            img = 100*(squeeze(cube1(:,slice,:))-squeeze(cube2(:,slice,:)))./max(max(cube1(:,slice,:)));
    elseif plane ==3
         img = 100*(cube1(:,:,slice)- cube2(:,:,slice))./max(max(cube1(:,:,slice)));
    end
 
    imgMin = min(img(:));
    imgMax = max(img(:));

    imgRange = sort([linspace(imgMin,imgMax,63)]);
    [~,idx]  = min(abs(imgRange));
    idx2 = 64-idx-1;
    a = flip(imgRange(1:idx)/min(imgRange(1:idx)));
    b = flip(imgRange(idx+1:end)/max(imgRange(idx+1:end)));
  
 
    d1 = ones(1,idx);
    d2 = ones(1,idx2);
    blueRow  = [d1 1 b];
    greenRow = [a  1 b];
    redRow   = [a  1 d2];
    
    costumMap = [blueRow; greenRow; redRow]'; 

end


