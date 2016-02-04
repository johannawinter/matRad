% generates a red-blue colormap for difference plots
function costumMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,plane)
    if plane == 1
         img = 100*(squeeze(cube1(slice,:,:))-squeeze(cube2(slice,:,:)))./max(max(cube1(slice,:,:)));
    elseif plane == 2
            img = 100*(squeeze(cube1(:,slice,:))-squeeze(cube2(:,slice,:)))./max(max(cube1(:,slice,:)));
    elseif plane ==3
         img = 100*(cube1(:,:,slice)-cube2(:,:,slice))./max(max(cube1(:,:,slice)));
    end
 
    imgMin = min(img(:));
    imgMax = max(img(:));

    imgRange = linspace(imgMin,imgMax,62);
    [~,idx]  = min(abs(imgRange));
    idx2 = 62-idx;

    a = linspace(0,1,idx);
    b = linspace(1,0,idx2);
    d1 = ones(1,idx);
    d2 = ones(1,idx2);

    blueRow  = [d1 b];
    greenRow = [a b];
    redRow   = [a d2];
    costumMap = [blueRow; greenRow; redRow]'; 
end


