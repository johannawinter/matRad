%% update weights, focus distance, sigma1, sigma2 and LET in APM base data after getting new base data
% test which fields are different
% change machines and fields
for i = 1:255
    A(i) = isequal(machine_APM_gantry.machine.data(i).LET(:),machine_gantry.machine.data(i).LET(:));
end
sum(A)

% update APM (fixed BL) base data
machine_APM.machine.meta.SAD = machine.meta.SAD;
for i = 1:255
    machine_APM.machine.data(i).weight = machine.data(i).weight;
end
for i = 1:255
    machine_APM.machine.data(i).initFocus.dist = machine.data(i).initFocus.dist;
end
for i = 1:255
    machine_APM.machine.data(i).sigma1 = machine.data(i).sigma1;
end
for i = 1:255
    machine_APM.machine.data(i).sigma2 = machine.data(i).sigma2;
end
for i = 1:255
    machine_APM.machine.data(i).LET = machine.data(i).LET;
end
machine_APM.machine.meta.updated_on = '05-Dec-2017';
machine_APM.machine.meta.updated_by = 'winterjo';
machine = machine_APM.machine;
save('C:\Matlab\matrad\protons_HIT_APM','machine')

% update APM gantry base data
machine_APM_gantry.machine.meta.SAD = machine_gantry.machine.meta.SAD;
machine_APM_gantry.machine.meta.BAMStoIsoDist = machine_gantry.machine.meta.BAMStoIsoDist;
for i = 1:255
    machine_APM_gantry.machine.data(i).weight = machine_gantry.machine.data(i).weight;
end
for i = 1:255
    machine_APM_gantry.machine.data(i).initFocus.dist = machine_gantry.machine.data(i).initFocus.dist;
end
for i = 1:255
    machine_APM_gantry.machine.data(i).initFocus.sigma = machine_gantry.machine.data(i).initFocus.sigma;
end
for i = 1:255
    machine_APM_gantry.machine.data(i).sigma1 = machine_gantry.machine.data(i).sigma1;
end
for i = 1:255
    machine_APM_gantry.machine.data(i).sigma2 = machine_gantry.machine.data(i).sigma2;
end
for i = 1:255
    machine_APM_gantry.machine.data(i).LET = machine_gantry.machine.data(i).LET;
end
machine_APM_gantry.machine.meta.updated_on = '05-Dec-2017';
machine_APM_gantry.machine.meta.updated_by = 'winterjo';
machine = machine_APM_gantry.machine;
save('C:\Matlab\matrad\protons_HIT_APMgantry','machine')


%% set all weights >1 to 1
for i = 1:255
    largerOne(i) = any(machine.data(i).weight>1);
    if largerOne(i)
        largerOnePosition(i,1) = find(machine.data(i).weight>1,1);
        largerOnePosition(i,2) = find(machine.data(i).weight>1,1,'last');
    end
    smallerZero(i) = any(machine.data(i).weight<=0);
end

sum(smallerZero)

for i = 1:255
    machine.data(i).weight(machine.data(i).weight>1) = 1;
end

