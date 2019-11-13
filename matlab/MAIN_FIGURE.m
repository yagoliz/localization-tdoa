figure(1)

load('delay_time_RS.mat')
load('delay_time_US.mat')
load('delay_time_RS_ceck.mat')
load('PPM_final_def.mat')

sgtitle(sprintf('ANALYZE DISTRIBUTION OF D\\_US AND D\\_RSC | 50 TESTS\nEstimate fo [PPM]:  | RTL1  %.3f | RTL2  %.3f | RTL3  %.3f |',PPM_final_def(1),PPM_final_def(2),PPM_final_def(3)),'fontweight','bold');

subplot(3,2,1)
hist(delay_time_US(1,:),-1:3)
grid; ylim([0 50])
title(sprintf('RTL1 and RTL2 - D\\_US'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');

subplot(3,2,2)
hist(delay_time_RS_ceck(1,:),-1:3)
grid; ylim([0 50])
title(sprintf('RTL1 and RTL2 - D\\_RSC'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');

subplot(3,2,3)
hist(delay_time_US(2,:),-1:3)
grid; ylim([0 50])
title(sprintf('RTL1 and RTL3 - D\\_US'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');

subplot(3,2,4)
hist(delay_time_RS_ceck(2,:),-1:3)
grid; ylim([0 50])
title(sprintf('RTL1 and RTL3 - D\\_RSC'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');

subplot(3,2,5)
hist(delay_time_US(3,:),-1:3)
grid; ylim([0 50])
title(sprintf('RTL2 and RTL3 - D\\_US'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');

subplot(3,2,6)
hist(abs(delay_time_RS_ceck(3,:)),-1:3)
grid; ylim([0 50])
title(sprintf('RTL2 and RTL3 - D\\_RSC'))
ylabel('# Experiments','fontweight','bold');
xlabel('Samples','fontweight','bold');
h = findobj(gca,'Type','patch');
