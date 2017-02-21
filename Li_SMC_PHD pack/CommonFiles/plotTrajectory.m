   
figure
if ~isempty(X1)
    [a,b]=cart2pol(X1(1,:),X1(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'g.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X2)
    [a,b]=cart2pol(X2(1,:),X2(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'g.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end

if ~isempty(X3)
    [a,b]=cart2pol(X3(1,:),X3(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'r.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X4)
    [a,b]=cart2pol(X4(1,:),X4(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'r.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X5)
    [a,b]=cart2pol(X5(1,:),X5(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'m.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X6)
    [a,b]=cart2pol(X6(1,:),X6(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'m.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X7)
    [a,b]=cart2pol(X7(1,:),X7(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'m.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X8)
    [a,b]=cart2pol(X8(1,:),X8(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'b.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X9)
    [a,b]=cart2pol(X9(1,:),X9(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'b.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
if ~isempty(X10)
    [a,b]=cart2pol(X10(1,:),X10(3,:));
    polar(a(1),b(1),'k^'),hold on
        polar(a,b,'b.-'),hold on
    polar(a(end),b(end),'ks'),hold on
end
ylim([0 2300])
set(gca,'FontSize',14,'FontName','Times new Roman'); set(gcf,'Color','White'); 