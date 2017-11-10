ospa_vals= zeros(N,3);
ospa_c= 100;
ospa_p= 1;
for k=1:N
    trueX = [];
    for i=1:TrueTracks
        if(~isnan(x1(k,i))||~isnan(y1(k,i)))
            trueX = [trueX, [x1(k,i);y1(k,i)]];
        end
    end
    estX = [];
    for i=1:numel(Logs)
        if(~isnan(Logs{i}.xV(1:2,k)))
            estX = [estX, Logs{i}.xV(1:2,k)];
        end
    end
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(trueX,estX,ospa_c,ospa_p);
end

figure; ospa= gcf; hold on;
subplot(3,1,1); plot(1:N,ospa_vals(:,1),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
subplot(3,1,2); plot(1:N,ospa_vals(:,2),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
subplot(3,1,3); plot(1:N,ospa_vals(:,3),'k'); grid on; set(gca, 'XLim',[1 N]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
xlabel('Time');