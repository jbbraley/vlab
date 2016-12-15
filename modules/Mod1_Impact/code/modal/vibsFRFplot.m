function fh = vibsFRFplot(Hs,hh,in,out,w)
%% Vibs Example - Construct FRF - Mode Contributions
%
% logmag plot of Hpq with dashed mode contributions
% complex conjugate solution
%
% jdv 06232015; 07232015; 08162015; 11112015

% plot frf
titl = sprintf('H%d%d Magnitude',out,in);
fprintf(strcat(titl,'\n'));

% create figure
fh = figure; 
ah = axes; 
fh.PaperPositionMode = 'auto';
fh.Position = [100 100 1300 600];

% plot frf
plot(ah,w,mag2db(abs(Hs)),'o-k','linewidth',2,'markersize',2);
hold(ah,'all')
plot(ah,w,mag2db(abs(hh)),'--','linewidth',1.5);
hold(ah,'off')

% form legend
ne = size(hh,2);
lg = {['H' num2str(out) num2str(in)]};
for ii = 2:ne+1
    lg{ii} = ['Mode: ' num2str(ii-1)];
end
legend(lg,'location','northeast');

% format
grid(ah,'on'); 
grid(ah,'minor');
xlabel(ah,'Radial Frequency [rad/sec]');
ylabel(ah,'[db]');
set(ah,'fontsize',18,'fontname','times new roman');





