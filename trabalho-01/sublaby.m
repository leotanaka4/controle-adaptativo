function sublaby(myylabel)
% FUNCTION SUBLABELS - coloca labels genéricos em subplots 2x1.
%                      Útil para usar juntamente com o pacote PSFRAG.

hy = ylabel(myylabel);

opy = get(hy);
ha  = gca;
opa = get(ha);

% Nova posicao para o ylabel
xmin = opa.XLim(1);
xmax = opa.XLim(2);
ymin = opa.YLim(1);
ymax = opa.YLim(2);
py = opy.Position;
py(2) = ymax; 
py(1) = -0.8*(xmax/10); % ajuste
%set(hy,'FontName','Helvetica');
set(hy,'FontSize',12);
set(hy,'Position',py);
set(hy,'Rotation',0);
set(hy,'HorizontalAlignment','left');
set(hy,'VerticalAlignment','middle'); %top cap middle baseline bottom

%---x---