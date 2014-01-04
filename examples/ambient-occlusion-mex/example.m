[V,F] = load_mesh('../shared/cheburashka.off');
N = per_vertex_normals(V,F);
S = ambient_occlusion(V,F,V,N,1000);
L = cotmatrix(V,F);
M = massmatrix(V,F,'voronoi');
[EV,~] = eigs(-L,M,15,'sm');
Z = EV(:,7);
qZ = round(255*(Z - min(Z))/(max(Z)-min(Z)));
C = ind2rgb(qZ,jet(256));
%Z = connected_components(F);
%C = ind2rgb(Z(:),jet(max(Z)));

nsp = 3;
fs = 20;
subplot(nsp,1,1);
set(tsurf(F,V),'CData',Z,'EdgeColor','none',fphong);
title('Pseudo-color  ','FontSize',fs);
axis equal;
view(2);
colormap(jet(256));

subplot(nsp,1,2);
set(tsurf(F,V),'CData',permute(1-[S S S],[1 3 2]),'EdgeColor','none',fphong);
title('Inverted ambient occlusion  ','FontSize',fs);
axis equal;
view(2);

subplot(nsp,1,3);
set(tsurf(F,V),'CData',bsxfun(@times,1-S,C),'EdgeColor','none',fphong);
title('Product  ','FontSize',fs);
axis equal;
view(2);
