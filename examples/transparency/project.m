function W = project(V,MV,P,VP)
  V = V * MV' * P';
  % WTF! Perspective division is not in the documentation for gluProject: 
  % http://www.opengl.org/sdk/docs/man2/xhtml/gluProject.xml
  V = bsxfun(@rdivide,V,V(:,4));
  W = [ ...
    VP(1) + (VP(3) * (V(:,1)+1))/2, ...
    VP(2) + (VP(4) * (V(:,2)+1))/2, ...
    (V(:,3)+1)/2];
end
