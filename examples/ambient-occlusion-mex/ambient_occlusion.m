  % AMBIENT_OCCLUSION Compute ambient occlusion per given point
  %
  % S = ambient_occlusion(V,F,P,N,num_samples)
  %
  % Inputs:
  %    V  #V by 3 list of mesh vertex positions
  %    F  #F by 3 list of mesh triangle facet indices into V
  %    P  #P by 3 list of origin points
  %    N  #P by 3 list of origin normals
  %    num_samples  number of samples
  % Outputs:
  %    S  #P list of ambient occlusion values between 1 (fully occluded) and 0
  %      (not occluded)
  %
