#!/bin/bash

echo Converting ISPC tutorial $1 to CPP tutorial $2
cp $1 $2
sed -i.backup  's/.isph\"/.h\"/g' $2
sed -i.backup  's/RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING/RTC_INTERSECT1/g' $2
sed -i.backup  's/RTC_INTERSECT_VARYING/RTC_INTERSECT1/g' $2
sed -i.backup  's/print(/printf(/g' $2
sed -i.backup  's/uniform //g' $2
sed -i.backup  's/ uniform//g' $2
sed -i.backup  's/varying //g' $2
sed -i.backup  's/ varying//g' $2
sed -i.backup  's/programIndex/0/g' $2
sed -i.backup  's/extern/extern \"C\"/g' $2
sed -i.backup  's/export/extern \"C\"/g' $2
sed -i.backup  's/launch\[numTilesX\*numTilesY\] renderTile(/launch_renderTile(numTilesX\*numTilesY,/g' $2
sed -i.backup  's/launch\[numPhi+1\] animateSphere(/launch_animateSphere(animateSphere,numPhi+1,/g' $2
sed -i.backup  's/M_PI/float(pi)/g' $2
sed -i.backup  's/\*pi\*/\*float(pi)\*/g' $2
sed -i.backup  's/\*pi\//\*float(pi)\//g' $2
sed -i.backup  's/one_over_pi/float(one_over_pi)/g' $2
sed -i.backup  's/one_over_two_pi/float(one_over_two_pi)/g' $2
sed -i.backup  's/one_over_four_pi/float(one_over_four_pi)/g' $2
sed -i.backup  's/[^_]two_pi/float(two_pi)/g' $2
#sed -i.backup  's/RTC_MATRIX_COLUMN_MAJOR/RTC_MATRIX_COLUMN_MAJOR_ALIGNED16/g' $2
sed -i.backup  's/sync;//g' $2
sed -i.backup  's/make_Vec2f/Vec2f/g' $2
sed -i.backup  's/make_Vec3f/Vec3f/g' $2
sed -i.backup  's/make_Vec3fa/Vec3fa/g' $2
sed -i.backup  's/make_Sample3f/Sample3f/g' $2
sed -i.backup  's/\#if 0 \/\/ FIXME: pointer gather/\#if 1 \/\/ FIXME: pointer gather/g' $2
sed -i.backup  's/foreach (i=0 ... N)/for (size_t i = 0; i<N; i++)/g' $2
sed -i.backup  's/foreach (y = y0 ... y1, x = x0 ... x1)/for (int y = y0; y<y1; y++) for (int x = x0; x<x1; x++)/g' $2
sed -i.backup  's/foreach (phi = 0 ... numPhi+1, theta = 0 ... numTheta)/for (int phi = 0; phi <numPhi+1; phi++) for (int theta = 0; theta<numTheta; theta++)/g' $2
sed -i.backup  's/foreach (theta = 0 ... numTheta)/for (int theta = 0; theta<numTheta; theta++)/g' $2
sed -i.backup  's/task void renderTile(int\* pixels,/void renderTile(int taskIndex, int\* pixels,/g' $2
sed -i.backup  's/task void animateSphere (Vertex\* vertices,/void animateSphere (int taskIndex, Vertex\* vertices,/g' $2
sed -i.backup  's/Vec3f renderPixelStandard(float x, float y, const Vec3f\& vx, const Vec3f\& vy, const Vec3f\& vz, const Vec3f\& p)/Vec3fa renderPixelStandard(float x, float y, const Vec3fa\& vx, const Vec3fa\& vy, const Vec3fa\& vz, const Vec3fa\& p)/g' $2
sed -i.backup  's/RTCIntersectFuncVarying/RTCIntersectFunc/g' $2
sed -i.backup  's/RTCOccludedFuncVarying/RTCOccludedFunc/g' $2
#sed -i.backup  's/\#if 1 \/\/ enables parallel execution/\#if 0/g' $2
sed -i.backup  's/RTCFilterFuncVarying/RTCFilterFunc/g' $2
sed -i.backup  's/Vec3f\([^a]\)/Vec3fa\1/g' $2

sed -i.backup  's/new Vec3fa\[12\]/(Vec3fa\*) alignedMalloc(12\*sizeof(Vec3fa))/g' $2
sed -i.backup  's/delete\[\] colors/alignedFree(colors)/g' $2

sed -i.backup  's/new Vec3fa\[width\*height\]/(Vec3fa\*) alignedMalloc(width\*height\*sizeof(Vec3fa))/g' $2
sed -i.backup  's/delete\[\] g_accu/alignedFree(g_accu)/g' $2

sed -i.backup  's/if (id < 0 || id >= numMaterials) continue;//g' $2
sed -i.backup  's/foreach_unique (id in materialID)//g' $2
sed -i.backup  's/ISPCMaterial\* material = \&materials\[id\];/ISPCMaterial\* material = \&materials\[materialID\];/g' $2
sed -i.backup  's/\#define __device__//g' $2
sed -i.backup  's/__device__//g' $2

sed -i.backup  's/make_Ray/RTCRay/g' $2

sed -i.backup 's/\#define PARALLEL_COMMIT//g' $2
 
