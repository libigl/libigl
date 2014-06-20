rm -fr buildXcode
mkdir buildXcode
cd buildXcode
cmake -GXcode ../
open buildXcode/*.xcodeproj
