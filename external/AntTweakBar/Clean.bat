RMDIR /S /Q src\debug32
RMDIR /S /Q src\debug64
RMDIR /S /Q src\release32
RMDIR /S /Q src\release64
RMDIR /S /Q src\ipch
CD src
DEL *.ncb *.aps *.o *.bak *.user *.sdf
DEL /A:h *.suo
CD ..
RMDIR /S /Q lib\debug
RMDIR /S /Q examples\debug32
RMDIR /S /Q examples\debug64
RMDIR /S /Q examples\tmp
DEL lib\*.exp
CD examples
DEL *.ncb *.aps *.o *.bak *.user *.sdf
DEL /A:h *.suo
DEL /S BuildLog.htm
CD ..

PAUSE


