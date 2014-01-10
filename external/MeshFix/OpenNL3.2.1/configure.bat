echo OFF

rem GENERATOR is the tool you want to compile the sources with. See below for possible targets.
set GENERATOR="Visual Studio 9 2008"

rem set here the path to your installation of CMake.
set CMAKE_PATH="C:\Program Files (x86)\CMake 2.8"


rem Possible values for the generator (-G) are as follow. Ues the one suitable for you
rem
rem 	"Borland Makefiles"
rem 	"MSYS Makefiles"
rem 	"MinGW Makefiles"
rem 	"NMake Makefiles"
rem 	"Unix Makefiles"
rem 	"Visual Studio 6"
rem 	"Visual Studio 7"
rem 	"Visual Studio 7 .NET 2003"
rem 	"Visual Studio 8 2005"
rem 	"Visual Studio 8 2005 Win64"
rem 	"Visual Studio 9 2008"
rem 	"Visual Studio 9 2008 Win64"
rem     "Visual Studio 10"
rem     "Visual Studio 10 Win64"
rem 	"Watcom WMake"
rem 	"CodeBlocks - MinGW Makefiles"
rem     "CodeBlocks - NMake Makefiles"
rem 	"CodeBlocks - Unix Makefiles"
rem 	"Eclipse CDT4 - MinGW Makefiles"
rem 	"Eclipse CDT4 - NMake Makefiles"
rem 	"Eclipse CDT4 - Unix Makefiles"

rem *****************************************************************************

rem Creating directories if they are missing
if NOT EXIST build (
	echo creating build directory...
	md build
)
if NOT EXIST build\Windows (
	echo creating build\Windows directory...
	md build\Windows
)
chdir build\Windows\

echo ***** Creating %GENERATOR% project *****

rem Reg QUERY "HKLM\SOFTWARE\Wow6432Node\Kitware\Cmake 2.6.1" /ve 
	
rem for /f "tokens=2* delims= " %%i in ('reg query "HKLM\SOFTWARE\Wow6432Node\Kitware\Cmake 2.6.1" /ve') do (
rem 	set CMAKE_PATH=%%j
rem )
rem echo CMake found at %CMAKE_PATH%

%CMAKE_PATH%\bin\cmake.exe -G %GENERATOR% ..\..\  
pause
