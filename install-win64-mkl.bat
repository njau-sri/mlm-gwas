
rmdir /s /q x64

set REDIST=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Redist\MSVC\14.16.27012

set VCOMPDLL=%REDIST%\x64\Microsoft.VC141.OPENMP\vcomp140.dll

call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat"

MSBuild.exe /p:Configuration=Release /p:Platform=x64

cd x64\Release
mkdir mlm-gwas-win64
copy mlm-gwas.exe mlm-gwas-win64\
copy "%VCOMPDLL%" mlm-gwas-win64\

pause
