@echo off
pushd %cd%
cd /d "%~dp0"
if not exist build (
    mkdir build
) else (
    @del /S /Q /F build
)
pushd build
cmake .. -G "MinGW Makefiles" 
rem -DCMAKE_CXX_COMPILER:FILEPATH=C:/msys64/mingw64/bin/g++.exe -DCMAKE_C_COMPILER:FILEPATH=C:/msys64/mingw64/bin/gcc.exe -DCMAKE_EXPORT_COMPILE_COMMANDS=on
make
popd
popd