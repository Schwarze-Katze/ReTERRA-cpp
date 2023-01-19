@echo off
pushd %cd%
cd /d "%~dp0"
if not exist build (
    mkdir build
) else (
    @del /S /Q /F build
)
pushd build
cmake .. -DCMAKE_CXX_COMPILER:FILEPATH=C:/msys64/mingw64/bin/g++.exe -DCMAKE_C_COMPILER:FILEPATH=C:/msys64/mingw64/bin/gcc.exe -DCMAKE_EXPORT_COMPILE_COMMANDS=on -G "MinGW Makefiles"
make
popd
popd