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
make
popd
popd