
@echo off
 
 echo. > listecpp.txt

 for %%i in (*.cpp) do echo %%~ni.cpp >> listecpp.txt 
 
 
 

set liste=

for /f "delims=" %%a in ('type listecpp.txt') do (


rem echo %%a

set liste=%liste%%a%




 )

 echo !%liste%!

rem echo. > listecpp.txt