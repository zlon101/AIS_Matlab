rem 自动执行MainCSA.exe
@echo off
set exePath=E:\astego\CSA\MainCSA_EXE\redistribution_files_only\MainCSA.exe
taskkill /F /IM MainCSA.exe
timeout /t 60

start %exePath% E:\astego\Images\BOSS_ALL\ 1 10
start %exePath% E:\astego\Images\BOSS_ALL\ 1 10 E:\astego\CSA2\

echo on




rem timeout /t 60
rem start /min cmd.exe /K echo start run exe!
rem start cmd.exe /K %exePath% E:\astego\Images\BOSS_ALL\ 1 2

rem 当cmd窗口数量超过2时关闭全部窗口
rem set n=0
rem for /F "tokens=2" %%p in ('tasklist /NH /FI "imagename eq cmd.exe"') do (set /a n+=1)
rem if %n% gtr 2 (taskkill /F /IM cmd.exe)