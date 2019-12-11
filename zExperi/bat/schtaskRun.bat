@echo off
schtasks /create /sc hourly /mo 1 /tn MainCsaTask /tr E:\astego\bat\restartcsa.bat
schtasks /run /tn MainCsaTask
::schtasks /create /sc minute /mo 2 /sd 03/01/2002 /tn "MainCsaTask" /tr E:\astego\bat\restartcsa.bat


@echo on
::schtasks /end /tn "My App"
::schtasks /run /tn MainCsaTask