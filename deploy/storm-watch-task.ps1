# Register storm-watch to start at logon via Windows Task Scheduler.
#   powershell -ExecutionPolicy Bypass -File deploy\storm-watch-task.ps1
$repo = (Resolve-Path "$PSScriptRoot\..").Path
$py   = (Get-Command python).Source
$env:PYTHONPATH = "$repo\src"
$action  = New-ScheduledTaskAction -Execute $py `
    -Argument "-m atmospheric_data storm-watch start config\storm_watch.yaml" -WorkingDirectory $repo
$trigger = New-ScheduledTaskTrigger -AtLogOn
$settings = New-ScheduledTaskSettingsSet -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 1) `
    -StartWhenAvailable
Register-ScheduledTask -TaskName "met_h2o storm-watch" -Action $action -Trigger $trigger `
    -Settings $settings -Description "NWS alert monitor -> auto real_case ingestion" -Force
Write-Host "Registered 'met_h2o storm-watch' (starts at logon). Start now with: Start-ScheduledTask -TaskName 'met_h2o storm-watch'"
