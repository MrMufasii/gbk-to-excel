<#
    Uninstall-RightClick.ps1
    Removes the right-click "Extract annotations to Excel" option added by
    Install-RightClick.ps1. Per-user only; leaves gbk_to_excel.py untouched.
#>

$ErrorActionPreference = "SilentlyContinue"

foreach ($ext in @(".gbk", ".gbff", ".gb", ".genbank")) {
    Remove-Item -Path "HKCU:\Software\Classes\SystemFileAssociations\$ext\shell\GbkToExcel" -Recurse -Force
}
Remove-Item -Path "HKCU:\Software\Classes\Directory\shell\GbkToExcel" -Recurse -Force

$lnk = Join-Path $env:APPDATA "Microsoft\Windows\SendTo\Extract GBK to Excel.lnk"
Remove-Item -Path $lnk -Force

Write-Host "Removed the right-click / Send-to options. (gbk_to_excel.py and gbk_to_excel.cmd were left in place.)"
