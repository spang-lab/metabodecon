$currentScriptDir = Split-Path -Parent -Path $MyInvocation.MyCommand.Definition
$source = "$currentScriptDir/example_datasets"
$destination = "$currentScriptDir/example_datasets.zip"
Write-Host "Starting the zipping process for $source..."
if(Test-Path $destination) {
    Write-Host "$destination already exists. Removing..."
    Remove-Item $destination
}
Compress-Archive -Path $source -DestinationPath $destination
Write-Host "$destination successfully created."
