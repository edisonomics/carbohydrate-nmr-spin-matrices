# Download the pinned Spinach release into this repository's ignored lib\ folder.
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$version = if ($env:SPINACH_VERSION) { $env:SPINACH_VERSION } else { "2.10.1" }
$target = Join-Path $repoRoot "lib\Spinach-$version"
$archiveUrl = "https://github.com/IlyaKuprov/Spinach/archive/refs/tags/$version.zip"

if (Test-Path (Join-Path $target "kernel")) {
    Write-Host "Spinach $version is already installed at $target"
    exit 0
}

if (Test-Path $target) {
    throw "Found an incomplete installation at $target. Remove it manually, then run this installer again."
}

$temporaryRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("spinach-" + [guid]::NewGuid())
$archivePath = Join-Path $temporaryRoot "spinach.zip"
$extractedRoot = Join-Path $temporaryRoot "Spinach-$version"

try {
    New-Item -ItemType Directory -Path $temporaryRoot | Out-Null
    Write-Host "Downloading Spinach $version..."
    Invoke-WebRequest -Uri $archiveUrl -OutFile $archivePath
    Expand-Archive -Path $archivePath -DestinationPath $temporaryRoot

    if (-not (Test-Path (Join-Path $extractedRoot "kernel"))) {
        throw "The downloaded archive did not contain the expected Spinach kernel."
    }

    New-Item -ItemType Directory -Force -Path (Join-Path $repoRoot "lib") | Out-Null
    Move-Item -Path $extractedRoot -Destination $target

    foreach ($requiredDirectory in @("kernel", "etc", "experiments", "interfaces")) {
        if (-not (Test-Path (Join-Path $target $requiredDirectory))) {
            throw "Spinach installation is missing $requiredDirectory\."
        }
    }

    Write-Host "Spinach $version installed at $target"
    Write-Host "Next: set SPINACH_ROOT=$target before launching MATLAB."
}
finally {
    if (Test-Path $temporaryRoot) {
        Remove-Item -Recurse -Force $temporaryRoot
    }
}
