# Optional Windows native-open verification. Never modifies an existing user deck.
param([string]$AnalysisRoot = 'analysis_development_v2')
$ErrorActionPreference = 'Stop'
$analysisDirectory = (Resolve-Path -LiteralPath $AnalysisRoot).Path
$presentationPath = Join-Path $analysisDirectory 'results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pptx'
$pdfPath = Join-Path $analysisDirectory 'results/BLA_development_E18_P0_P10_P21_2026-10-06_final.pdf'
$sourceHashBefore = (Get-FileHash -LiteralPath $presentationPath -Algorithm SHA256).Hash
$hadPowerPoint = @(Get-Process POWERPNT -ErrorAction SilentlyContinue).Count -gt 0
$application = $null; $deck = $null; $oldSecurity = $null
try {
    $application = New-Object -ComObject PowerPoint.Application
    $oldSecurity = $application.AutomationSecurity
    $application.AutomationSecurity = 3
    $deck = $application.Presentations.Open($presentationPath, -1, 0, 0)
    if ($deck.Slides.Count -ne 18) { throw 'Unexpected native slide count' }
    $charts = 0; $tables = 0; $nativeValues = @()
    foreach ($slide in $deck.Slides) {
        foreach ($shape in $slide.Shapes) {
            if ($shape.HasChart -eq -1) {
                $charts++
                $chart = $shape.Chart
                $seriesCount = $chart.SeriesCollection().Count
                if ($seriesCount -lt 1) { throw 'Native chart has no series' }
                for ($seriesIndex = 1; $seriesIndex -le $seriesCount; $seriesIndex++) {
                    $series = $chart.SeriesCollection($seriesIndex)
                    $nativeValues += [pscustomobject]@{slide = $slide.SlideIndex; series = $series.Name; values = @($series.Values)}
                }
            }
            if ($shape.HasTable -eq -1) { $tables++ }
        }
    }
    if ($charts -ne 5 -or $tables -ne 6) { throw 'Unexpected native chart/table count' }
    # Export a separate PDF; source PPTX remains read-only and unchanged.
    $deck.SaveAs($pdfPath, 32)
    $sourceUnchanged = (Get-FileHash -LiteralPath $presentationPath -Algorithm SHA256).Hash -eq $sourceHashBefore
    if (!$sourceUnchanged) { throw 'Source PPTX changed during native verification' }
    $valuesReadable = @($nativeValues | Where-Object { $null -eq $_.series -or $null -eq $_.values[0] }).Count -eq 0
    $report = [pscustomobject]@{native_open_verified = $true; native_series_readable = $valuesReadable; slides = $deck.Slides.Count; charts = $charts; tables = $tables;
      source_pptx_unchanged = $sourceUnchanged; powerpoint_version = $application.Version; series = $nativeValues; pdf = $pdfPath}
    $report | ConvertTo-Json -Depth 8 | Out-File -LiteralPath (Join-Path $analysisDirectory 'results/validation/native_powerpoint.json') -Encoding utf8
    Write-Output ($report | Select-Object native_open_verified,slides,charts,tables,powerpoint_version | ConvertTo-Json -Compress)
} finally {
    if ($deck) { $deck.Close(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck) }
    if ($application) {
        if ($null -ne $oldSecurity) { $application.AutomationSecurity = $oldSecurity }
        if (!$hadPowerPoint) { $application.Quit() }
        [void][Runtime.InteropServices.Marshal]::ReleaseComObject($application)
    }
}
