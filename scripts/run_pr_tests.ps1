param(
    [Parameter(Mandatory=$true)][string]$TestExecutable,
    [Parameter(Mandatory=$true)][string]$VirtualLinacRoot,
    [Parameter(Mandatory=$true)][string]$TemporaryOutputFolder
)
$SUCCESS = 0
& $TestExecutable --test-suite-exclude="System" -- --root $VirtualLinacRoot
$SUCCESS += $LASTEXITCODE

$system_tests = @(
    'UniformPointSource',
    'WaterCylinderAndGantry',
    'CollimatorAndGantry'
)
Remove-Item -Path $TemporaryOutputFolder -Include *.dose,*.phsp
Foreach ($test in $system_tests) {
    & $TestExecutable --test-case="$test" -- --root $VirtualLinacRoot --out $TemporaryOutputFolder
    $SUCCESS += $LASTEXITCODE
    Remove-Item -Path $TemporaryOutputFolder -Include *.dose,*.phsp
}
Write-Output ("Exit code: {0}" -f $SUCCESS)
exit $SUCCESS