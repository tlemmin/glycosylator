set fp [open "dihe.dat" r]
set dihe [read $fp]
close $fp
[atomselect top all] set beta 0

foreach d $dihe {
    set a [atomselect top "serial $d"]
    $a set beta 1.0
    $a delete
}
