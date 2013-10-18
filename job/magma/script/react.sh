#!/bin/sh
#
# Name of smirksfile should not contain '>' character !
# 
set im [molfile read [molfile open stdin r hydrogens add] "" all]
set mode [lindex $argv 0]
set argv [lreplace $argv 0 0]
foreach arg $argv {
    set smirksfile [open $arg r]
    set smirkslist [list]
    while {[gets $smirksfile line] != -1} {
        if {$line != "" && ![string match "#*" $line]} {
        	set options [list removeh appendpathname checkaro rebuildaro]
        	if {[lindex $line 3] == "y"} {
        		lappend options preservecoordinates
        		}
            lappend smirkslist [list [concat [lindex $line 0] [lindex $line 2]] 1 forward $options]
            }
        }
    }
# puts stderr $smirkslist
set fhandle [molfile open stdout w hydrogens stripall]
foreach react $im {
    if {$mode == "parallel"} {
        molfile write $fhandle $react
        }
    # puts stderr [ens new "$react" E_SMILES {}]
    set products [ens transform $react $smirkslist "" "singlestep" $mode {
            removeh appendpathname checkaro rebuildaro
            } "" "" "" "" "" 1]
    if {$products != ""} {
        #set all_products [ens create]
        #foreach product $products {
        #    ens merge $all_products $product
        #     }
        #set unique_products [ens weed $all_products {duplicates}]
        foreach product $products {
            foreach molecule [ens split [ens weed $product {duplicates}]] {
                match ss -align besteffort -fuzz 2 $react $molecule
                molfile write $fhandle $molecule
                }
            }
        }
    }
