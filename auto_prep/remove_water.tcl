proc remove_water { x_center y_center} {
	set sel [atomselect top "not water or same fragment as (water and (x-$x_center)^2+(y-$y_center)^2>250)"]
	$sel writepdb removed_water.pdb
}
