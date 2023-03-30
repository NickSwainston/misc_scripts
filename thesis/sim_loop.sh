for i in 1257010784,5:4:48,1.0,P2C 1151680472,18:20,0.3,P1 1274683736,7:00,0.07,P2E; do
    IFS=","
    set -- $i
    obsid=$1
    ra=$2
    size=$3
    phase=$4
    echo ""
    echo $phase:  obs: $obsid  ra: $ra  size: $size
    for freq in 80 155 300; do
        echo $freq MHz
        mwa_tied_array_beam_psf -m ${obsid}_metafits_ppds.fits -r $ra -d -26:42:0 -f $freq -x $size -y $3 -o ${obsid}_${freq}MHz.txt
        python fwhm_from_vcsbeam_psf.py --psf ${obsid}_${freq}MHz.txt --output ${obsid}_${freq}
    done
done