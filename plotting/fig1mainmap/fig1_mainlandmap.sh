#!/usr/bin/env bash
# Fig. 1(a): epicentres of the converted CENC catalog (1970 to August 2023), the three analysis
# regions B-D, and the largest historical earthquakes in and near them.
# Run from this folder in the gmtenv conda environment:  bash fig1_mainlandmap.sh  ->  Fig1a_map.pdf
# Events are selected by date and magnitude, not by line number (the catalog file has changed).
# The 2008-10-21 "ML 8.5" entry (34.758N, 81.498E) is a typing error in the CENC catalog: no global
# agency lists an event there and no aftershocks follow; it is excluded.
set -e
CAT=../../../catalogs/CENCcat_Mw_final.txt
TMP=$(mktemp -d)

# epicentres, without the false 2008-10-21 entry
awk '!($1==2008 && $2==10 && $3==21 && $11>=8) && $11>0 {print $8, $7}' $CAT > $TMP/epi.txt
# largest historical earthquakes (Tangshan 1976, Kunlun 2001, Wenchuan 2008)
awk '($1==1976&&$2==7&&$3==28&&$11>=7.5) || ($1==2001&&$2==11&&$3==14&&$11>=7.5) ||
     ($1==2008&&$2==5&&$3==12&&$11>=7.5) {printf "%s %s %s %04d-%02d-%02d\n", $8, $7, $11, $1, $2, $3}' $CAT > $TMP/stars.txt
test $(wc -l < $TMP/stars.txt) -eq 3

# epicentres as a raster layer: number of events per 0.05-degree cell, drawn where > 0
gmt xyz2grd $TMP/epi.txt -R60/160/0/60 -I0.05 -An -G$TMP/epi.nc
echo "0 125/0/0 1e9 125/0/0" > $TMP/epi.cpt

gmt begin Fig1a_map pdf
    gmt set MAP_TICK_LENGTH 0p MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 10p
    gmt basemap -R80.1/13/142.5/52+r -JA111.3/32.5/16c -Baf
    gmt grdimage @earth_relief_15s -I+d -Cterra -t30 -E300
    gmt coast -Dl -W0.5p,black -N1/0.7p,- -A1000
    gmt grdimage $TMP/epi.nc -C$TMP/epi.cpt -Q -E600 -t10

    # analysis regions (lon/lat rectangles; drawn along parallels and meridians)
    for box in "73 90 35 45" "97 107 22 34" "110 117 34 42"; do
        set -- $box
        printf "%s %s\n%s %s\n%s %s\n%s %s\n%s %s\n" $1 $3 $2 $3 $2 $4 $1 $4 $1 $3 | gmt plot -Ap -W1.2p,white
    done
    printf "81.5 40 B\n102 28 C\n113.5 38 D\n" | gmt text -F+f22p,Helvetica,white+jCM

    # stars and labels (magnitude and date)
    awk '{print $1, $2}' $TMP/stars.txt | gmt plot -Sa0.5c -Gyellow -W0.5p,black
    awk '{printf "> %s %s 11p 2.1c l\nMw %s\n%s\n", $1, $2, $3, $4}' $TMP/stars.txt |
        gmt text -M -F+f9p,Helvetica-Bold+jLM -D0.35c/0 -Gwhite@30 -C2p/2p
gmt end
rm -rf $TMP
