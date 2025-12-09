REM  !/usr/bin/env -S bash -e
  SET export GMT_SESSION_NAME=$$
 REM SET MAP_FRAME_TYPE=plain
  SET MAP_FRAME_TYPE=fancy
  SET MAP_FRAME_WIDTH=2p
  SET FORMAT_GEO_MAP=ddd.F
  SET R=117.0/122.0/36.5/41.5
  SET R=113.0/123.5/34/42.5
  SET R=117.0/121.9/36.5/40.0
  SET cpt=terra
 
 gmt begin fig1a_geo png,pdf
    gmt basemap -JM15c -R%R% -BWeSN -Bxa2  -Bya2
    gmt grdimage @earth_relief_15s.grd -R%R%   -C%cpt%   -E300   -I+d  -Q

  rem  ---Border map----
  rem  gmt coast -A4 -Df E+L -I0/white -W0.2p,black
  REM     gmt coast -A+r -Di -W0.2p,black -Swhite -R%R%

rem ---------------Tectonic line------------------ 
     gmt plot -W1.2p,blue,-. Basin_M.dat
rem     gmt plot -W1.2p,-, F_BH.dat

REM ---------------Arrow------------------ 
rem    echo 120.7  39.7  235    2 |  gmt plot -Sv0.5c+bAl  -W1.5p  -Gred
rem     echo 120.6  39.3    55    2  |  gmt plot -Sv0.5c+bAl  -W1.5p  -Gred

REM ---------------AB Line------------------ 
echo 118.18   39.6305 > Line.d
echo 121.0   37.0 >> Line.d
gmt plot -W1.8p,red Line.d

rem ------ sedimentsry layer thickness  -----
gmt contour Hc.dat -W0.25p -C2.5,3.5,4.5,5.5,6.5,7.5  -A2.5,3.5,4.5,5.5,6.5,7.5,+f12p -Gd5i

rem  ---heat points map----
    gmt plot Q_LT60.dat -Sc0.20c -W0.8p -Gwhite
    gmt plot Q_LT60.dat -Sc0.10c -Gwhite

    gmt plot Q_60_90.dat -Sc0.25c -W0.3p -Gwhite
    gmt plot Q_60_90.dat -Sc0.15c -Gred

    gmt plot Q_GT90.dat  -Sc0.32c -W0.3p -Gwhite
    gmt plot Q_GT90.dat  -Sc0.22c  -G255/0/255

rem ---------------City name------------------ 
    gmt plot city1.dat -Sc0.35c -W0.8p,black -Gwhite
    gmt plot city1.dat -Sc0.23c -Gblack

    gmt plot city2.dat -Sc0.25c -W0.8p,black -Gwhite
    gmt plot city2.dat -Sc0.13c -Gblack

  rem  ---- Index cover box_Q0 ----
     echo 120.7  39.2 > box.d
     echo 121.9  39.2 >> box.d
     echo 121.9  40.0 >> box.d
     echo 120.7  40.0 >> box.d    
     gmt plot -Gwhite -L -W0.01p,white box.d

REM ----index/Q points---- 
     echo 120.9   39.7 |gmt plot -Sc0.32c -W0.3p -Gwhite 
     echo 120.9   39.7 |gmt plot -Sc0.22c -G255/0/255 

     echo 120.9   39.5 |gmt plot -Sc0.25c -W0.1p -Gwhite 
     echo 120.9   39.5 |gmt plot -Sc0.15c -Gred

     echo 120.9   39.3 |gmt plot -Sc0.20c -W0.8p -Gwhite 
     echo 120.9   39.3 |gmt plot -Sc0.10c -Gwhite

rem  ---- index cover  ----
     echo 117.0  36.5 > box.d
     echo 118.2  36.5 >> box.d
     echo 118.2  37.5 >> box.d
     echo 117.0  37.5 >> box.d    
     gmt plot -Gwhite -L -W0.01p,white box.d

REM ----index/Sediment thickness line---- 
echo 117.2     37.2 > Line.d
echo 117.5     37.2 >> Line.d
gmt plot -W0.25p  Line.d
echo 117.7     37.2 > Line.d
echo 118.0     37.2 >> Line.d
gmt plot -W0.25p  Line.d

REM ----index/Basin boundry---- 
echo 117.2     36.9 > Line.d
echo 118.0     36.9 >> Line.d
gmt plot -W1.2p,blue,-.  Line.d

REM ---- index/Fualt ---- 
echo 117.2     36.6 > Line.d
echo 118.0     36.6 >> Line.d
rem gmt plot -W1.2p,-  Line.d
gmt plot -W1.8p,red Line.d

gmt end show
 
 pause
