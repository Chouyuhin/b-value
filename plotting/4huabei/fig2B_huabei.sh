gmt begin region2_huabei png
    # --- 1. 设置范围与投影 ---
    # 范围：110/34/124/42
    # 投影：中心设在渤海附近
    gmt set MAP_TICK_LENGTH 0p
    gmt set MAP_FRAME_TYPE plain
    gmt basemap -R110/120/34/42 -JM117/38/10c -Baf -BWSen

    # --- 2. 绘制地形 ---
    gmt grdimage @earth_relief_30s -I+d -Cterra -t30

    # --- 3. 绘制行政界线 ---
    gmt coast -W0.5p,black -N1/0.5p,black -N2/0.2p,gray30 -A500 -Slightblue

    # --- 4. 绘制地震点 ---
    awk '$11>0 && $8>=110 && $8<=124 && $7>=34 && $7<=42 {print $8, $7}' \
    /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat-format.txt | \
    gmt plot -Sc0.02c -Gred -t60 -W0p

    # --- 5. 标注关键城市 ---
    echo "116.40 39.90 Beijing" | gmt text -F+f12p,Helvetica-Bold,black -D0/0.3c
    echo "117.20 39.13 Tianjin" | gmt text -F+f12p,Helvetica,black -D0.2c/0

    # --- 6. 比例尺 ---
    gmt coast -Lg122/35+c38+w200k+l"200 km"+f

gmt end show