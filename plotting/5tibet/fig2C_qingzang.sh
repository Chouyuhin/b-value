gmt begin region3_tibet svg
    # --- 1. 设置范围与投影 ---
    # 范围：72/100/32/44 (覆盖紫色框区域)
    # 投影：墨卡托投影 (-JM)，宽度15cm (因为是宽条带，设宽一点更好看)
    gmt set MAP_TICK_LENGTH 0p
    gmt set MAP_FRAME_TYPE plain
    gmt basemap -R73/90/35/45 -JM15c -Baf -BWSen

    # --- 2. 绘制地形 ---
    # 使用浮雕效果
    gmt grdimage @earth_relief_01m -I+d -Cterra -t30

    # --- 3. 绘制行政界线 ---
    # -N1: 国界 (粗), -N2: 省界 (细)
    gmt coast -W0.5p,black -N1/0.8p,black -N2/0.3p,gray30 -A1000 

    # --- 4. 绘制地震点 (自动筛选该范围内的数据) ---
    awk '$11>0 && $8>=73 && $8<=100 && $7>=32 && $7<=45 {print $8, $7}' \
    /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat-format.txt | \
    gmt plot -Sc0.02c -Gred -t50 -W0p

    # --- 5. 标注关键地名 (辅助定位) ---
    # 标注一些该区域的大城市或地标，方便确认位置
    echo "75.99 39.47 Kashgar" | gmt text -F+f12p,Helvetica-Bold,black -D0/0.2c
    echo "94.93 36.40 Golmud" | gmt text -F+f12p,Helvetica-Bold,black -D0/0.2c
    echo "80.11 37.11 Hotan" | gmt text -F+f12p,Helvetica-Bold,black -D0/0.2c
    echo "94.66 40.14 Dunhuang" | gmt text -F+f12p,Helvetica-Bold,black -D0/0.2c

    # --- 6. 比例尺 ---
    # 在左下角绘制比例尺
    gmt coast -Lg76/34+c38+w400k+l"400 km"+f
gmt end show