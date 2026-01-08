gmt begin region1_chuandian png
    # --- 1. 设置范围与投影 ---
    # 范围：基于你之前的设定 97/22/107/34
    # 投影：墨卡托投影 (-JM)，宽度设为 8cm (适合子图)
    # 或者继续用兰伯特 (-JA)，这里为了方便查看局部细节推荐用 -JM
    gmt set MAP_TICK_LENGTH 0p
    gmt set MAP_FRAME_TYPE plain
    gmt basemap -R97/107/22/34 -JM103/28/8c -Baf -BWSen

    # --- 2. 绘制地形 ---
    gmt grdimage @earth_relief_30s -I+d -Cterra -t30

    # --- 3. 绘制行政界线与水系 ---
    # -N1: 国界, -N2: 省界 (川滇地区省界很重要)
    gmt coast -W0.5p,black -N1/0.5p,black -N2/0.2p,gray30 -A1000 

    # --- 4. 绘制地震点 ---
    # 筛选该范围内的点进行绘制
    awk '$11>0 && $8>=97 && $8<=107 && $7>=22 && $7<=34 {print $8, $7}' \
    /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat-format.txt | \
    gmt plot -Sc0.02c -Gred -t60 -W0p

    # --- 5. (可选) 标注关键地名 ---
    echo "104.06 30.67 Chengdu" | gmt text -F+f12p,Helvetica-Bold,black -D0.2c/0

    # --- 6. 比例尺 ---
    gmt coast -Lg98/23+c28+w200k+l"200 km"+f

gmt end show