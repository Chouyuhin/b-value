gmt begin longmenshan_basemap pdf
    

    # --- 1. 设置基础底图框 ---
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 10p,Helvetica,black
    gmt basemap -R102/106/30/33.5 -JM15c -Baf -BWSen

    # ========================================================
    # 核心修改区域开始
    # ========================================================
    gmt makecpt -Cgray -G0.99/1 -T-10000/10000 -D
    # --- 2.1 绘制地形底图 (底层) ---
    gmt grdimage @earth_relief_15s -R102/106/30/33.5 -I+a315+nt0.5 -C 
    # --- 5. 绘制地震点 ---
    EQ_FILE="/Users/chouyuhin/_Harvard/b-value/b_value/experiments/1chuandian/1.1龙门山/longmenshan.csv"
    
    # awk 参数详解:
    # -F,     : 告诉 awk 数据是用"逗号"分隔的 (关键!)
    # NR>1    : 跳过第 1 行 (表头 year,month...)
    # $11>0   : 震级(第11列) > 0
    # print   : 输出 经度($8), 纬度($7), 年份($1 用于颜色映射)
    gmt makecpt -C'#274753,#297270,#299d8f,#8ab07c,#e7c66b,#f3a361,#e66d50' -T1970/2023
    awk -F, 'NR>1 && $11>0 && $8>=102 && $8<=106 && $7>=30 && $7<=33.5 {print $8, $7, $1}' \
    "$EQ_FILE" | sort -k3n | \
    gmt plot -Sc0.03c -C -t20
    gmt colorbar -C -DjLB+w4c/0.3c+o0.5c/0.5c -B10+l"Year" -F+gwhite+p0.5p
    

    # ========================================================
    # 核心修改区域结束
    # ========================================================

    # --- 3. 绘制断层 ---
    F_FILE=/Users/chouyuhin/_Harvard/b-value/b_value/plotting/faults/processed_fault_data.txt
    awk '{ if ($1 != p) { print "> " $1 }; p=$1; print $2, $3 }' "$F_FILE" | \
        gmt plot -W1.0p,firebrick

    # --- 4. 绘制行政界线 ---
    # 去掉了 -Dl，用默认精度即可；-W 画岸线（内陆可能不明显）；-N1国界 -N2省界
    gmt coast -N1/0.5p,black -N2/0.3p,gray30
    
    
    # 5. 【新增】绘制最大地震 (黄色五角星)
    # ========================================================
    # 自动找出震级($11)最大的一行并绘制
    awk -F, 'NR>1 {print $0}' "$EQ_FILE" | sort -t, -k11nr | head -n 1 | \
    awk -F, '{print $8, $7}' | \
    gmt plot -Sa0.5c -Gyellow -W0.5p,black
    # 标注震级文字 (向右偏移 0.5cm)
    # awk -F, 'NR>1 {print $0}' "$EQ_FILE" | sort -t, -k11nr | head -n 1 | \
    # awk -F, '{print $8, $7, "M"$11, "("$1")"}' | \
    # gmt text -F+f11p,Helvetica-Bold,black -D0.5c/0
    # 标注成都
    echo "104.06 30.67" | gmt plot -Ss0.2c -Gblack -W0.5p,black
    # echo "104.06 30.67 Chengdu" | gmt text -F+f12p,Helvetica-Bold,black -D0.3c/0

    
gmt end show