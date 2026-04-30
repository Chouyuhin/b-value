gmt begin xiaojiang_basemap pdf

    # --- 1. 设置基础底图框 ---
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 10p,Helvetica,black
    gmt basemap -R101.8/104.4/26/27.7 -JM15c -Baf -BWSen

    # --- 2. 绘制地形底图 ---
    gmt makecpt -Cgray -G0.99/1 -T-10000/10000 -D
    gmt grdimage @earth_relief_01s -R101.8/104.4/26/27.7 -I+a315+nt0.5 -C 

    # --- 3. 绘制地震点 ---
    EQ_FILE="/Users/chouyuhin/_Harvard/b-value/b_value/experiments/1chuandian/1.4小江/xiaojiang.csv"
    gmt makecpt -C'#274753,#297270,#299d8f,#8ab07c,#e7c66b,#f3a361,#e66d50' -T1970/2023
    
    awk -F, 'NR>1 && $11>0 && $8>=101.8 && $8<=104.4 && $7>=26 && $7<=27.7 {print $8, $7, $1}' \
    "$EQ_FILE" | sort -k3n | \
    gmt plot -Sc0.03c -C -t20
    
    gmt colorbar -C -DjLB+w4c/0.3c+o0.5c/0.5c -B10+l"Year" -F+gwhite+p0.5p

    # ========================================================
    # 4. 绘制断层线条 (红色线条)
    # ========================================================
    F_FILE="/Users/chouyuhin/_Harvard/b-value/b_value/plotting/faults/processed_fault_data.txt"
    
    if [ -f "$F_FILE" ]; then
        # 画线逻辑：当名字($1)变化时，插入 GMT 分隔符 ">"
        awk '{ if ($1 != p) { print "> " $1 }; p=$1; print $2, $3 }' "$F_FILE" | \
            gmt plot -W1.0p,firebrick

        # ========================================================
        # 4.1 【简化版】绘制断层名字
        # 逻辑：只要发现名字($1)跟上一行不一样，就在当前点($2, $3) 写上名字
        # ========================================================
        # ========================================================
        # 4.1 【修复版】绘制断层名字
        # ========================================================
        # 逻辑改进：
        # 1. ($2+0) > 90: 强制检查第2列(经度)必须大于90。
        #    这会自动过滤掉表头(Lon=0)、空行(0)或非数字行。
        # 2. $1 != p: 名字发生变化时才标注。
        
    #     awk '{
    #     # 1. 过滤：只允许经度 > 100 的行通过 (防止读到表头变成0,0)
    #     if ( ($2+0) > 100 && $1 != p ) {
            
    #         # 2. 定位：打印这一行的经纬度 -> 这就是"该在的地方"
    #         print $2, $3, $1
    #         p=$1
    #     }
    # }' "$F_FILE" | \
    # # 3. 微调：在上述坐标基础上，向右挪动 0.2cm (-D)
    # gmt text -F+f7p,Helvetica-Bold,firebrick -D0.2c/0 -Gwhite@40 -N
    fi

    # --- 5. 绘制行政界线 ---
    gmt coast -N1/0.5p,black -N2/0.3p,gray30

    # ========================================================
    # 6. 绘制最大地震 (包含年月日)
    # ========================================================
    
    # 6.1 画五角星 (黄色)
    awk -F, 'NR>1 {print $0}' "$EQ_FILE" | sort -t, -k11nr | head -n 3 | \
    awk -F, '{print $8, $7}' | \
    gmt plot -Sa0.5c -Gyellow -W0.5p,black
    
    # 6.2 标注文字 (震级 + 年月日)
    # print内容: 经度, 纬度, "M"+震级, 年"/"月"/"日
    # 结果示例: M6.8 2022/9/5
    awk -F, 'NR>1 {print $0}' "$EQ_FILE" | sort -t, -k11nr | head -n 3 | \
    awk -F, '{print $8, $7, "M"$11, $1"/"$2"/"$3}' | \
    gmt text -F+f11p,Helvetica-Bold,black -D0.5c/0

gmt end show