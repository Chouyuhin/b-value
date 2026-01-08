gmt begin bvalue_dist png
    # --- 1. 全局设置 ---
    gmt set MAP_TICK_LENGTH 0p
    gmt set MAP_FRAME_TYPE plain
    
    # --- 2. 绘制底图框架 ---
    # -JA: 兰伯特方位等面积投影，中心(111.3, 32.5)，宽度10cm
    gmt basemap -R80.1/13/142.5/52+r -JA111.3/32.5/10c -Baf 
    
    # --- 3. 绘制地形 ---
    # -t30: 30%透明度，让地形退后，突显地震点
    gmt grdimage @earth_relief_15s -I+d -Cterra -t30
    
    # --- 4. 绘制海岸线和国界 ---
    gmt coast -Dl -W0.5p,black -N1/0.7p,- -A1000 
 
    # --- 5. 绘制地震点 (核心修改部分) ---
    # 修改点：
    # 1. 去掉了 -F, (因为现在是空格分隔)
    # 2. 去掉了 NR>1 (假设这种格式通常没有标题行，如果有标题行，awk的数值判断也能自动过滤掉非数字文本)
    # 3. 逻辑：如果第11列(震级)大于0，则打印 第8列(经度) 和 第7列(纬度)
    
    awk '$11>0 {print $8, $7}' /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat-format.txt | \
    gmt plot -Sc0.01c -Gred -t60 -W0p

    awk 'NR==456523 || NR==46526 || NR==435990 || NR==296238 {print $8, $7}' \
    /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat_Mw_final.txt | \
    gmt plot -Sa0.4c -Gyellow -W0.5p,black
    # awk 'NR==456523 || NR==46526  || NR==435990 || NR==296238 \
    # {print $8, $7, "M"$11}' \
    # /Users/chouyuhin/_Harvard/b-value/b_value/catalogs/CENCcat_Mw_final.txt | \
    # gmt text -F+f8p,Helvetica-Bold,black -D0.3c/0 -Gwhite@40
    # --- 6. (可选) 绘制中国国界高亮 ---
    # gmt plot CN-border-L1.gmt -W0.5p,black
    # --- 6. 绘制矩形框 (新增部分) ---
    # 设定矩形范围 (例如：经度 102-106，纬度 30-34)
    # -W1.5p,blue  : 蓝色实线边框
    # -L           : 确保闭合
    # 如果你想给框里加个半透明底色，可以加上 -Gwhite@80 (白色80%透明)
    # 龙门山地区
    # echo "102 30
    # 106 30
    # 106 34
    # 102 34
    # 102 30" | gmt plot -W1.5p,blue -L
    # 华北平原地区
    # echo "110 34
    # 120 34
    # 120 42
    # 110 42
    # 110 34" | gmt plot -W1p,white -L
    # # 川滇
    # echo "97 22 107 34" | gmt plot -Sr+s -W1p,white
    # # --- 绘制青藏高原主体范围 ---
    # echo "73 35 90 45" | gmt plot -Sr+s -W1p,white
gmt end show