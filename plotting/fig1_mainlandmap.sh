gmt begin bvalue_dist png
    # --- 1. 全局设置 ---
    gmt set MAP_TICK_LENGTH 0p
    gmt set MAP_FRAME_TYPE plain
    
    # --- 2. 绘制底图框架 ---
    # -JA: 兰伯特方位等面积投影，中心(111.3, 32.5)，宽度10cm
    gmt basemap -R80.1/13/142.5/52+r -JA111.3/32.5/10c -Baf 
    
    # --- 3. 绘制地形 ---
    # -t30: 30%透明度，让地形退后，突显地震点
    gmt grdimage @earth_relief_01m -I+d -Cterra -t30
    
    # --- 4. 绘制海岸线和国界 ---
    gmt coast -Dl -W0.5p,black -N1/0.7p,- -A1000 
 
    # --- 5. 绘制地震点 (核心修改部分) ---
    # 修改点：
    # 1. 去掉了 -F, (因为现在是空格分隔)
    # 2. 去掉了 NR>1 (假设这种格式通常没有标题行，如果有标题行，awk的数值判断也能自动过滤掉非数字文本)
    # 3. 逻辑：如果第11列(震级)大于0，则打印 第8列(经度) 和 第7列(纬度)
    
    awk '$11>0 {print $8, $7}' /Users/chouyuhin/_Harvard/b-value/catalogs/CENCcat-format.txt | \
    gmt plot -Sc0.01c -Gred -t60 -W0p

    # --- 6. (可选) 绘制中国国界高亮 ---
    # gmt plot CN-border-L1.gmt -W0.5p,black

gmt end show