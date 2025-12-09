gmt begin bvalue_dist eps
    # --- 1. 全局设置 ---
    gmt set MAP_FRAME_TYPE plain
    gmt set MAP_TICK_LENGTH 2p
    gmt set FONT_ANNOT_PRIMARY 8p
    gmt set FONT_LABEL 8p
    
    # --- 2. 准备数据 (假设文件有标题行，需要跳过第1行) ---
    # 提取格式: 经度 纬度 震级 震级(用于控制颜色)
    # 震级 > 0 的才绘制
    awk -F, 'NR>1 && $11>0 {print $8, $7, $11, $11}' /Users/chouyuhin/_Harvard/b-value/catalogs/standard_catalog_QJ_alltime.csv > temp_quakes.txt

    # --- 3. 制作地震的色标 (CPT) ---
    # 使用 inferno 色标，从 3级到 8级，颜色从黄到黑/紫，视觉反差大
    gmt makecpt -Cinferno -T3/8/1 -Z -H > quake.cpt

    # --- 4. 绘制底图与地形 ---
    # -R 范围建议稍微大一点以包住整个中国
    # -JA 投影中心设置在从中国几何中心附近 (105E, 35N)
    gmt basemap -R70/140/15/55 -JA105/35/16c -Baf 
    
    # 绘制地形：使用 gray 灰度图作为底色，让彩色的地震点更明显
    # -I+d 增加光照效果
    gmt grdimage @earth_relief_02m -I+d -Cgray -t40

    # 绘制海岸线和国界
    # -N1 绘制国界，灰色
    gmt coast -W0.2p,gray30 -N1/0.5p,gray40 -A1000 

    # --- 5. 绘制地震点 ---
    # 关键参数解释：
    # -Sc0.05c+s0.02 : 基础圆圈大小 0.05cm，+s0.02 表示根据第3列(震级)乘以0.02来缩放大小
    #                  或者直接用 -Sc 接受第3列作为直径。
    #                  这里建议：使用 -SE (Simulate) 或者手动缩放。
    #                  为了简单，我们用 awk 算好直径放在第3列，或者用 -Sc 指定统一单位
    
    # 修正方案：让震级决定大小。
    # 我们可以用 gmt plot 的 -Scc 选项，这需要输入数据为: 经度 纬度 大小(cm) 颜色值
    # 下面这行 awk 把震级转换成圆圈直径: 直径 = (震级^2) / 40 (这是一个经验公式，让大震级指数级变大)
    awk -F, 'NR>1 && $11>0 {print $8, $7, ($11*$11)/50, $11}' /Users/chouyuhin/_Harvard/b-value/catalogs/standard_catalog_QJ_alltime.csv | \
    gmt plot -Sc -Cred -t40 -W0.1p,black
    
    # 解释：
    # -Sc : 读取第3列作为圆的直径(cm)
    # -Cquake.cpt : 根据第4列的数据在 CPT 中查找颜色
    # -t40 : 40% 透明度，处理重叠问题
    # -W0.1p,black : 给每个圆加一个极细的黑边，增加对比度

    # --- 6. 绘制中国国界高亮 (可选) ---
    # gmt plot CN-border-L1.gmt -W0.8p,black

    # --- 7. 添加图例和色标 ---
    # 绘制色标
    gmt colorbar -Cquake.cpt -DjBR+w4c/0.3c+o0.5c/0.5c -Bxa1+l"Magnitude" -G3/8
    
    # 整理临时文件
    rm temp_quakes.txt quake.cpt

gmt end show