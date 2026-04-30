# Copilot 使用指南（仓库特定）

目的：帮助 AI 编码代理快速理解此代码库的结构、主要工作流与常用命令，便于在不需要人工背景说明的情况下完成开发任务。

概览与大图景
- 这是一个面向地震学的 b-value 计算与制图仓库。主要工作流由两类工件驱动：
  - 数据与分析：`experiments/` 下的 Jupyter Notebook（`bvalue4*.ipynb`）负责读取地震目录、计算统计量并生成中间结果（CSV / 图）。
  - 绘图与发布：`plotting/` 与各 `experiments/*/subregion_mainmap_gmt.sh` 脚本使用 GMT（Generic Mapping Tools，见 `plotting/gmt.conf`）生成高质量地图与论文图。

关键文件（优先阅读）
- 数据样例：`bigEQs.csv`（仓库根目录）。
- 分析入口：[experiments/bvalue4mainland.ipynb](experiments/bvalue4mainland.ipynb)（演示了 `pandas`、`numpy`、`matplotlib` 与 `seismostats` 的用法）。
- 区域分析 notebooks：`experiments/*/bvalue4*.ipynb`（每个子区域一套分析流程）。
- 区域地图脚本：`experiments/*/subregion_mainmap_gmt.sh`（调用 GMT 绘制子区域主图）。
- 论文/出版图脚本：`plotting/*/fig2*.sh`、`plotting/fig1mainmap/fig1_mainlandmap.sh`。
- GMT 配置：`plotting/gmt.conf`（脚本依赖的 GMT 默认设置）。

项目约定与模式（对代理重要）
- Notebook-first：大部分分析逻辑以交互式笔记本为主。修改代码时优先在对应 notebook 中定位相关单元。要重现结果，可使用 `jupyter nbconvert --execute` 无头运行。
- Shell-first 绘图：地图与最终图使用 Bash + GMT 脚本；这些脚本通常假定当前工作目录与相对路径，且依赖 `plotting/gmt.conf`。
- 绝对路径注意：部分 notebook（例如 `experiments/bvalue4mainland.ipynb`）引用了本地绝对路径（如 `/Users/chouyuhin/_Harvard/b-value/catalogs/...`）。在异地运行或 CI 中执行前必须替换或参数化这些路径。
- 文件/目录名字：仓库包含中文目录名（例如 `1.1龙门山`、`1.4小江`），代理应以 UTF-8 安全方式处理路径，并优先使用相对路径。

运行与复现（具体示例）
- 运行某个区域的 GMT 地图（在仓库根目录）：
  ```bash
  bash experiments/1chuandian/1.4小江/subregion_mainmap_gmt.sh
  ```
- 无头运行主 Notebook：
  ```bash
  jupyter nbconvert --execute --inplace experiments/bvalue4mainland.ipynb
  ```
- 生成论文图（示例）：
  ```bash
  bash plotting/3chuandian/fig2A_chuandian.sh
  ```

依赖与环境提示
- 必备：GMT（版本仓库中 `plotting/gmt.conf` 为 6.x 语法，建议 >=6.5）、Python（Jupyter）、常见科学包：`pandas, numpy, matplotlib`。
- 特殊：notebook 中使用 `seismostats`（见 `experiments/bvalue4mainland.ipynb`），在运行前确认该包可用或改为等价实现。

调试建议（面向自动代理）
- 若 GMT 脚本失败：逐步运行脚本内部命令（在脚本同目录下手动执行关键行），确认 `gmt.conf` 生效、输入 CSV 存在以及坐标列顺序正确。
- 若 notebook 出错：先在交互式 Jupyter 中逐单元运行定位失败单元；若是路径问题，搜索全仓库 `grep -R "/Users/chouyuhin"` 并替换为相对路径或参数。

变更与合并注意事项
- 修改分析逻辑时，优先修改对应 notebook 并将关键数据导出为 CSV；不要仅修改 notebook 的交互单元而不产生可复用的脚本/函数。
- 对绘图参数（颜色、投影、边界）改动应同时更新 `plotting/gmt.conf` 或相关脚本内的注释，以保持图形可再现性。

要点总结（快速检查清单）
- 查找数据入口：`bigEQs.csv` 与 `experiments/*/*.csv`。
- 分析逻辑：`experiments/*/bvalue4*.ipynb`。
- 绘图逻辑：`experiments/*/subregion_mainmap_gmt.sh` 与 `plotting/*.sh` + `plotting/gmt.conf`。
- 注意绝对路径和中文目录名；确认 `seismostats` 与 GMT 已安装。

若有遗漏或你希望包含 CI/环境安装脚本、确切依赖列表（requirements.txt），我可以把这些内容补充到本文件中。请告诉我接下来要优先完善的部分。
