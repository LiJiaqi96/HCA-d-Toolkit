细胞类型自动注释工具 v0.1

**1.功能说明**  
输入参考数据集seuratobj路径、查询数据集seuratobj路径、output路径，完成细胞类型注释，并输出csv文件记录注释结果  

**2.依赖及安装方法**  
注意：最好安装R>=4.0或创建新的anaconda环境，否则SingleR的安装可能会失败。
Seurat >= 3.X.X  
BiocManager >= 3.12  
SingleR >= 1.4.1  
SingleCellExperiment >= 1.12.0  

```bash
install.packages("Seurat")  
install.packages("BiocManager")  
BiocManager::install("SingleR")  
BiocManager::install("SingleCellExperiment")  
```

**3.使用方法**  
在满足依赖包的环境中运行：  

Rscript CellTypeAnno_Rtoolit.R ref_dir query_dir out_dir  
  
其中：  
ref_dir: 参考数据集seuratobj路径  
query_dir: 查询数据集seuratobj路径  
out_dir: output路径  

**4.输出文件**    
query_dir+".csv"，csv中存储细胞名称和自动注释的细胞类型  
  
===========================
不断完善中，欢迎提意见！！  
- 是否需要in-place修改seuratobj，增加新的metadata  
- 完整的数据集作参考数据集可能比较慢，是否增设降采样功能  
- 待参考数据集完善后，考虑将输入argument中ref_dir替换为organ  
- ...

