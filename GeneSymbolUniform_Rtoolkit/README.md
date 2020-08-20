# GeneSymbolUniform_Rtoolkit
An R script that helps uniform gene symbols of different single-cell RNA sequencing data

欢迎使用Gene Symbol Uniform R toolkit!

# 使用方法

0.确保dependency的R包可用，括号内为参考版本：Seurat(3.1.4), stringr(1.4.0), dplyr(0.8.5)    
dplyr版本需要低于1.0.1（不含），否则可能会有问题

1.解压后进入该文件夹: cd path_to_this_dir，先解压"GeneSymbolRef_SelectAll_upd0731.csv.zip"（zip是因为github文件上传大小限制），文件夹内其他文件及其名称保持不变  

2.输入下面的命令并执行：

**Rscript RToolkit_GeneSymbolUniform.R [query_obj_path] [output_dir]**

query_obj_path: 要处理的Seurat Object的绝对路径  

output_dir: 输出处理后的matrix以及report的路径  

3.运行过程中会输出一些阶段性的提示  

4.运行结束后，会在output_dir下写入两个csv文件：  

UniformedExpression.csv，统一gene symbol list的表达矩阵，目前是43878个gene，cell数目和输入Seurat Object一致  

ModificationReport.csv，gene symbol修改记录，gene list与输入的Seurat Object一致  


**2020.08.05，v0.3 released**  

加入辅助函数"as_matrix"用于将大规模（一般是细胞数大于4w）的系数矩阵转换为matrix进而转换为data.frame。当前版本中，如果细胞数大于3w则使用as_matrix，否则使用R内置的as.matrix函数。

**2020.08.04， v0.2 released**  

将SeuratObj使用的slot从“counts”调整为“data”，即使用normalized data作为处理对象，方便最终的上传。

**2020.08.02， v0.1 released**  

Gene Symbol Uniform R toolkit实现统一基因表达矩阵中基因名称的统一，基于R实现，更好地支持SeuratObj。同时使用命令行完成整套工作，方便简洁。


# 致谢
感谢Qiuchen Meng和Ziheng Zou在total gene list确定过程中的工作，以及后续的测试工作。感谢Yixin Chen和Minsheng Hao在优化工具中的建议和工作。感谢Sijie Chen全过程的帮助。
