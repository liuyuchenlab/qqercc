# qqercc 
###### 对添加ERCC spike-in的RNA-SEQ数据进行标准化
##### 安装

```
install.packages("devtools")
library(devtools)
```

###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

```
devtools::install_github('liuyuchenlab/qqercc')  
library(qqercc)  
```
###### 数据需要包含ERCC-基因，格式如下

![image](https://github.com/user-attachments/assets/e63b865a-7cc6-441a-be71-c73ad5575ad0)


###### 只可以输入XLSX文件：  
###### 运行函数会自动创建ercc文件夹，并保存ERCC normalized文件，默认使用deseq2方法
###### control_group指的是根据某一组进行relative计算得到总体reads数的比较结果
###### 图片自动保存到ercc文件夹中

```
result <- qqercc("ercc.xlsx", control_group = "control")  
```
##### normalized total reads plot  

![image](https://github.com/user-attachments/assets/4e754b12-9064-4a58-bfb2-4ab76a14ba27)

##### 可以添加自定义颜色和统计比较

```
#定义比较组
comparisons_list <- list(
  c("DMSO-2-cell", "PlaB-2-cell"))
# 定义颜色，与比较组对应
colors <- c("#FF5733", "#33FF57")  # 每个颜色对应comparisons_list中的一组比较
#运行函数
result <- qqercc("plab_ercc.xlsx", control_group = "DMSO-2-cell", 
                 comparisons = comparisons_list, method = "deseq2", 
                 format = "dot", 
                 color = colors
                )
```
![image](https://github.com/user-attachments/assets/7241afa4-73af-4550-944f-343e8919dab8)


##### 可以通过format选择其他可视化形式：点图

```
result <- qqercc("ercc.xlsx", control_group = "control", 
                 comparisons = comparisons_list,
                 format = "dot"
)
```

![image](https://github.com/user-attachments/assets/9910f215-20d4-4eb0-bedd-6053516420c0)

###### 或者箱线图

```
result <- qqercc("ercc.xlsx", control_group = "control", 
                 comparisons = comparisons_list,
                 format = "box"
)
```

![image](https://github.com/user-attachments/assets/e5e39733-85fb-4dbd-a474-7bdda9f61eba)


##### ERCC文件夹

![image](https://github.com/user-attachments/assets/eb67441a-397c-4820-98ff-32f3832acfca)

### 快去试试吧！






