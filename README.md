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
###### normalized total reads plot  

![image](https://github.com/user-attachments/assets/8a8aac0f-0e22-4315-be0f-546eab93a4de)


###### 可以添加自定义颜色和统计比较

```
#定义比较组
comparisons_list <- list(
  c("control", "treat"))
# 定义颜色，与比较组对应
colors <- c("#CE1212", "#810000")  # 每个颜色对应comparisons_list中的一组比较
#运行函数
result <- qqercc("ercc.xlsx", control_group = "control", 
                 comparisons = comparisons_list, method = "deseq2", 
                 format = "bar", 
                 color = colors
)
```
![image](https://github.com/user-attachments/assets/e4c531e2-463d-4382-b099-9904efba9a38)



###### 可以通过format选择其他可视化形式：点图

```
result <- qqercc("ercc.xlsx", control_group = "control", 
                 comparisons = comparisons_list, method = "deseq2", 
                 format = "dot", 
                 color = colors
)
```

![image](https://github.com/user-attachments/assets/433a7be2-61cb-4d94-923e-350a59d3a5d2)


###### 或者箱线图

```
result <- qqercc("ercc.xlsx", control_group = "control", 
                 comparisons = comparisons_list, method = "deseq2", 
                 format = "box", 
                 color = colors
)
```

![image](https://github.com/user-attachments/assets/14618b6c-01ec-42c0-a9c6-3ff45de79048)



##### ERCC文件夹

![image](https://github.com/user-attachments/assets/eb67441a-397c-4820-98ff-32f3832acfca)

### 快去试试吧！






