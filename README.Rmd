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



```
data(PCR)
#PCR <- read.csv('test.csv')  
```
###### Gapdh可以换成其他基因，这一步可以得到相对表达量和三个图片：  

###### 图片在plot里，第一张是普通的柱形图，第二张加上了p值，第三张加上了显著性标记 
###### p值和显著性标记的位置是根据分组设置高度的，可以去AI或者PS里平移到合适位置
###### p值计算使用近似T检验，可以使用relative数据去其他知名软件计算并绘图
```
results <- qqpcr(PCR, reference_gene = "Gapdh", control_group = "control")  
```
##### 柱形图  
![image](https://github.com/user-attachments/assets/e88aaf11-bacf-46d7-bea7-db0787323342)

##### p值柱形图  
![image](https://github.com/user-attachments/assets/6498d295-17fe-4fc0-ad6c-c0fdf9a4f439)

##### 显著性标记柱形图  
![image](https://github.com/user-attachments/assets/0bd87dce-9cad-432c-8626-0d473a13eb75)

###### 函数运行时会自动保存行列转换的CT值、相对表达值、统计数据和p值还有一张basic图片到自动创建的文件夹中，文件夹名为当天的日期。
![image](https://github.com/user-attachments/assets/5f8052c7-84cc-4bd9-92d3-ca7012261c1b)

##### 默认只有16个颜色，也可以添加自定义颜色和截断
```
custom_colors <- c("#FF0000", "#00FF00", "#0000FF")  # 用户自定义的颜色

results <- qqpcr(PCR, 'Gapdh', 'control',  breaks = list(c(4, 10), c(14, 20),c(25,30)), custom_colors = custom_colors)

```
##### 示例展示
![plot_zoom_png](https://github.com/user-attachments/assets/79e13b2c-06e6-41a3-a904-d706bf6d915e)

##### 也可以自定义保存
![image](https://github.com/user-attachments/assets/b661ec38-b609-431c-b799-658d3f0577bc)

### 快去试试吧！






