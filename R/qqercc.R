
#' qqercc
#'
#' @param ercc_file raw counts
#' @param control_group control group
#' @param comparisons t.test
#' @param method mean,median,deseq2
#' @param format bar,box,dot
#' @param color diy
#'
#' @return ercc_normalized file,relative_to_control_picture
#' @export
#'
#' @examples
#' result <- qqercc("ercc.xlsx", control_group = "control")
qqercc <- function(ercc_file, control_group, comparisons = NULL, method = "deseq2", format = "bar", color = NULL) {

  # 检查并加载必要的包
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Please install it.")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please install it.")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required. Please install it.")
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required. Please install it.")
  }

  # 内部加载必要的包
  suppressMessages({
    suppressWarnings({
      library(openxlsx)
      library(ggpubr)
      library(dplyr)
      library(RColorBrewer)
      library(DESeq2)
      library(ggplot2)
    })
  })

  # 创建 'ercc' 文件夹（如果不存在）
  if (!dir.exists("ercc")) {
    dir.create("ercc")
  }

  # 读取数据
  ercc <- openxlsx::read.xlsx(ercc_file, rowNames = TRUE)

  # 提取ERCC数据（仅包含ERCC基因的所有行）
  ercc_data <- ercc[grep("ERCC-", rownames(ercc)), ]

  # 去掉全零行
  ercc_data <- ercc_data[rowSums(ercc_data != 0) > 0, ]

  # 根据method选择标准化基准
  if (method == "mean") {
    # 计算每个样本的ERCC总和
    ercc_totals <- colSums(ercc_data, na.rm = TRUE)

    # 计算ERCC总和的平均值（作为基准）
    ercc_average <- mean(ercc_totals)

    # 计算标准化因子：每个样本的ERCC总和与平均值的比值
    normalization_factors <- ercc_totals / ercc_average

  } else if (method == "median") {
    # 计算每个样本的ERCC中位数
    ercc_medians <- apply(ercc_data, 2, median)  # 计算每个样本的中位数

    # 计算所有样本的中位数的平均值
    total_median_average <- mean(ercc_medians)

    # 计算标准化因子：每个样本的中位数 / 所有样本中位数的平均值
    normalization_factors <- ercc_medians / total_median_average

  } else if (method == "deseq2") {
    # DESeq2方法
    countData <- as.matrix(ercc_data)  # 确保输入是矩阵
    sampleInfo <- data.frame(row.names = colnames(countData), condition = gsub("-\\d+$", "", colnames(countData)))

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = sampleInfo, design = ~1)

    # 提取ERCC基因的行索引
    ercc_rows <- grep("ERCC-", rownames(countData))

    # 基于ERCC基因子集创建一个单独的DESeqDataSet用于估算大小因子
    dds_ercc <- dds[ercc_rows, ]
    dds_ercc <- DESeq2::estimateSizeFactors(dds_ercc)

    # 获取基于ERCC估算的大小因子
    size_factors_ercc <- DESeq2::sizeFactors(dds_ercc)

    # 应用到原始DESeqDataSet
    DESeq2::sizeFactors(dds) <- size_factors_ercc[colnames(dds)]

    # 使用DESeq2进行标准化
    normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))

  } else {
    stop("Invalid method. Please choose 'mean', 'median', or 'deseq2'.")
  }

  # 如果不是 DESeq2，则使用标准化因子对原始数据进行标准化
  if (method != "deseq2") {
    normalized_counts <- sweep(ercc_data, 2, normalization_factors, FUN = "*", check.margin = FALSE)
  }

  # 保存标准化后的结果
  write.csv(normalized_counts, "ercc/ERCC_normalized_counts.csv")

  # ----------------------------------------------
  # 提取分组和重复编号信息
  coldata <- data.frame(
    condition = gsub("-\\d+$", "", colnames(normalized_counts)),  # 提取分组名
    replicate = factor(gsub("^.+-", "", colnames(normalized_counts))),  # 提取重复编号
    row.names = colnames(normalized_counts)
  )

  # ----------------------------------------------
  # 提取控制组的样本数据
  control_samples <- normalized_counts[, grep(control_group, colnames(normalized_counts))]

  # 计算控制组标准化后的总reads数量
  control_totals <- colSums(control_samples, na.rm = TRUE)

  # 计算控制组的标准化后的总reads数的平均值
  control_average <- mean(control_totals)

  # ----------------------------------------------
  # 使用控制组的平均值作为基准对所有样本进行标准化
  final_normalized_counts <- sweep(normalized_counts, 2, control_average, FUN = "/")

  # 计算每个样本的标准化后总reads数量
  plot_data <- data.frame(
    condition = coldata$condition,
    total_reads = colSums(final_normalized_counts, na.rm = TRUE),  # 总reads
    replicate = coldata$replicate
  )

  # ----------------------------------------------
  # 绘图数据准备
  plot_summary <- dplyr::group_by(plot_data, condition, replicate) %>%
    dplyr::summarize(
      mean_reads = mean(total_reads),
      sd_reads = sd(total_reads),
      .groups = "drop"
    )

  # 计算分组平均值
  plot_average <- dplyr::group_by(plot_data, condition) %>%
    dplyr::summarize(
      mean_reads = mean(total_reads),  # 计算每个分组的平均值
      sd_reads = sd(total_reads),
      .groups = "drop"
    )

  # 打印每个分组的平均值（默认）
  print("Average values for each group:")
  print(plot_average)  # 打印每个分组的平均值

  # 使用 RColorBrewer 生成 Paired 调色板颜色
  color_palette <- if (is.null(color)) {
    RColorBrewer::brewer.pal(12, "Paired")[1:length(unique(plot_summary$condition))]
  } else {
    color  # 如果用户提供了自定义颜色，使用用户提供的颜色
  }

  # 全局设置y轴间隔为0.2
  y_axis_scale <- scale_y_continuous(
    breaks = seq(0,  max(plot_summary$mean_reads) * 1.2, by = 0.2),  # 设置y轴的刻度为0.2的间隔
    limits = c(0, max(plot_summary$mean_reads) * 1.2),  # 确保y轴从0开始
    expand = c(0, 0)  # 不让y轴超出下方的0
  )

  # 设置文字字体和大小
  theme_settings <- ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12, face = "bold"),  # 调整字体大小为12，加粗
    axis.title = ggplot2::element_text(size = 14, face = "bold"), # 调整字体大小为14，加粗
    plot.title = ggplot2::element_text(size = 16, face = "bold")  # 调整字体大小为16，加粗
  )

  # 根据format选择图形类型（箱线图、点图或条形图）
  if (format == "box") {
    p_orig <- ggplot2::ggplot(plot_summary, aes(x = condition, y = mean_reads, fill = condition)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
      ggplot2::theme_classic() +
      scale_fill_manual(values = color_palette) +
      y_axis_scale +  # 使用全局设置的y轴刻度
      theme_settings
  } else if (format == "dot") {
    p_orig <- ggplot2::ggplot(plot_summary, aes(x = condition, y = mean_reads, color = condition)) +
      ggplot2::geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
      ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
      ggplot2::theme_classic() +
      scale_color_manual(values = color_palette) +
      y_axis_scale +  # 使用全局设置的y轴刻度
      theme_settings
    # 添加每个分组的平均值线
    p_orig <- p_orig + geom_errorbar(data=plot_average, aes(y=mean_reads, ymax=mean_reads, ymin=mean_reads),
                                     colour="black",width = 0.2,size = 0.8)
  } else {  # 默认使用条形图
    p_orig <- ggpubr::ggbarplot(
      plot_summary,
      x = "condition",
      y = "mean_reads",
      add = "mean_se",              # 添加均值和标准误
      color = "condition",
      fill = "condition",
      palette = color_palette  # 使用随机颜色
    ) +
      ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
      ggplot2::theme_classic() +
      scale_y_continuous(limits = c(0, max(plot_summary$mean_reads) * 1.2), expand = expansion(mult = c(0, 0.05))) +
      y_axis_scale +  # 使用全局设置的y轴刻度
      theme_settings
  }

  # 保存没有比较的图
  ggplot2::ggsave("ercc/ERCC_Normalization_Plot_No_Comparison.pdf", p_orig, width = 10, height = 10)
  print(p_orig)

  # ----------------------------------------------
  # 如果需要进行分组比较
  p_comp <- NULL
  if (!is.null(comparisons)) {
    if (format == "box") {
      p_comp <- ggplot2::ggplot(plot_summary, aes(x = condition, y = mean_reads, fill = condition)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
        ggpubr::stat_compare_means(
          comparisons = comparisons,  # 用户传入的分组比较列表
          method = "t.test",
          label = "p.signif",
          size = 5,
          label.y = max(plot_summary$mean_reads) * 1.1  # 调整显著性标记线的位置
        ) +
        scale_fill_manual(values = color_palette) +
        ggplot2::theme_classic() +
        y_axis_scale +  # 使用全局设置的y轴刻度
        theme_settings
    } else if (format == "dot") {
      p_comp <- ggplot2::ggplot(plot_summary, aes(x = condition, y = mean_reads, color = condition)) +
        ggplot2::geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
        ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
        ggpubr::stat_compare_means(
          comparisons = comparisons,  # 用户传入的分组比较列表
          method = "t.test",
          label = "p.signif",
          size = 5,
          label.y = max(plot_summary$mean_reads) * 1.1  # 调整显著性标记线的位置
        ) +
        scale_color_manual(values = color_palette) +
        ggplot2::theme_classic() +
        y_axis_scale +  # 使用全局设置的y轴刻度
        theme_settings
      # 添加每个分组的平均值线
      p_comp <- p_comp + geom_errorbar(data=plot_average, aes(y=mean_reads, ymax=mean_reads, ymin=mean_reads),
                                       colour="black",width = 0.2,size = 0.8)
    } else {  # 默认使用条形图
      p_comp <- ggpubr::ggbarplot(
        plot_summary,
        x = "condition",
        y = "mean_reads",
        add = "mean_se",              # 添加均值和标准误
        color = "condition",
        fill = "condition",
        palette = color_palette  # 使用随机颜色
      ) +
        ggplot2::labs(y = "Total Reads", x = "", title = "ERCC Spike-in Normalization") +
        ggpubr::stat_compare_means(
          comparisons = comparisons,  # 用户传入的分组比较列表
          method = "t.test",
          label = "p.signif",
          size = 5,
          label.y = max(plot_summary$mean_reads) * 1.1  # 调整显著性标记线的位置
        ) +
        ggplot2::theme_classic() +
        y_axis_scale +  # 使用全局设置的y轴刻度
        theme_settings
    }

    # 保存有比较的图
    ggplot2::ggsave("ercc/ERCC_Normalization_Plot_With_Comparison.pdf", p_comp, width = 10, height = 10)
    print(p_comp)
  }

  # 返回结果
  return(list(
    normalized_counts = normalized_counts,
    plot_no_comparison = p_orig,
    plot_with_comparison = p_comp
  ))
}
