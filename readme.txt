2015/12/30 SLIC超像素分割程序SLICSuperpixelSegmentation.exe：
调用方式为命令行命令：
SLICSuperpixelSegmentation image_name 20 superpixel_number output_folder

2016/3/18 程序内容概述：
主程序Main.m分为以下几个部分：
路径和参数设置；
图像的超像素分割；
超像素的代表值求解（初步用的是超像素块的均值）；
代表值的梯度值（为了符合结构张量形式）；
超像素等级结构张量构建，显著性图生成；
超像素级到像素级的映射；

2016/3/18 Function文件夹内部分函数说明：
gray2bins.m	图像量化函数，输入灰度图像和量化等级，输出量化后的图像
bins2hist.m	求解区域的直方图，输入为图像区域和灰度等级，输出为该区域的直方图
sp2pixel.m	超像素到像素等级的映射，输入为原始图像、超像素数目、超像素索引图、超像素邻接矩阵、和超像素等级的显著图，输出为像素级显著图