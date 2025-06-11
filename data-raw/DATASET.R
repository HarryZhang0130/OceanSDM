library(devtools)
# 创建license
usethis::use_gpl_license(version = 3, include_future = TRUE)
# 创建函数，并插入注释：
# 在Rstudio中创建一个新的函数fit_temp，并将其保存为fit_temp.R，保存到R包Thermal.thresholds目录下“R”子文件夹
use_r("fit_temp")
# 在Rproj中打开此函数将光标定位在函数内部，然后依次点击菜单栏的 Code、Insert Roxygen Skeleton，就插入了注释框架。
# 对注释框架进行修改，并保存。
# 加载并测试函数功能
load_all()
ws<-read.csv("C:/shark/whaleshark/whaleshark_tmp2.csv",header = T,sep=",");str(ws)
fit_temp(ws,30)
# 讲数据写入包内
usethis::use_data(ws)
# 加入一个对数据ws的描述性数据data.R
use_r("data")
# 函数的注释写完了之后，它并不能被R自动识别，需要将其转义才可以。
devtools::document()
# 添加说明书Vignette，在myutils.Rmd中生成说明书的网页
usethis::use_vignette("myutils")
# 保存说明书
devtools::build_vignettes()
# 修改Description, 手动设置title, version, author, description,并加入 Imports,Suggests（需要或建议使用的其他R包），修改完后直接在Rproj中打开修改保存即可。
check()
# 保存建包的源代码
use_data_raw()
# 封装包
devtools::build()
# 上传至github
git config
use_git()
