##後�?�クラスタリング等�?�ためにpackeage用�?
install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library("cluster")
library("factoextra")
library("gplots")

#�?散の高いm個�?�遺伝子を抽出し、m×mの行�?�作製
data<-read.table("coadExp.tumor.tab", row.names=1, header=T, sep="\t", quote="")　#col.name=1ではheader�?定できな�?(多�??空�?�?から)
data2<-data　#コピ�?�して�?�?ータを温�?
data2$"variance"<-apply(data2, 1, var) #dataの�?行にvarを演算し、varianceリストを追�?
head(data2) #結果を一部表示
data2_order<-data2[order(data2$variance, decreasing=T), ] #varianceを降�?で並び替�?
head(data2_order$variance) #varianceと値と降�??であることを確�?
data3<-head(data2_order, 1000) #variance上�?1000遺伝子を抽出

#K-meansでクラスタリング
data3_gene_km<-kmeans(data3, center=4, iter.max=100, nstart=5) #4群�?け、最大100回収束させる過程�?5�?
data3_case_km<-kmeans(t(data3), center=4, iter.max=100, nstart=5)

fviz_cluster(data3_gene_km, data = data3,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "t", 　　　　　　　# 楕�??を追�? 
             star.plot = TRUE, 　　　　　　　　# クラスターのcentroidから�?�?子に線を引く
             repel = TRUE, 　　　　　　　　　　# ラベルが重ならな�?ようにする
             ggtheme = theme_minimal()
)

# kmeansクラスタリングを�?�ロ�?トする（症例�?
fviz_cluster(data3_case_km, data =t(data3),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "t", 　　　　　　　# 楕�??を追�?
             star.plot = TRUE, 　　　　　　　　# クラスターのcentroidから�?�?子に線を引く 
             repel = TRUE, 　　　　　　　　　　# ラベルが重ならな�?ようにする
             ggtheme = theme_minimal()
)

#階層�?クラスタリングを適用、ヒート�?�ップで可視化
dist.gene<-dist(data3, method="euclidean") #ユークリ�?ド距離の算�?�
dist.case<-dist(t(data3), method="euclidean") #ユークリ�?ド距離の算�?�

fviz_dist(dist.gene) # 遺伝子間距離を視覚化
fviz_dist(dist.case, lab_size = 3) # �?例間距離を視覚化

c1<-hclust(d=dist.gene, method="ward.D2")
c2<-hclust(d=dist.case, method="ward.D2")

rowcolor <- c("1" = "firebrick2", "2" = "dodgerblue2", "3" = "palegreen3", "4" = "gold") # look up tableで�?clusterに配色
colcolor <- c("1" = "firebrick2", "2" = "dodgerblue2", "3" = "palegreen3", "4" = "gold")

heatmap.2(as.matrix(data3), scale = "none", 
          Rowv = as.dendrogram(c1),  　　　　
          Colv = as.dendrogram(c2),　　　　　
          col = bluered(100),　     # heatmap本体�?�色
          trace = "none", 　　　　　# �?cell�?の�?(中央からのズレで発現レベル表現)を消去
          density.info = "none",　　
          cexRow = 0.6, cexCol = 0.6, # �?字�?�大きさ調整
          RowSideColors = rowcolor[data3_gene_km$cluster], # 遺伝子�?�Kmeans clusteringの配色look up
          ColSideColors = colcolor[data3_case_km$cluster]) # �?例�?�kmeans clusteringの配色look up