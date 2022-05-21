from pymol.cgo import *
from math import *
from pymol import cmd
from pymol.vfont import plain
cmd.bg_color("white")
dendrogram=[]
labels_all = []

# label color
Magenta = [1,0,1]

cluster=[ CYLINDER, 350,1400,0,350,2000,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 50,-50,0,50,1400,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_50_0=[]
axes=[[  0.0,0.0,72.5],[0.0,72.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,312.5,75.2679],'1',9.0625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,72.5],[0.0,72.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,221.875,75.2679],'undef',9.0625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,72.5],[0.0,72.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,403.125,75.2679],'4',9.0625,color=Magenta,axes=axes)
labels_all.extend( size_girth_50_0)
center__50_0=[]
axes=[[  0.0,0.0,72.5],[0.0,72.5,0.0],[0.0,0.0,0.0]]
cyl_text(center__50_0,plain,[50,1037.5,75.2679],'d',9.0625,color=Magenta,axes=axes)
labels_all.extend( center__50_0)

cluster=[ CYLINDER, 450,950,0,450,1400,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_450_9_5=[]
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_9_5,plain,[450,1062.5,62.7679],'3',2.8125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_9_5,plain,[450,1034.38,62.7679],'9.5',2.8125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_9_5,plain,[450,1090.62,62.7679],'7',2.8125,color=Magenta,axes=axes)
labels_all.extend( size_girth_450_9_5)
center__450_9_5=[]
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(center__450_9_5,plain,[450,1287.5,62.7679],'a',2.8125,color=Magenta,axes=axes)
labels_all.extend( center__450_9_5)

cluster=[ CYLINDER, 50,1400,0,450,1400,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 50,1400,0,50]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 450,1400,0,50]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 350,500,0,350,950,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_350_5=[]
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_350_5,plain,[350,612.5,62.7679],'2',2.8125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_350_5,plain,[350,584.375,62.7679],'5',2.8125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_350_5,plain,[350,640.625,62.7679],'6',2.8125,color=Magenta,axes=axes)
labels_all.extend( size_girth_350_5)
center__350_5=[]
axes=[[  0.0,0.0,22.5],[0.0,22.5,0.0],[0.0,0.0,0.0]]
cyl_text(center__350_5,plain,[350,837.5,62.7679],'a',2.8125,color=Magenta,axes=axes)
labels_all.extend( center__350_5)

cluster=[ CYLINDER, 650,-50,0,650,950,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_650_0=[]
axes=[[  0.0,0.0,50],[0.0,50,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_650_0,plain,[650,200,69.6429],'1',6.25,color=Magenta,axes=axes)
axes=[[  0.0,0.0,50],[0.0,50,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_650_0,plain,[650,137.5,69.6429],'undef',6.25,color=Magenta,axes=axes)
axes=[[  0.0,0.0,50],[0.0,50,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_650_0,plain,[650,262.5,69.6429],'3',6.25,color=Magenta,axes=axes)
labels_all.extend( size_girth_650_0)
center__650_0=[]
axes=[[  0.0,0.0,50],[0.0,50,0.0],[0.0,0.0,0.0]]
cyl_text(center__650_0,plain,[650,700,69.6429],'c',6.25,color=Magenta,axes=axes)
labels_all.extend( center__650_0)

cluster=[ CYLINDER, 350,950,0,650,950,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 350,950,0,50]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 650,950,0,50]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 250,-50,0,250,500,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_250_0=[]
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_250_0,plain,[250,87.5,64.0179],'1',3.4375,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_250_0,plain,[250,53.125,64.0179],'undef',3.4375,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_250_0,plain,[250,121.875,64.0179],'1',3.4375,color=Magenta,axes=axes)
labels_all.extend( size_girth_250_0)
center__250_0=[]
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(center__250_0,plain,[250,362.5,64.0179],'a',3.4375,color=Magenta,axes=axes)
labels_all.extend( center__250_0)

cluster=[ CYLINDER, 450,-50,0,450,500,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
size_girth_450_0=[]
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_0,plain,[450,87.5,64.0179],'1',3.4375,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_0,plain,[450,53.125,64.0179],'undef',3.4375,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_450_0,plain,[450,121.875,64.0179],'2',3.4375,color=Magenta,axes=axes)
labels_all.extend( size_girth_450_0)
center__450_0=[]
axes=[[  0.0,0.0,27.5],[0.0,27.5,0.0],[0.0,0.0,0.0]]
cyl_text(center__450_0,plain,[450,362.5,64.0179],'b',3.4375,color=Magenta,axes=axes)
labels_all.extend( center__450_0)

cluster=[ CYLINDER, 250,500,0,450,500,0,50,1,0.5,0,1,0.5,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 250,500,0,50]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.5,0, SPHERE, 450,500,0,50]
dendrogram.extend( cluster)
cmd.load_cgo( dendrogram, 'dendrogram')
cmd.load_cgo( labels_all, 'labels_all')

axes=[[ 80 ,0.0,0.0],[0.0,40,0.0],[0.0,0.0,0.0]]
labels=[]
# girth label color
Black = [0,0,0]

cyl_text(labels,plain,[1010,1980,0],'20.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,1780,0],'18.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,1580,0],'16.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,1380,0],'14.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,1180,0],'12.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,980,0],'10.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,780,0],'8.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,580,0],'6.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,380,0],'4.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,180,0],'2.000',5,color=Black,axes=axes)
cyl_text(labels,plain,[1010,-20,0],'0.000',5,color=Black,axes=axes)
cmd.load_cgo(labels,'labels')

