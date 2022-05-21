from pymol.cgo import *
from math import *
from pymol import cmd
from pymol.vfont import plain
cmd.bg_color("white")
dendrogram=[]
labels_all = []

# label color
Magenta = [1,0,1]

cluster=[ CYLINDER, 45,163.118,0,45,163.118,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 5,33.1947,0,5,163.118,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/2.mol", "2_3")
cmd.translate( [5, 98.1564, 17.407], "2_3")
size_girth_2_3=[]
axes=[[  0.0,0.0,6.49616],[0.0,6.49616,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_3,plain,[5,65.6756,7.33833],'1',0.81202,color=Magenta,axes=axes)
axes=[[  0.0,0.0,6.49616],[0.0,6.49616,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_3,plain,[5,57.5554,7.33833],'undef',0.81202,color=Magenta,axes=axes)
axes=[[  0.0,0.0,6.49616],[0.0,6.49616,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_3,plain,[5,73.7958,7.33833],'3',0.81202,color=Magenta,axes=axes)
labels_all.extend( size_girth_2_3)
center_2=[]
axes=[[  0.0,0.0,6.49616],[0.0,6.49616,0.0],[0.0,0.0,0.0]]
cyl_text(center_2,plain,[5,130.637,7.33833],'2',0.81202,color=Magenta,axes=axes)
labels_all.extend( center_2)

cluster=[ CYLINDER, 55,148.902,0,55,163.118,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/1.mol", "1_9")
cmd.translate( [55, 156.01, 17.407], "1_9")
size_girth_1_9=[]
axes=[[  0.0,0.0,0.71081],[0.0,0.71081,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_9,plain,[55,152.456,5.89199],'4',0.0888512,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.71081],[0.0,0.71081,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_9,plain,[55,151.567,5.89199],'0.744509',0.0888512,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.71081],[0.0,0.71081,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_9,plain,[55,153.344,5.89199],'9',0.0888512,color=Magenta,axes=axes)
labels_all.extend( size_girth_1_9)
center_1=[]
axes=[[  0.0,0.0,0.71081],[0.0,0.71081,0.0],[0.0,0.0,0.0]]
cyl_text(center_1,plain,[55,159.564,5.89199],'1',0.0888512,color=Magenta,axes=axes)
labels_all.extend( center_1)

cluster=[ CYLINDER, 5,163.118,0,55,163.118,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 5,163.118,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 55,163.118,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 35,48.6486,0,35,148.902,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/0.mol", "0_8")
cmd.translate( [35, 98.7752, 17.407], "0_8")
size_girth_0_8=[]
axes=[[  0.0,0.0,5.01266],[0.0,5.01266,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_8,plain,[35,73.7119,6.96745],'2',0.626583,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.01266],[0.0,5.01266,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_8,plain,[35,67.4461,6.96745],'0.243243',0.626583,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.01266],[0.0,5.01266,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_8,plain,[35,79.9777,6.96745],'8',0.626583,color=Magenta,axes=axes)
labels_all.extend( size_girth_0_8)
center_0=[]
axes=[[  0.0,0.0,5.01266],[0.0,5.01266,0.0],[0.0,0.0,0.0]]
cyl_text(center_0,plain,[35,123.839,6.96745],'0',0.626583,color=Magenta,axes=axes)
labels_all.extend( center_0)

cluster=[ CYLINDER, 75,36.3636,0,75,148.902,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/1.mol", "1_7")
cmd.translate( [75, 92.6327, 17.407], "1_7")
size_girth_1_7=[]
axes=[[  0.0,0.0,5.62691],[0.0,5.62691,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_7,plain,[75,64.4981,7.12101],'2',0.703364,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.62691],[0.0,5.62691,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_7,plain,[75,57.4645,7.12101],'0.181818',0.703364,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.62691],[0.0,5.62691,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_7,plain,[75,71.5318,7.12101],'7',0.703364,color=Magenta,axes=axes)
labels_all.extend( size_girth_1_7)
center_1=[]
axes=[[  0.0,0.0,5.62691],[0.0,5.62691,0.0],[0.0,0.0,0.0]]
cyl_text(center_1,plain,[75,120.767,7.12101],'1',0.703364,color=Magenta,axes=axes)
labels_all.extend( center_1)

cluster=[ CYLINDER, 35,148.902,0,75,148.902,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 35,148.902,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 75,148.902,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 25,33.1947,0,25,48.6486,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/0.mol", "0_1")
cmd.translate( [25, 40.9217, 17.407], "0_1")
size_girth_0_1=[]
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_1,plain,[25,37.0582,5.90746],'1',0.0965866,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_1,plain,[25,36.0923,5.90746],'undef',0.0965866,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_0_1,plain,[25,38.0241,5.90746],'1',0.0965866,color=Magenta,axes=axes)
labels_all.extend( size_girth_0_1)
center_0=[]
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(center_0,plain,[25,44.7851,5.90746],'0',0.0965866,color=Magenta,axes=axes)
labels_all.extend( center_0)

cluster=[ CYLINDER, 45,33.1947,0,45,48.6486,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/3.mol", "3_4")
cmd.translate( [45, 40.9217, 17.407], "3_4")
size_girth_3_4=[]
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_3_4,plain,[45,37.0582,5.90746],'1',0.0965866,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_3_4,plain,[45,36.0923,5.90746],'undef',0.0965866,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_3_4,plain,[45,38.0241,5.90746],'4',0.0965866,color=Magenta,axes=axes)
labels_all.extend( size_girth_3_4)
center_3=[]
axes=[[  0.0,0.0,0.772693],[0.0,0.772693,0.0],[0.0,0.0,0.0]]
cyl_text(center_3,plain,[45,44.7851,5.90746],'3',0.0965866,color=Magenta,axes=axes)
labels_all.extend( center_3)

cluster=[ CYLINDER, 25,48.6486,0,45,48.6486,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 25,48.6486,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 45,48.6486,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 65,33.1947,0,65,36.3636,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/1.mol", "1_2")
cmd.translate( [65, 34.7792, 17.407], "1_2")
size_girth_1_2=[]
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_2,plain,[65,33.987,5.7539],'1',0.0198054,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_2,plain,[65,33.7889,5.7539],'undef',0.0198054,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_1_2,plain,[65,34.185,5.7539],'2',0.0198054,color=Magenta,axes=axes)
labels_all.extend( size_girth_1_2)
center_1=[]
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(center_1,plain,[65,35.5714,5.7539],'1',0.0198054,color=Magenta,axes=axes)
labels_all.extend( center_1)

cluster=[ CYLINDER, 85,33.1947,0,85,36.3636,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/output/cluster/4.mol", "4_5")
cmd.translate( [85, 34.7792, 17.407], "4_5")
size_girth_4_5=[]
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_4_5,plain,[85,33.987,5.7539],'1',0.0198054,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_4_5,plain,[85,33.7889,5.7539],'undef',0.0198054,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_4_5,plain,[85,34.185,5.7539],'5',0.0198054,color=Magenta,axes=axes)
labels_all.extend( size_girth_4_5)
center_4=[]
axes=[[  0.0,0.0,0.158443],[0.0,0.158443,0.0],[0.0,0.0,0.0]]
cyl_text(center_4,plain,[85,35.5714,5.7539],'4',0.0198054,color=Magenta,axes=axes)
labels_all.extend( center_4)

cluster=[ CYLINDER, 65,36.3636,0,85,36.3636,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 65,36.3636,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 85,36.3636,0,5]
dendrogram.extend( cluster)
cmd.load_cgo( dendrogram, 'dendrogram')
cmd.load_cgo( labels_all, 'labels_all')

axes=[[ 5.07018 ,0.0,0.0],[0.0,2.53509,0.0],[0.0,0.0,0.0]]
labels=[]
# girth label color
Black = [0,0,0]

cyl_text(labels,plain,[120.634,161.85,0],'0.816',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,149.175,0],'0.752',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,136.5,0],'0.689',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,123.824,0],'0.625',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,111.149,0],'0.562',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,98.4733,0],'0.499',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,85.7978,0],'0.435',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,73.1224,0],'0.372',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,60.4469,0],'0.309',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,47.7715,0],'0.245',0.316886,color=Black,axes=axes)
cyl_text(labels,plain,[120.634,35.0961,0],'0.182',0.316886,color=Black,axes=axes)
cmd.load_cgo(labels,'labels')

