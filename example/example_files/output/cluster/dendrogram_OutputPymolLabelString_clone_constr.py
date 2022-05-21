from pymol.cgo import *
from math import *
from pymol import cmd
from pymol.vfont import plain
cmd.bg_color("white")
dendrogram=[]
labels_all = []

# label color
Magenta = [1,0,1]

cluster=[ CYLINDER, 35,2800,0,35,2800,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 5,955,0,5,2800,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_5_0=[]
axes=[[  0.0,0.0,92.25],[0.0,92.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_5_0,plain,[5,1416.25,28.7768],'1',11.5312,color=Magenta,axes=axes)
axes=[[  0.0,0.0,92.25],[0.0,92.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_5_0,plain,[5,1300.94,28.7768],'undef',11.5312,color=Magenta,axes=axes)
axes=[[  0.0,0.0,92.25],[0.0,92.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_5_0,plain,[5,1531.56,28.7768],'4',11.5312,color=Magenta,axes=axes)
labels_all.extend( size_girth_5_0)
center__5_0=[]
axes=[[  0.0,0.0,92.25],[0.0,92.25,0.0],[0.0,0.0,0.0]]
cyl_text(center__5_0,plain,[5,2338.75,28.7768],'d',11.5312,color=Magenta,axes=axes)
labels_all.extend( center__5_0)

cluster=[ CYLINDER, 45,1900,0,45,2800,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_45_9_5=[]
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_9_5,plain,[45,2125,16.9643],'3',5.625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_9_5,plain,[45,2068.75,16.9643],'9.5',5.625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_9_5,plain,[45,2181.25,16.9643],'7',5.625,color=Magenta,axes=axes)
labels_all.extend( size_girth_45_9_5)
center__45_9_5=[]
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(center__45_9_5,plain,[45,2575,16.9643],'a',5.625,color=Magenta,axes=axes)
labels_all.extend( center__45_9_5)

cluster=[ CYLINDER, 5,2800,0,45,2800,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 5,2800,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 45,2800,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 35,1000,0,35,1900,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_35_5=[]
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_35_5,plain,[35,1225,16.9643],'2',5.625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_35_5,plain,[35,1168.75,16.9643],'5',5.625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_35_5,plain,[35,1281.25,16.9643],'6',5.625,color=Magenta,axes=axes)
labels_all.extend( size_girth_35_5)
center__35_5=[]
axes=[[  0.0,0.0,45],[0.0,45,0.0],[0.0,0.0,0.0]]
cyl_text(center__35_5,plain,[35,1675,16.9643],'a',5.625,color=Magenta,axes=axes)
labels_all.extend( center__35_5)

cluster=[ CYLINDER, 65,955,0,65,1900,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_65_0=[]
axes=[[  0.0,0.0,47.25],[0.0,47.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_65_0,plain,[65,1191.25,17.5268],'1',5.90625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,47.25],[0.0,47.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_65_0,plain,[65,1132.19,17.5268],'undef',5.90625,color=Magenta,axes=axes)
axes=[[  0.0,0.0,47.25],[0.0,47.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_65_0,plain,[65,1250.31,17.5268],'3',5.90625,color=Magenta,axes=axes)
labels_all.extend( size_girth_65_0)
center__65_0=[]
axes=[[  0.0,0.0,47.25],[0.0,47.25,0.0],[0.0,0.0,0.0]]
cyl_text(center__65_0,plain,[65,1663.75,17.5268],'c',5.90625,color=Magenta,axes=axes)
labels_all.extend( center__65_0)

cluster=[ CYLINDER, 35,1900,0,65,1900,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 35,1900,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 65,1900,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 25,955,0,25,1000,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_25_0=[]
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_25_0,plain,[25,966.25,6.27679],'1',0.28125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_25_0,plain,[25,963.438,6.27679],'undef',0.28125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_25_0,plain,[25,969.062,6.27679],'1',0.28125,color=Magenta,axes=axes)
labels_all.extend( size_girth_25_0)
center__25_0=[]
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(center__25_0,plain,[25,988.75,6.27679],'a',0.28125,color=Magenta,axes=axes)
labels_all.extend( center__25_0)

cluster=[ CYLINDER, 45,955,0,45,1000,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
size_girth_45_0=[]
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_0,plain,[45,966.25,6.27679],'1',0.28125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_0,plain,[45,963.438,6.27679],'undef',0.28125,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_45_0,plain,[45,969.062,6.27679],'2',0.28125,color=Magenta,axes=axes)
labels_all.extend( size_girth_45_0)
center__45_0=[]
axes=[[  0.0,0.0,2.25],[0.0,2.25,0.0],[0.0,0.0,0.0]]
cyl_text(center__45_0,plain,[45,988.75,6.27679],'b',0.28125,color=Magenta,axes=axes)
labels_all.extend( center__45_0)

cluster=[ CYLINDER, 25,1000,0,45,1000,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 25,1000,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 45,1000,0,5]
dendrogram.extend( cluster)
cmd.load_cgo( dendrogram, 'dendrogram')
cmd.load_cgo( labels_all, 'labels_all')

axes=[[ 72 ,0.0,0.0],[0.0,36,0.0],[0.0,0.0,0.0]]
labels=[]
# girth label color
Black = [0,0,0]

cyl_text(labels,plain,[109,2782,0],'14.000',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,2602,0],'13.100',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,2422,0],'12.200',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,2242,0],'11.300',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,2062,0],'10.400',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,1882,0],'9.500',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,1702,0],'8.600',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,1522,0],'7.700',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,1342,0],'6.800',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,1162,0],'5.900',4.5,color=Black,axes=axes)
cyl_text(labels,plain,[109,982,0],'5.000',4.5,color=Black,axes=axes)
cmd.load_cgo(labels,'labels')

