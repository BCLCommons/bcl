from pymol.cgo import *
from math import *
from pymol import cmd
from pymol.vfont import plain
cmd.bg_color("white")
dendrogram=[]
labels_all = []

# label color
Magenta = [1,0,1]

cluster=[ CYLINDER, 45,2221.88,0,45,2221.88,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 15,1388.05,0,15,2221.88,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0000_final.pdb", "cluster_0000_final_8")
cmd.translate( [-9.76571, 1785.2, 24.7657], "cluster_0000_final_8")
cmd.show("cartoon", "cluster_0000_final_8" )
cmd.hide("lines", "cluster_0000_final_8" )
util.chainbow("cluster_0000_final_8")
size_girth_current_label=[]
axes=[[  0.0,0.0,41.6915],[0.0,41.6915,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[15,1596.51,16.1372],'2',5.21144,color=Magenta,axes=axes)
axes=[[  0.0,0.0,41.6915],[0.0,41.6915,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[15,1544.39,16.1372],'6.94025',5.21144,color=Magenta,axes=axes)
axes=[[  0.0,0.0,41.6915],[0.0,41.6915,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[15,1648.62,16.1372],'8',5.21144,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0000_final=[]
axes=[[  0.0,0.0,41.6915],[0.0,41.6915,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0000_final,plain,[15,2013.42,16.1372],'cluster_0000_final.pdb',5.21144,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0000_final)

cluster=[ CYLINDER, 65,1646.01,0,65,2221.88,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0003_final.pdb", "cluster_0003_final_9")
cmd.translate( [40.1574, 1914.1, 24.8426], "cluster_0003_final_9")
cmd.show("cartoon", "cluster_0003_final_9" )
cmd.hide("lines", "cluster_0003_final_9" )
util.chainbow("cluster_0003_final_9")
size_girth_current_label=[]
axes=[[  0.0,0.0,28.7935],[0.0,28.7935,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1789.98,12.9127],'3',3.59919,color=Magenta,axes=axes)
axes=[[  0.0,0.0,28.7935],[0.0,28.7935,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1753.99,12.9127],'8.23005',3.59919,color=Magenta,axes=axes)
axes=[[  0.0,0.0,28.7935],[0.0,28.7935,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1825.97,12.9127],'9',3.59919,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0003_final=[]
axes=[[  0.0,0.0,28.7935],[0.0,28.7935,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0003_final,plain,[65,2077.91,12.9127],'cluster_0003_final.pdb',3.59919,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0003_final)

cluster=[ CYLINDER, 15,2221.88,0,65,2221.88,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 15,2221.88,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 65,2221.88,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 5,1098.86,0,5,1388.05,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0000_final.pdb", "cluster_0000_final_1")
cmd.translate( [-19.7657, 1223.69, 24.7657], "cluster_0000_final_1")
cmd.show("cartoon", "cluster_0000_final_1" )
cmd.hide("lines", "cluster_0000_final_1" )
util.chainbow("cluster_0000_final_1")
size_girth_current_label=[]
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[5,1171.15,9.32922],'1',1.80747,color=Magenta,axes=axes)
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[5,1153.08,9.32922],'undef',1.80747,color=Magenta,axes=axes)
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[5,1189.23,9.32922],'1',1.80747,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0000_final=[]
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0000_final,plain,[5,1315.75,9.32922],'cluster_0000_final.pdb',1.80747,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0000_final)

cluster=[ CYLINDER, 25,1098.86,0,25,1388.05,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0004_final.pdb", "cluster_0004_final_5")
cmd.translate( [0.316583, 1223.77, 24.6834], "cluster_0004_final_5")
cmd.show("cartoon", "cluster_0004_final_5" )
cmd.hide("lines", "cluster_0004_final_5" )
util.chainbow("cluster_0004_final_5")
size_girth_current_label=[]
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[25,1171.15,9.32922],'1',1.80747,color=Magenta,axes=axes)
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[25,1153.08,9.32922],'undef',1.80747,color=Magenta,axes=axes)
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[25,1189.23,9.32922],'5',1.80747,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0004_final=[]
axes=[[  0.0,0.0,14.4597],[0.0,14.4597,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0004_final,plain,[25,1315.75,9.32922],'cluster_0004_final.pdb',1.80747,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0004_final)

cluster=[ CYLINDER, 5,1388.05,0,25,1388.05,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 5,1388.05,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 25,1388.05,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 55,1126.25,0,55,1646.01,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0001_final.pdb", "cluster_0001_final_7")
cmd.translate( [30.2732, 1366.4, 24.7268], "cluster_0001_final_7")
cmd.show("cartoon", "cluster_0001_final_7" )
cmd.hide("lines", "cluster_0001_final_7" )
util.chainbow("cluster_0001_final_7")
size_girth_current_label=[]
axes=[[  0.0,0.0,25.9882],[0.0,25.9882,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[55,1256.19,12.2113],'2',3.24852,color=Magenta,axes=axes)
axes=[[  0.0,0.0,25.9882],[0.0,25.9882,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[55,1223.7,12.2113],'5.63123',3.24852,color=Magenta,axes=axes)
axes=[[  0.0,0.0,25.9882],[0.0,25.9882,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[55,1288.67,12.2113],'7',3.24852,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0001_final=[]
axes=[[  0.0,0.0,25.9882],[0.0,25.9882,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0001_final,plain,[55,1516.07,12.2113],'cluster_0001_final.pdb',3.24852,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0001_final)

cluster=[ CYLINDER, 85,1098.86,0,85,1646.01,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0002_final.pdb", "cluster_0002_final_3")
cmd.translate( [60.4413, 1352.87, 24.5587], "cluster_0002_final_3")
cmd.show("cartoon", "cluster_0002_final_3" )
cmd.hide("lines", "cluster_0002_final_3" )
util.chainbow("cluster_0002_final_3")
size_girth_current_label=[]
axes=[[  0.0,0.0,27.3577],[0.0,27.3577,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[85,1235.64,12.5537],'1',3.41972,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.3577],[0.0,27.3577,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[85,1201.45,12.5537],'undef',3.41972,color=Magenta,axes=axes)
axes=[[  0.0,0.0,27.3577],[0.0,27.3577,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[85,1269.84,12.5537],'3',3.41972,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0002_final=[]
axes=[[  0.0,0.0,27.3577],[0.0,27.3577,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0002_final,plain,[85,1509.22,12.5537],'cluster_0002_final.pdb',3.41972,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0002_final)

cluster=[ CYLINDER, 55,1646.01,0,85,1646.01,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 55,1646.01,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 85,1646.01,0,5]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 45,1098.86,0,45,1126.25,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0001_final.pdb", "cluster_0001_final_2")
cmd.translate( [20.2732, 1092.82, 24.7268], "cluster_0001_final_2")
cmd.show("cartoon", "cluster_0001_final_2" )
cmd.hide("lines", "cluster_0001_final_2" )
util.chainbow("cluster_0001_final_2")
size_girth_current_label=[]
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[45,1105.7,6.05667],'1',0.171193,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[45,1103.99,6.05667],'undef',0.171193,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[45,1107.41,6.05667],'2',0.171193,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0001_final=[]
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0001_final,plain,[45,1119.4,6.05667],'cluster_0001_final.pdb',0.171193,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0001_final)

cluster=[ CYLINDER, 65,1098.86,0,65,1126.25,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cmd.load( "example/example_files/input/cluster//cluster_0003_final.pdb", "cluster_0003_final_4")
cmd.translate( [40.1574, 1092.71, 24.8426], "cluster_0003_final_4")
cmd.show("cartoon", "cluster_0003_final_4" )
cmd.hide("lines", "cluster_0003_final_4" )
util.chainbow("cluster_0003_final_4")
size_girth_current_label=[]
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1105.7,6.05667],'1',0.171193,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1103.99,6.05667],'undef',0.171193,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_current_label,plain,[65,1107.41,6.05667],'4',0.171193,color=Magenta,axes=axes)
labels_all.extend( size_girth_current_label)
center_cluster_0003_final=[]
axes=[[  0.0,0.0,1.36954],[0.0,1.36954,0.0],[0.0,0.0,0.0]]
cyl_text(center_cluster_0003_final,plain,[65,1119.4,6.05667],'cluster_0003_final.pdb',0.171193,color=Magenta,axes=axes)
labels_all.extend( center_cluster_0003_final)

cluster=[ CYLINDER, 45,1126.25,0,65,1126.25,0,5,0,0,0,0,0,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 45,1126.25,0,5]
dendrogram.extend( cluster)
cluster=[ COLOR, 0,0,0, SPHERE, 65,1126.25,0,5]
dendrogram.extend( cluster)
cmd.load_cgo( dendrogram, 'dendrogram')
cmd.load_cgo( labels_all, 'labels_all')

axes=[[ 43.8254 ,0.0,0.0],[0.0,21.9127,0.0],[0.0,0.0,0.0]]
labels=[]
# girth label color
Black = [0,0,0]

cyl_text(labels,plain,[125.478,2210.92,0],'11.109',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,2101.36,0],'10.562',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1991.8,0],'10.014',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1882.23,0],'9.466',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1772.67,0],'8.918',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1663.11,0],'8.370',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1553.54,0],'7.822',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1443.98,0],'7.275',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1334.42,0],'6.727',2.73909,color=Black,axes=axes)
cyl_text(labels,plain,[125.478,1224.85,0],'6.179',2.73909,color=Black,axes=axes)
cmd.load_cgo(labels,'labels')

