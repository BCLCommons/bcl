from pymol.cgo import *
from math import *
from pymol import cmd
from pymol.vfont import plain
cmd.bg_color("white")
dendrogram=[]
labels_all = []

# label color
Magenta = [1,0,1]

cluster=[ CYLINDER, 110,654.591,0,110,654.591,0,2,1,0.658154,0,1,0.658154,0]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 2,144.153,0,2,654.591,0,2,1,0.580617,0,1,0.580617,0]
dendrogram.extend( cluster)
size_girth_2_0=[]
axes=[[  0.0,0.0,25.5219],[0.0,25.5219,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_0,plain,[2,271.763,8.66619],'1',3.19024,color=Magenta,axes=axes)
axes=[[  0.0,0.0,25.5219],[0.0,25.5219,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_0,plain,[2,239.86,8.66619],'undef',3.19024,color=Magenta,axes=axes)
axes=[[  0.0,0.0,25.5219],[0.0,25.5219,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_2_0,plain,[2,303.665,8.66619],'7',3.19024,color=Magenta,axes=axes)
labels_all.extend( size_girth_2_0)
center__2_0=[]
axes=[[  0.0,0.0,25.5219],[0.0,25.5219,0.0],[0.0,0.0,0.0]]
cyl_text(center__2_0,plain,[2,526.982,8.66619],'1000_0006',3.19024,color=Magenta,axes=axes)
labels_all.extend( center__2_0)

cluster=[ CYLINDER, 122,631,0,122,654.591,0,2,1,0.666769,0,1,0.666769,0]
dendrogram.extend( cluster)
size_girth_122_6_31=[]
axes=[[  0.0,0.0,1.17955],[0.0,1.17955,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_6_31,plain,[122,636.898,2.5806],'9',0.147444,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.17955],[0.0,1.17955,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_6_31,plain,[122,635.423,2.5806],'6.31',0.147444,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.17955],[0.0,1.17955,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_6_31,plain,[122,638.372,2.5806],'19',0.147444,color=Magenta,axes=axes)
labels_all.extend( size_girth_122_6_31)
center__122_6_31=[]
axes=[[  0.0,0.0,1.17955],[0.0,1.17955,0.0],[0.0,0.0,0.0]]
cyl_text(center__122_6_31,plain,[122,648.693,2.5806],'1000_0004',0.147444,color=Magenta,axes=axes)
labels_all.extend( center__122_6_31)

cluster=[ CYLINDER, 2,654.591,0,122,654.591,0,2,1,0.658154,0,1,0.658154,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.658154,0, SPHERE, 2,654.591,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.658154,0, SPHERE, 122,654.591,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 110,539.709,0,110,631,0,2,1,0.742391,0,1,0.742391,0]
dendrogram.extend( cluster)
size_girth_110_5_39709=[]
axes=[[  0.0,0.0,4.56455],[0.0,4.56455,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_110_5_39709,plain,[110,562.532,3.42685],'8',0.570569,color=Magenta,axes=axes)
axes=[[  0.0,0.0,4.56455],[0.0,4.56455,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_110_5_39709,plain,[110,556.826,3.42685],'5.39709',0.570569,color=Magenta,axes=axes)
axes=[[  0.0,0.0,4.56455],[0.0,4.56455,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_110_5_39709,plain,[110,568.237,3.42685],'18',0.570569,color=Magenta,axes=axes)
labels_all.extend( size_girth_110_5_39709)
center__110_5_39709=[]
axes=[[  0.0,0.0,4.56455],[0.0,4.56455,0.0],[0.0,0.0,0.0]]
cyl_text(center__110_5_39709,plain,[110,608.177,3.42685],'1000_0004',0.570569,color=Magenta,axes=axes)
labels_all.extend( center__110_5_39709)

cluster=[ CYLINDER, 218,144.153,0,218,631,0,2,1,0.0617913,0,1,0.0617913,0]
dendrogram.extend( cluster)
size_girth_218_0=[]
axes=[[  0.0,0.0,24.3423],[0.0,24.3423,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_218_0,plain,[218,265.865,8.3713],'1',3.04279,color=Magenta,axes=axes)
axes=[[  0.0,0.0,24.3423],[0.0,24.3423,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_218_0,plain,[218,235.437,8.3713],'undef',3.04279,color=Magenta,axes=axes)
axes=[[  0.0,0.0,24.3423],[0.0,24.3423,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_218_0,plain,[218,296.293,8.3713],'4',3.04279,color=Magenta,axes=axes)
labels_all.extend( size_girth_218_0)
center__218_0=[]
axes=[[  0.0,0.0,24.3423],[0.0,24.3423,0.0],[0.0,0.0,0.0]]
cyl_text(center__218_0,plain,[218,509.288,8.3713],'1000_0003',3.04279,color=Magenta,axes=axes)
labels_all.extend( center__218_0)

cluster=[ CYLINDER, 110,631,0,218,631,0,2,1,0.666769,0,1,0.666769,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.666769,0, SPHERE, 110,631,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.666769,0, SPHERE, 218,631,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 98,482.608,0,98,539.709,0,2,1,0.75598,0,1,0.75598,0]
dendrogram.extend( cluster)
size_girth_98_4_82608=[]
axes=[[  0.0,0.0,2.85505],[0.0,2.85505,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_4_82608,plain,[98,496.883,2.99948],'7',0.356881,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.85505],[0.0,2.85505,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_4_82608,plain,[98,493.314,2.99948],'4.82608',0.356881,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.85505],[0.0,2.85505,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_4_82608,plain,[98,500.452,2.99948],'17',0.356881,color=Magenta,axes=axes)
labels_all.extend( size_girth_98_4_82608)
center__98_4_82608=[]
axes=[[  0.0,0.0,2.85505],[0.0,2.85505,0.0],[0.0,0.0,0.0]]
cyl_text(center__98_4_82608,plain,[98,525.434,2.99948],'1000_0004',0.356881,color=Magenta,axes=axes)
labels_all.extend( center__98_4_82608)

cluster=[ CYLINDER, 194,144.153,0,194,539.709,0,2,1,0.64727,0,1,0.64727,0]
dendrogram.extend( cluster)
size_girth_194_0=[]
axes=[[  0.0,0.0,19.7778],[0.0,19.7778,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_194_0,plain,[194,243.042,7.23016],'1',2.47222,color=Magenta,axes=axes)
axes=[[  0.0,0.0,19.7778],[0.0,19.7778,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_194_0,plain,[194,218.32,7.23016],'undef',2.47222,color=Magenta,axes=axes)
axes=[[  0.0,0.0,19.7778],[0.0,19.7778,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_194_0,plain,[194,267.764,7.23016],'1',2.47222,color=Magenta,axes=axes)
labels_all.extend( size_girth_194_0)
center__194_0=[]
axes=[[  0.0,0.0,19.7778],[0.0,19.7778,0.0],[0.0,0.0,0.0]]
cyl_text(center__194_0,plain,[194,440.82,7.23016],'1000_0000',2.47222,color=Magenta,axes=axes)
labels_all.extend( center__194_0)

cluster=[ CYLINDER, 98,539.709,0,194,539.709,0,2,1,0.742391,0,1,0.742391,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.742391,0, SPHERE, 98,539.709,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.742391,0, SPHERE, 194,539.709,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 86,386.311,0,86,482.608,0,2,1,0.761129,0,1,0.761129,0]
dendrogram.extend( cluster)
size_girth_86_3_86311=[]
axes=[[  0.0,0.0,4.81485],[0.0,4.81485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_86_3_86311,plain,[86,410.385,3.48943],'6',0.601856,color=Magenta,axes=axes)
axes=[[  0.0,0.0,4.81485],[0.0,4.81485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_86_3_86311,plain,[86,404.367,3.48943],'3.86311',0.601856,color=Magenta,axes=axes)
axes=[[  0.0,0.0,4.81485],[0.0,4.81485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_86_3_86311,plain,[86,416.404,3.48943],'16',0.601856,color=Magenta,axes=axes)
labels_all.extend( size_girth_86_3_86311)
center__86_3_86311=[]
axes=[[  0.0,0.0,4.81485],[0.0,4.81485,0.0],[0.0,0.0,0.0]]
cyl_text(center__86_3_86311,plain,[86,458.534,3.48943],'1000_0004',0.601856,color=Magenta,axes=axes)
labels_all.extend( center__86_3_86311)

cluster=[ CYLINDER, 170,144.153,0,170,482.608,0,2,1,0.725087,0,1,0.725087,0]
dendrogram.extend( cluster)
size_girth_170_0=[]
axes=[[  0.0,0.0,16.9227],[0.0,16.9227,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_170_0,plain,[170,228.767,6.5164],'1',2.11534,color=Magenta,axes=axes)
axes=[[  0.0,0.0,16.9227],[0.0,16.9227,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_170_0,plain,[170,207.614,6.5164],'undef',2.11534,color=Magenta,axes=axes)
axes=[[  0.0,0.0,16.9227],[0.0,16.9227,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_170_0,plain,[170,249.92,6.5164],'3',2.11534,color=Magenta,axes=axes)
labels_all.extend( size_girth_170_0)
center__170_0=[]
axes=[[  0.0,0.0,16.9227],[0.0,16.9227,0.0],[0.0,0.0,0.0]]
cyl_text(center__170_0,plain,[170,397.994,6.5164],'1000_0002',2.11534,color=Magenta,axes=axes)
labels_all.extend( center__170_0)

cluster=[ CYLINDER, 86,482.608,0,170,482.608,0,2,1,0.75598,0,1,0.75598,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.75598,0, SPHERE, 86,482.608,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.75598,0, SPHERE, 170,482.608,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 74,359.632,0,74,386.311,0,2,1,0.787497,0,1,0.787497,0]
dendrogram.extend( cluster)
size_girth_74_3_59632=[]
axes=[[  0.0,0.0,1.33395],[0.0,1.33395,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_3_59632,plain,[74,366.302,2.6192],'5',0.166744,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.33395],[0.0,1.33395,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_3_59632,plain,[74,364.634,2.6192],'3.59632',0.166744,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.33395],[0.0,1.33395,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_3_59632,plain,[74,367.969,2.6192],'15',0.166744,color=Magenta,axes=axes)
labels_all.extend( size_girth_74_3_59632)
center__74_3_59632=[]
axes=[[  0.0,0.0,1.33395],[0.0,1.33395,0.0],[0.0,0.0,0.0]]
cyl_text(center__74_3_59632,plain,[74,379.641,2.6192],'1000_0004',0.166744,color=Magenta,axes=axes)
labels_all.extend( center__74_3_59632)

cluster=[ CYLINDER, 146,144.153,0,146,386.311,0,2,1,0.629287,0,1,0.629287,0]
dendrogram.extend( cluster)
size_girth_146_0=[]
axes=[[  0.0,0.0,12.1079],[0.0,12.1079,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_146_0,plain,[146,204.693,5.31269],'1',1.51349,color=Magenta,axes=axes)
axes=[[  0.0,0.0,12.1079],[0.0,12.1079,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_146_0,plain,[146,189.558,5.31269],'undef',1.51349,color=Magenta,axes=axes)
axes=[[  0.0,0.0,12.1079],[0.0,12.1079,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_146_0,plain,[146,219.828,5.31269],'10',1.51349,color=Magenta,axes=axes)
labels_all.extend( size_girth_146_0)
center__146_0=[]
axes=[[  0.0,0.0,12.1079],[0.0,12.1079,0.0],[0.0,0.0,0.0]]
cyl_text(center__146_0,plain,[146,325.772,5.31269],'1000_0009',1.51349,color=Magenta,axes=axes)
labels_all.extend( center__146_0)

cluster=[ CYLINDER, 74,386.311,0,146,386.311,0,2,1,0.761129,0,1,0.761129,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.761129,0, SPHERE, 74,386.311,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.761129,0, SPHERE, 146,386.311,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 62,324.795,0,62,359.632,0,2,1,0.825811,0,1,0.825811,0]
dendrogram.extend( cluster)
size_girth_62_3_24795=[]
axes=[[  0.0,0.0,1.74185],[0.0,1.74185,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_62_3_24795,plain,[62,333.504,2.72118],'4',0.217731,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.74185],[0.0,1.74185,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_62_3_24795,plain,[62,331.327,2.72118],'3.24795',0.217731,color=Magenta,axes=axes)
axes=[[  0.0,0.0,1.74185],[0.0,1.74185,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_62_3_24795,plain,[62,335.682,2.72118],'14',0.217731,color=Magenta,axes=axes)
labels_all.extend( size_girth_62_3_24795)
center__62_3_24795=[]
axes=[[  0.0,0.0,1.74185],[0.0,1.74185,0.0],[0.0,0.0,0.0]]
cyl_text(center__62_3_24795,plain,[62,350.923,2.72118],'1000_0004',0.217731,color=Magenta,axes=axes)
labels_all.extend( center__62_3_24795)

cluster=[ CYLINDER, 122,144.153,0,122,359.632,0,2,1,0.634243,0,1,0.634243,0]
dendrogram.extend( cluster)
size_girth_122_0=[]
axes=[[  0.0,0.0,10.7739],[0.0,10.7739,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_0,plain,[122,198.023,4.9792],'1',1.34674,color=Magenta,axes=axes)
axes=[[  0.0,0.0,10.7739],[0.0,10.7739,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_0,plain,[122,184.556,4.9792],'undef',1.34674,color=Magenta,axes=axes)
axes=[[  0.0,0.0,10.7739],[0.0,10.7739,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_122_0,plain,[122,211.49,4.9792],'6',1.34674,color=Magenta,axes=axes)
labels_all.extend( size_girth_122_0)
center__122_0=[]
axes=[[  0.0,0.0,10.7739],[0.0,10.7739,0.0],[0.0,0.0,0.0]]
cyl_text(center__122_0,plain,[122,305.762,4.9792],'1000_0005',1.34674,color=Magenta,axes=axes)
labels_all.extend( center__122_0)

cluster=[ CYLINDER, 62,359.632,0,122,359.632,0,2,1,0.787497,0,1,0.787497,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.787497,0, SPHERE, 62,359.632,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.787497,0, SPHERE, 122,359.632,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 50,275.821,0,50,324.795,0,2,1,1,0.0597391,1,1,0.0597391]
dendrogram.extend( cluster)
size_girth_50_2_75821=[]
axes=[[  0.0,0.0,2.4487],[0.0,2.4487,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_2_75821,plain,[50,288.064,2.89789],'3',0.306088,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.4487],[0.0,2.4487,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_2_75821,plain,[50,285.004,2.89789],'2.75821',0.306088,color=Magenta,axes=axes)
axes=[[  0.0,0.0,2.4487],[0.0,2.4487,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_2_75821,plain,[50,291.125,2.89789],'13',0.306088,color=Magenta,axes=axes)
labels_all.extend( size_girth_50_2_75821)
center__50_2_75821=[]
axes=[[  0.0,0.0,2.4487],[0.0,2.4487,0.0],[0.0,0.0,0.0]]
cyl_text(center__50_2_75821,plain,[50,312.552,2.89789],'1000_0004',0.306088,color=Magenta,axes=axes)
labels_all.extend( center__50_2_75821)

cluster=[ CYLINDER, 98,144.153,0,98,324.795,0,2,1,0.124026,0,1,0.124026,0]
dendrogram.extend( cluster)
size_girth_98_0=[]
axes=[[  0.0,0.0,9.03209],[0.0,9.03209,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_0,plain,[98,189.314,4.54374],'1',1.12901,color=Magenta,axes=axes)
axes=[[  0.0,0.0,9.03209],[0.0,9.03209,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_0,plain,[98,178.024,4.54374],'undef',1.12901,color=Magenta,axes=axes)
axes=[[  0.0,0.0,9.03209],[0.0,9.03209,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_98_0,plain,[98,200.604,4.54374],'2',1.12901,color=Magenta,axes=axes)
labels_all.extend( size_girth_98_0)
center__98_0=[]
axes=[[  0.0,0.0,9.03209],[0.0,9.03209,0.0],[0.0,0.0,0.0]]
cyl_text(center__98_0,plain,[98,279.635,4.54374],'1000_0001',1.12901,color=Magenta,axes=axes)
labels_all.extend( center__98_0)

cluster=[ CYLINDER, 50,324.795,0,98,324.795,0,2,1,0.825811,0,1,0.825811,0]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.825811,0, SPHERE, 50,324.795,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,0.825811,0, SPHERE, 98,324.795,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 38,156.603,0,38,275.821,0,2,1,1,0.347683,1,1,0.347683]
dendrogram.extend( cluster)
size_girth_38_1_56603=[]
axes=[[  0.0,0.0,5.9609],[0.0,5.9609,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_38_1_56603,plain,[38,186.408,3.77594],'2',0.745112,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.9609],[0.0,5.9609,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_38_1_56603,plain,[38,178.956,3.77594],'1.56603',0.745112,color=Magenta,axes=axes)
axes=[[  0.0,0.0,5.9609],[0.0,5.9609,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_38_1_56603,plain,[38,193.859,3.77594],'12',0.745112,color=Magenta,axes=axes)
labels_all.extend( size_girth_38_1_56603)
center__38_1_56603=[]
axes=[[  0.0,0.0,5.9609],[0.0,5.9609,0.0],[0.0,0.0,0.0]]
cyl_text(center__38_1_56603,plain,[38,246.016,3.77594],'1000_0004',0.745112,color=Magenta,axes=axes)
labels_all.extend( center__38_1_56603)

cluster=[ CYLINDER, 74,144.153,0,74,275.821,0,2,1,0.483852,0,1,0.483852,0]
dendrogram.extend( cluster)
size_girth_74_0=[]
axes=[[  0.0,0.0,6.58338],[0.0,6.58338,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_0,plain,[74,177.07,3.93156],'1',0.822923,color=Magenta,axes=axes)
axes=[[  0.0,0.0,6.58338],[0.0,6.58338,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_0,plain,[74,168.841,3.93156],'undef',0.822923,color=Magenta,axes=axes)
axes=[[  0.0,0.0,6.58338],[0.0,6.58338,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_74_0,plain,[74,185.299,3.93156],'9',0.822923,color=Magenta,axes=axes)
labels_all.extend( size_girth_74_0)
center__74_0=[]
axes=[[  0.0,0.0,6.58338],[0.0,6.58338,0.0],[0.0,0.0,0.0]]
cyl_text(center__74_0,plain,[74,242.904,3.93156],'1000_0008',0.822923,color=Magenta,axes=axes)
labels_all.extend( center__74_0)

cluster=[ CYLINDER, 38,275.821,0,74,275.821,0,2,1,1,0.0597391,1,1,0.0597391]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,1,0.0597391, SPHERE, 38,275.821,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,1,0.0597391, SPHERE, 74,275.821,0,2]
dendrogram.extend( cluster)
cluster=[ CYLINDER, 26,144.153,0,26,156.603,0,2,1,1,0.386061,1,1,0.386061]
dendrogram.extend( cluster)
size_girth_26_0=[]
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_26_0,plain,[26,147.266,2.44134],'1',0.0778106,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_26_0,plain,[26,146.488,2.44134],'undef',0.0778106,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_26_0,plain,[26,148.044,2.44134],'5',0.0778106,color=Magenta,axes=axes)
labels_all.extend( size_girth_26_0)
center__26_0=[]
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(center__26_0,plain,[26,153.491,2.44134],'1000_0004',0.0778106,color=Magenta,axes=axes)
labels_all.extend( center__26_0)

cluster=[ CYLINDER, 50,144.153,0,50,156.603,0,2,1,1,0.309304,1,1,0.309304]
dendrogram.extend( cluster)
size_girth_50_0=[]
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,147.266,2.44134],'1',0.0778106,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,146.488,2.44134],'undef',0.0778106,color=Magenta,axes=axes)
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(size_girth_50_0,plain,[50,148.044,2.44134],'8',0.0778106,color=Magenta,axes=axes)
labels_all.extend( size_girth_50_0)
center__50_0=[]
axes=[[  0.0,0.0,0.622485],[0.0,0.622485,0.0],[0.0,0.0,0.0]]
cyl_text(center__50_0,plain,[50,153.491,2.44134],'1000_0007',0.0778106,color=Magenta,axes=axes)
labels_all.extend( center__50_0)

cluster=[ CYLINDER, 26,156.603,0,50,156.603,0,2,1,1,0.347683,1,1,0.347683]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,1,0.347683, SPHERE, 26,156.603,0,2]
dendrogram.extend( cluster)
cluster=[ COLOR, 1,1,0.347683, SPHERE, 50,156.603,0,2]
dendrogram.extend( cluster)
cmd.load_cgo( dendrogram, 'dendrogram')
cmd.load_cgo( labels_all, 'labels_all')

axes=[[ 19.9195 ,0.0,0.0],[0.0,9.95976,0.0],[0.0,0.0,0.0]]
labels=[]
# girth label color
Black = [0,0,0]

cyl_text(labels,plain,[250.49,649.611,0],'6.546',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,599.812,0],'6.048',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,550.014,0],'5.550',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,500.215,0],'5.052',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,450.416,0],'4.554',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,400.617,0],'4.056',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,350.818,0],'3.558',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,301.02,0],'3.060',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,251.221,0],'2.562',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,201.422,0],'2.064',1.24497,color=Black,axes=axes)
cyl_text(labels,plain,[250.49,151.623,0],'1.566',1.24497,color=Black,axes=axes)
cmd.load_cgo(labels,'labels')

