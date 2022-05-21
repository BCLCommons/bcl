# BCL generated heatmap
set terminal png enhanced transparent font "Arial,12" size 1080,800
set output "example/example_files/output/biol/1nsj_1thf_combi/donor1_score_heatmap.png"
set encoding iso_8859_1
set view map
set title "1THF->1NSJ-> 1THF"
unset key

set xlabel "donor crossover 1NSJ"
set xrange [ -0.5 : 3.5 ]
set xtics rotate by 0 ("    6 " 0, "    7 " 1, "    8 " 2, "    9 " 3)
set ylabel "scaffold crossover 1THF"
set yrange [ -0.5 : 6.5 ]
set ytics ("   14 " 0, "   15 " 1, "   16 " 2, "   17 " 3, "   18 " 4, "   19 " 5, "   20 " 6)
set cblabel "per residue score"
set cbrange [ * : * ]
#set cbtics 1
#set format cb "%3.1f"

set palette rgbformulae 22, 13, -31
#set palette rgbformulae 3, 11, 6 # green-red-violet
#set palette rgbformulae 30,31,32 # color printable on gray (black-blue-violet-yellow-white)
#set palette rgbformulae 33,13,10 # rainbow (blue-green-yellow-red)
set object rect from 2.5,6.5 to 3.5,7.5  front fillcolor rgb "black" fs empty border 0 linewidth 0.9
splot '-' using 1:2:3 with image
# number x values 4
# number y values 7
0	0 	-4.56654
1	0 	-4.61907
2	0 	-4.57085
3	0 	-4.67009

0	1 	-4.65024
1	1 	-4.70225
2	1 	-4.65415
3	1 	-4.75351

0	2 	-4.64271
1	2 	-4.69481
2	2 	-4.64677
3	2 	-4.74609

0	3 	-4.63046
1	3 	-4.68186
2	3 	-4.63405
3	3 	-4.73312

0	4 	-4.58878
1	4 	-4.64079
2	4 	-4.58958
3	4 	-4.68945

0	5 	-4.5718
1	5 	-4.6233
2	5 	-4.57138
3	5 	-4.67105

0	6 	-4.58644
1	6 	-4.63812
2	6 	-4.58641
3	6 	-4.68612

e
