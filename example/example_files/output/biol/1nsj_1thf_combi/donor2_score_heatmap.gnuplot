# BCL generated heatmap
set terminal png enhanced transparent font "Arial,12" size 1080,800
set output "example/example_files/output/biol/1nsj_1thf_combi/donor2_score_heatmap.png"
set encoding iso_8859_1
set view map
set title "1THF->1NSJ-> 1THF"
unset key

set xlabel "donor crossover 1NSJ"
set xrange [ -0.5 : 3.5 ]
set xtics rotate by 0 ("   98 " 0, "   99 " 1, "  100 " 2, "  101 " 3)
set ylabel "scaffold crossover 1THF"
set yrange [ -0.5 : 8.5 ]
set ytics ("   87 " 0, "   88 " 1, "   89 " 2, "   90 " 3, "   91 " 4, "   92 " 5, "   93 " 6, "   94 " 7, "   95 " 8)
set cblabel "per residue score"
set cbrange [ * : * ]
#set cbtics 1
#set format cb "%3.1f"

set palette rgbformulae 22, 13, -31
#set palette rgbformulae 3, 11, 6 # green-red-violet
#set palette rgbformulae 30,31,32 # color printable on gray (black-blue-violet-yellow-white)
#set palette rgbformulae 33,13,10 # rainbow (blue-green-yellow-red)
set object rect from 0.5,5.5 to 1.5,6.5  front fillcolor rgb "black" fs empty border 0 linewidth 0.9
splot '-' using 1:2:3 with image
# number x values 4
# number y values 9
0	0 	-4.63324
1	0 	-4.62854
2	0 	-4.63351
3	0 	-4.63897

0	1 	-4.59483
1	1 	-4.59048
2	1 	-4.59529
3	1 	-4.60073

0	2 	-4.60029
1	2 	-4.59579
2	2 	-4.60058
3	2 	-4.60613

0	3 	-4.60032
1	3 	-4.59591
2	3 	-4.60064
3	3 	-4.60628

0	4 	-4.58575
1	4 	-4.58095
2	4 	-4.58582
3	4 	-4.59135

0	5 	-4.59292
1	5 	-4.5882
2	5 	-4.59292
3	5 	-4.59881

0	6 	-4.59467
1	6 	-4.58987
2	6 	-4.59453
3	6 	-4.60047

0	7 	-4.5794
1	7 	-4.57501
2	7 	-4.57911
3	7 	-4.58501

0	8 	-4.57777
1	8 	-4.57308
2	8 	-4.5775
3	8 	-4.58286

e
