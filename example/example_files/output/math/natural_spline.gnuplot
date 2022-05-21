# BCL generated heatmap
set terminal png enhanced transparent font "Arial,12" size 1080,800
set output "heatmap.png"
set encoding iso_8859_1
set view map
set title "Natural spline"
unset key

set xrange [ -0.5 : 35.5 ]
set xtics rotate by 0 (" -180 " 0, " -170 " 1, " -160 " 2, " -150 " 3, " -140 " 4, " -130 " 5, " -120 " 6, " -110 " 7, " -100 " 8, " -90 " 9, " -80 " 10, " -70 " 11, " -60 " 12, " -50 " 13, " -40 " 14, " -30 " 15, " -20 " 16, " -10 " 17, " 0 " 18, " 10 " 19, " 20 " 20, " 30 " 21, " 40 " 22, " 50 " 23, " 60 " 24, " 70 " 25, " 80 " 26, " 90 " 27, " 100 " 28, " 110 " 29, " 120 " 30, " 130 " 31, " 140 " 32, " 150 " 33, " 160 " 34, " 170 " 35)
set yrange [ -0.5 : 1.5 ]
set noytics
set cbrange [ * : * ]
#set cbtics 1
#set format cb "%3.1f"

set palette rgbformulae 22, 13, -31
splot '-' using 1:2:3 with image
# number x values 36
# number y values 2
0	0 	26
1	0 	3
2	0 	1
3	0 	2
4	0 	1
5	0 	3
6	0 	6
7	0 	3
8	0 	8
9	0 	2
10	0 	7
11	0 	8
12	0 	3
13	0 	4
14	0 	2
15	0 	1
16	0 	2
17	0 	5
18	0 	30
19	0 	0
20	0 	2
21	0 	4
22	0 	6
23	0 	3
24	0 	4
25	0 	3
26	0 	3
27	0 	4
28	0 	11
29	0 	5
30	0 	8
31	0 	5
32	0 	2
33	0 	0
34	0 	2
35	0 	2

0	1 	26
1	1 	3
2	1 	1
3	1 	2
4	1 	1
5	1 	3
6	1 	6
7	1 	3
8	1 	8
9	1 	2
10	1 	7
11	1 	8
12	1 	3
13	1 	4
14	1 	2
15	1 	1
16	1 	2
17	1 	5
18	1 	30
19	1 	0
20	1 	2
21	1 	4
22	1 	6
23	1 	3
24	1 	4
25	1 	3
26	1 	3
27	1 	4
28	1 	11
29	1 	5
30	1 	8
31	1 	5
32	1 	2
33	1 	0
34	1 	2
35	1 	2

e
