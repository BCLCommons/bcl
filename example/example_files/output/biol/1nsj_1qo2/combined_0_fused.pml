#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1qo2/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1qo2/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1qo2/donor1.pdb, donor1_model
load example/example_files/output/biol/1nsj_1qo2/donor2.pdb, donor2_model
load example/example_files/output/biol/1nsj_1qo2/donor3.pdb, donor3_model
hide all
cd example/example_files/output/biol/1nsj_1qo2

#END HEADER
show cartoon, /fused_model//A
show cartoon, /scaffold_model//A
color scaffold_color, /fused_model//A
color scaffold_color, /scaffold_model//A
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/4-4
color overlap_color, /donor1_model//A/3-3
color donor_color, /donor1_model//A/4-6
color ignore_color, /scaffold_model//A/5-7
select donor1_site, donor1_site | /fused_model//A/5-7
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/7-13
select donor1_site, donor1_site | /fused_model//A/8-14
color donor_color, donor1_site
color ignore_color, /scaffold_model//A/8-30
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/14-21
color ignore_color, /scaffold_model//A/31-38
select donor1_site, donor1_site | /fused_model//A/15-22
color donor_color, donor1_site
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor1_site
rotate y, 90
#ray
png fused_site_donor1_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor1_model.png, dpi=300
disable all
enable donor1_model
#ray
png donor1_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor1_model.png, dpi=300
disable all
enable fused_model
orient donor1_site
rotate y, 90
#ray
png fused_site_donor1_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor1_model.png, dpi=300
disable all
enable donor1_model
#ray
png donor1_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor1_model.png, dpi=300
select donor2_site, donor2_model & scaffold_model
show cartoon, /donor2_model//A/
color ignore_color, /donor2_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/98-98
color overlap_color, /donor2_model//A/79-79
color donor_color, /donor2_model//A/80-82
color ignore_color, /scaffold_model//A/99-101
select donor2_site, donor2_site | /fused_model//A/83-85
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/83-88
select donor2_site, donor2_site | /fused_model//A/86-91
color donor_color, donor2_site
color ignore_color, /scaffold_model//A/102-109
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//A/89-92
color ignore_color, /scaffold_model//A/110-113
select donor2_site, donor2_site | /fused_model//A/92-95
color donor_color, donor2_site
color overlap_color, /scaffold_model//A/114-117
color overlap_color, /donor2_model//A/93-96
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor2_site
rotate y, 90
#ray
png fused_site_donor2_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor2_model.png, dpi=300
disable all
enable donor2_model
#ray
png donor2_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor2_model.png, dpi=300
disable all
enable fused_model
orient donor2_site
rotate y, 90
#ray
png fused_site_donor2_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor2_model.png, dpi=300
disable all
enable donor2_model
#ray
png donor2_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor2_model.png, dpi=300
select donor3_site, donor3_model & scaffold_model
show cartoon, /donor3_model//A/
color ignore_color, /donor3_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/160-160
color overlap_color, /donor3_model//A/124-124
color donor_color, /donor3_model//A/125-127
color ignore_color, /scaffold_model//A/161-163
select donor3_site, donor3_site | /fused_model//A/143-145
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//A/128-140
select donor3_site, donor3_site | /fused_model//A/146-158
color donor_color, donor3_site
color ignore_color, /scaffold_model//A/164-180
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//A/141-147
color ignore_color, /scaffold_model//A/181-187
select donor3_site, donor3_site | /fused_model//A/159-165
color donor_color, donor3_site
color overlap_color, /scaffold_model//A/188-188
color overlap_color, /donor3_model//A/148-148
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor3_site
rotate y, 90
#ray
png fused_site_donor3_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor3_model.png, dpi=300
disable all
enable donor3_model
#ray
png donor3_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor3_model.png, dpi=300
disable all
enable fused_model
orient donor3_site
rotate y, 90
#ray
png fused_site_donor3_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor3_model.png, dpi=300
disable all
enable donor3_model
#ray
png donor3_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor3_model.png, dpi=300
color scaffold_color, /fused_model//A/167-219
color scaffold_color, /scaffold_model//A/189-241
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//A/4/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//A/4/C
pseudoatom peptide_atom_n1l, selection=/fused_model//A/5/N
pseudoatom peptide_atom_car1l, selection=/fused_model//A/5/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//A/4/CA /fused_model//A/4/C
show sticks, /fused_model//A/5/N /fused_model//A/5/CA
color carbon, /fused_model//A/4/C
color nitrogen, /fused_model//A/5/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal2l, selection=/fused_model//A/82/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//A/82/C
pseudoatom peptide_atom_n2l, selection=/fused_model//A/83/N
pseudoatom peptide_atom_car2l, selection=/fused_model//A/83/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//A/82/CA /fused_model//A/82/C
show sticks, /fused_model//A/83/N /fused_model//A/83/CA
color carbon, /fused_model//A/82/C
color nitrogen, /fused_model//A/83/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal2r, selection=/fused_model//A/95/CA
pseudoatom peptide_atom_c2r, selection=/fused_model//A/95/C
pseudoatom peptide_atom_n2r, selection=/fused_model//A/96/N
pseudoatom peptide_atom_car2r, selection=/fused_model//A/96/CA
group atoms2r, peptide_atom_cal2r peptide_atom_c2r peptide_atom_n2r peptide_atom_car2r
group peptide_atoms, atoms2r
show sticks, /fused_model//A/95/CA /fused_model//A/95/C
show sticks, /fused_model//A/96/N /fused_model//A/96/CA
color carbon, /fused_model//A/95/C
color nitrogen, /fused_model//A/96/N
dihedral angle2r, peptide_atom_cal2r, peptide_atom_c2r, peptide_atom_n2r, peptide_atom_car2r
distance dist2r, peptide_atom_c2r, peptide_atom_n2r
group peptide_good, angle2r
group peptide_good, dist2r
pseudoatom peptide_atom_cal3l, selection=/fused_model//A/142/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//A/142/C
pseudoatom peptide_atom_n3l, selection=/fused_model//A/143/N
pseudoatom peptide_atom_car3l, selection=/fused_model//A/143/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//A/142/CA /fused_model//A/142/C
show sticks, /fused_model//A/143/N /fused_model//A/143/CA
color carbon, /fused_model//A/142/C
color nitrogen, /fused_model//A/143/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//A/165/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//A/165/C
pseudoatom peptide_atom_n3r, selection=/fused_model//A/166/N
pseudoatom peptide_atom_car3r, selection=/fused_model//A/166/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//A/165/CA /fused_model//A/165/C
show sticks, /fused_model//A/166/N /fused_model//A/166/CA
color carbon, /fused_model//A/165/C
color nitrogen, /fused_model//A/166/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
disable peptide_good
disable peptide_bad
disable peptide_atoms
