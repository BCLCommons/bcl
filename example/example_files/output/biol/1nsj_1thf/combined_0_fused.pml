#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1thf/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1thf/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1thf/donor1.pdb, donor1_model
load example/example_files/output/biol/1nsj_1thf/donor2.pdb, donor2_model
load example/example_files/output/biol/1nsj_1thf/donor3.pdb, donor3_model
hide all
cd example/example_files/output/biol/1nsj_1thf

#END HEADER
show cartoon, /fused_model//D
show cartoon, /scaffold_model//D
color scaffold_color, /fused_model//D
color scaffold_color, /scaffold_model//D
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//D/7-9
color overlap_color, /donor1_model//A/3-5
color donor_color, /donor1_model//A/6-6
color ignore_color, /scaffold_model//D/10-10
select donor1_site, donor1_site | /fused_model//D/10-10
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/7-14
select donor1_site, donor1_site | /fused_model//D/11-18
color donor_color, donor1_site
color ignore_color, /scaffold_model//D/11-31
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/15-21
color ignore_color, /scaffold_model//D/32-38
select donor1_site, donor1_site | /fused_model//D/19-25
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
color overlap_color, /scaffold_model//D/99-99
color overlap_color, /donor2_model//A/79-79
color donor_color, /donor2_model//A/80-82
color ignore_color, /scaffold_model//D/100-102
select donor2_site, donor2_site | /fused_model//D/87-89
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/83-87
select donor2_site, donor2_site | /fused_model//D/90-94
color donor_color, donor2_site
color ignore_color, /scaffold_model//D/103-109
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//A/88-93
color ignore_color, /scaffold_model//D/110-115
select donor2_site, donor2_site | /fused_model//D/95-100
color donor_color, donor2_site
color overlap_color, /scaffold_model//D/116-118
color overlap_color, /donor2_model//A/94-96
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
color overlap_color, /scaffold_model//D/167-168
color overlap_color, /donor3_model//A/124-125
color donor_color, /donor3_model//A/126-127
color ignore_color, /scaffold_model//D/169-170
select donor3_site, donor3_site | /fused_model//D/154-155
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//A/128-140
select donor3_site, donor3_site | /fused_model//D/156-168
color donor_color, donor3_site
color ignore_color, /scaffold_model//D/171-186
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//A/141-141
color ignore_color, /scaffold_model//D/187-187
select donor3_site, donor3_site | /fused_model//D/169-169
color donor_color, donor3_site
color overlap_color, /scaffold_model//D/188-193
color overlap_color, /donor3_model//A/142-147
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
color scaffold_color, /fused_model//D/176-235
color scaffold_color, /scaffold_model//D/194-253
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//D/9/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//D/9/C
pseudoatom peptide_atom_n1l, selection=/fused_model//D/10/N
pseudoatom peptide_atom_car1l, selection=/fused_model//D/10/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//D/9/CA /fused_model//D/9/C
show sticks, /fused_model//D/10/N /fused_model//D/10/CA
color carbon, /fused_model//D/9/C
color nitrogen, /fused_model//D/10/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal2l, selection=/fused_model//D/86/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//D/86/C
pseudoatom peptide_atom_n2l, selection=/fused_model//D/87/N
pseudoatom peptide_atom_car2l, selection=/fused_model//D/87/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//D/86/CA /fused_model//D/86/C
show sticks, /fused_model//D/87/N /fused_model//D/87/CA
color carbon, /fused_model//D/86/C
color nitrogen, /fused_model//D/87/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal2r, selection=/fused_model//D/100/CA
pseudoatom peptide_atom_c2r, selection=/fused_model//D/100/C
pseudoatom peptide_atom_n2r, selection=/fused_model//D/101/N
pseudoatom peptide_atom_car2r, selection=/fused_model//D/101/CA
group atoms2r, peptide_atom_cal2r peptide_atom_c2r peptide_atom_n2r peptide_atom_car2r
group peptide_atoms, atoms2r
show sticks, /fused_model//D/100/CA /fused_model//D/100/C
show sticks, /fused_model//D/101/N /fused_model//D/101/CA
color carbon, /fused_model//D/100/C
color nitrogen, /fused_model//D/101/N
dihedral angle2r, peptide_atom_cal2r, peptide_atom_c2r, peptide_atom_n2r, peptide_atom_car2r
distance dist2r, peptide_atom_c2r, peptide_atom_n2r
group peptide_good, angle2r
group peptide_good, dist2r
pseudoatom peptide_atom_cal3l, selection=/fused_model//D/153/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//D/153/C
pseudoatom peptide_atom_n3l, selection=/fused_model//D/154/N
pseudoatom peptide_atom_car3l, selection=/fused_model//D/154/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//D/153/CA /fused_model//D/153/C
show sticks, /fused_model//D/154/N /fused_model//D/154/CA
color carbon, /fused_model//D/153/C
color nitrogen, /fused_model//D/154/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//D/169/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//D/169/C
pseudoatom peptide_atom_n3r, selection=/fused_model//D/170/N
pseudoatom peptide_atom_car3r, selection=/fused_model//D/170/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//D/169/CA /fused_model//D/169/C
show sticks, /fused_model//D/170/N /fused_model//D/170/CA
color carbon, /fused_model//D/169/C
color nitrogen, /fused_model//D/170/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
disable peptide_good
disable peptide_bad
disable peptide_atoms
