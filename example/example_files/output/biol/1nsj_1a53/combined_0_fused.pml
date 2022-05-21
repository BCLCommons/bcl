#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1a53/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1a53/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1a53/donor1.pdb, donor1_model
load example/example_files/output/biol/1nsj_1a53/donor2.pdb, donor2_model
load example/example_files/output/biol/1nsj_1a53/donor3.pdb, donor3_model
hide all
cd example/example_files/output/biol/1nsj_1a53

#END HEADER
show cartoon, /fused_model//A
show cartoon, /scaffold_model//A
color scaffold_color, /fused_model//A
color scaffold_color, /scaffold_model//A
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/47-47
color overlap_color, /donor1_model//A/4-4
color donor_color, /donor1_model//A/5-6
color ignore_color, /scaffold_model//A/48-49
select donor1_site, donor1_site | /fused_model//A/48-49
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/7-11
select donor1_site, donor1_site | /fused_model//A/50-54
color donor_color, donor1_site
color ignore_color, /scaffold_model//A/50-65
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/12-16
color ignore_color, /scaffold_model//A/66-70
select donor1_site, donor1_site | /fused_model//A/55-59
color donor_color, donor1_site
color overlap_color, /scaffold_model//A/71-74
color overlap_color, /donor1_model//A/17-20
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
color overlap_color, /scaffold_model//A/130-131
color overlap_color, /donor2_model//A/79-80
color donor_color, /donor2_model//A/81-81
color ignore_color, /scaffold_model//A/132-132
select donor2_site, donor2_site | /fused_model//A/121-121
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/82-87
select donor2_site, donor2_site | /fused_model//A/122-127
color donor_color, donor2_site
color ignore_color, /scaffold_model//A/133-137
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//A/88-88
color ignore_color, /scaffold_model//A/138-138
select donor2_site, donor2_site | /fused_model//A/128-128
color donor_color, donor2_site
color overlap_color, /scaffold_model//A/139-146
color overlap_color, /donor2_model//A/89-96
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
color overlap_color, /scaffold_model//A/175-177
color overlap_color, /donor3_model//A/122-124
color donor_color, /donor3_model//A/125-126
color ignore_color, /scaffold_model//A/178-179
select donor3_site, donor3_site | /fused_model//A/168-169
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//A/127-140
select donor3_site, donor3_site | /fused_model//A/170-183
color donor_color, donor3_site
color ignore_color, /scaffold_model//A/180-189
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//A/141-144
color ignore_color, /scaffold_model//A/190-193
select donor3_site, donor3_site | /fused_model//A/184-187
color donor_color, donor3_site
color overlap_color, /scaffold_model//A/194-199
color overlap_color, /donor3_model//A/145-150
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
color scaffold_color, /fused_model//A/194-241
color scaffold_color, /scaffold_model//A/200-247
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//A/47/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//A/47/C
pseudoatom peptide_atom_n1l, selection=/fused_model//A/48/N
pseudoatom peptide_atom_car1l, selection=/fused_model//A/48/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//A/47/CA /fused_model//A/47/C
show sticks, /fused_model//A/48/N /fused_model//A/48/CA
color carbon, /fused_model//A/47/C
color nitrogen, /fused_model//A/48/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal1r, selection=/fused_model//A/59/CA
pseudoatom peptide_atom_c1r, selection=/fused_model//A/59/C
pseudoatom peptide_atom_n1r, selection=/fused_model//A/60/N
pseudoatom peptide_atom_car1r, selection=/fused_model//A/60/CA
group atoms1r, peptide_atom_cal1r peptide_atom_c1r peptide_atom_n1r peptide_atom_car1r
group peptide_atoms, atoms1r
show sticks, /fused_model//A/59/CA /fused_model//A/59/C
show sticks, /fused_model//A/60/N /fused_model//A/60/CA
color carbon, /fused_model//A/59/C
color nitrogen, /fused_model//A/60/N
dihedral angle1r, peptide_atom_cal1r, peptide_atom_c1r, peptide_atom_n1r, peptide_atom_car1r
distance dist1r, peptide_atom_c1r, peptide_atom_n1r
group peptide_good, angle1r
group peptide_good, dist1r
pseudoatom peptide_atom_cal2l, selection=/fused_model//A/120/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//A/120/C
pseudoatom peptide_atom_n2l, selection=/fused_model//A/121/N
pseudoatom peptide_atom_car2l, selection=/fused_model//A/121/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//A/120/CA /fused_model//A/120/C
show sticks, /fused_model//A/121/N /fused_model//A/121/CA
color carbon, /fused_model//A/120/C
color nitrogen, /fused_model//A/121/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal2r, selection=/fused_model//A/128/CA
pseudoatom peptide_atom_c2r, selection=/fused_model//A/128/C
pseudoatom peptide_atom_n2r, selection=/fused_model//A/129/N
pseudoatom peptide_atom_car2r, selection=/fused_model//A/129/CA
group atoms2r, peptide_atom_cal2r peptide_atom_c2r peptide_atom_n2r peptide_atom_car2r
group peptide_atoms, atoms2r
show sticks, /fused_model//A/128/CA /fused_model//A/128/C
show sticks, /fused_model//A/129/N /fused_model//A/129/CA
color carbon, /fused_model//A/128/C
color nitrogen, /fused_model//A/129/N
dihedral angle2r, peptide_atom_cal2r, peptide_atom_c2r, peptide_atom_n2r, peptide_atom_car2r
distance dist2r, peptide_atom_c2r, peptide_atom_n2r
group peptide_bad, angle2r
group peptide_bad, dist2r
pseudoatom peptide_atom_cal3l, selection=/fused_model//A/167/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//A/167/C
pseudoatom peptide_atom_n3l, selection=/fused_model//A/168/N
pseudoatom peptide_atom_car3l, selection=/fused_model//A/168/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//A/167/CA /fused_model//A/167/C
show sticks, /fused_model//A/168/N /fused_model//A/168/CA
color carbon, /fused_model//A/167/C
color nitrogen, /fused_model//A/168/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//A/187/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//A/187/C
pseudoatom peptide_atom_n3r, selection=/fused_model//A/188/N
pseudoatom peptide_atom_car3r, selection=/fused_model//A/188/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//A/187/CA /fused_model//A/187/C
show sticks, /fused_model//A/188/N /fused_model//A/188/CA
color carbon, /fused_model//A/187/C
color nitrogen, /fused_model//A/188/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
disable peptide_good
disable peptide_bad
disable peptide_atoms
