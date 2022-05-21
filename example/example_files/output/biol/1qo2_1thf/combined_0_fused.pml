#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1qo2_1thf/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1qo2_1thf/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1qo2_1thf/donor1.pdb, donor1_model
load example/example_files/output/biol/1qo2_1thf/donor2.pdb, donor2_model
load example/example_files/output/biol/1qo2_1thf/donor3.pdb, donor3_model
load example/example_files/output/biol/1qo2_1thf/donor4.pdb, donor4_model
load example/example_files/output/biol/1qo2_1thf/donor5.pdb, donor5_model
hide all
cd example/example_files/output/biol/1qo2_1thf

#END HEADER
show cartoon, /fused_model//D
show cartoon, /scaffold_model//D
color scaffold_color, /fused_model//D
color scaffold_color, /scaffold_model//D
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//D/6-9
color overlap_color, /donor1_model//A/3-6
color donor_color, /donor1_model//A/7-9
color ignore_color, /scaffold_model//D/10-12
select donor1_site, donor1_site | /fused_model//D/10-12
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/10-31
select donor1_site, donor1_site | /fused_model//D/13-34
color donor_color, donor1_site
color ignore_color, /scaffold_model//D/13-31
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/32-34
color ignore_color, /scaffold_model//D/32-34
select donor1_site, donor1_site | /fused_model//D/35-37
color donor_color, donor1_site
color overlap_color, /scaffold_model//D/35-42
color overlap_color, /donor1_model//A/35-42
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
color overlap_color, /scaffold_model//D/47-48
color overlap_color, /donor2_model//A/47-48
color donor_color, /donor2_model//A/49-51
color ignore_color, /scaffold_model//D/49-51
select donor2_site, donor2_site | /fused_model//D/52-54
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/52-61
select donor2_site, donor2_site | /fused_model//D/55-64
color donor_color, donor2_site
color ignore_color, /scaffold_model//D/52-57
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//A/62-63
color ignore_color, /scaffold_model//D/58-59
select donor2_site, donor2_site | /fused_model//D/65-66
color donor_color, donor2_site
color overlap_color, /scaffold_model//D/60-65
color overlap_color, /donor2_model//A/64-69
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
color overlap_color, /scaffold_model//D/77-79
color overlap_color, /donor3_model//A/76-78
color donor_color, /donor3_model//A/79-79
color ignore_color, /scaffold_model//D/80-80
select donor3_site, donor3_site | /fused_model//D/87-87
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//A/80-84
select donor3_site, donor3_site | /fused_model//D/88-92
color donor_color, donor3_site
color ignore_color, /scaffold_model//D/81-85
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//A/85-90
color ignore_color, /scaffold_model//D/86-91
select donor3_site, donor3_site | /fused_model//D/93-98
color donor_color, donor3_site
color overlap_color, /scaffold_model//D/92-94
color overlap_color, /donor3_model//A/91-93
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
select donor4_site, donor4_model & scaffold_model
show cartoon, /donor4_model//A/
color ignore_color, /donor4_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//D/125-130
color overlap_color, /donor4_model//A/122-127
color donor_color, /donor4_model//A/128-128
color ignore_color, /scaffold_model//D/131-131
select donor4_site, donor4_site | /fused_model//D/138-138
color donor_color, donor4_site
#cutpoint scaffold->donor END
color donor_color, /donor4_model//A/129-145
select donor4_site, donor4_site | /fused_model//D/139-155
color donor_color, donor4_site
color ignore_color, /scaffold_model//D/132-152
#cutpoint donor->scaffold BEGIN
color donor_color, /donor4_model//A/146-150
color ignore_color, /scaffold_model//D/153-157
select donor4_site, donor4_site | /fused_model//D/156-160
color donor_color, donor4_site
color overlap_color, /scaffold_model//D/158-162
color overlap_color, /donor4_model//A/151-155
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor4_site
rotate y, 90
#ray
png fused_site_donor4_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor4_model.png, dpi=300
disable all
enable donor4_model
#ray
png donor4_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor4_model.png, dpi=300
disable all
enable fused_model
orient donor4_site
rotate y, 90
#ray
png fused_site_donor4_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor4_model.png, dpi=300
disable all
enable donor4_model
#ray
png donor4_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor4_model.png, dpi=300
select donor5_site, donor5_model & scaffold_model
show cartoon, /donor5_model//A/
color ignore_color, /donor5_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//D/167-171
color overlap_color, /donor5_model//A/160-164
color donor_color, /donor5_model//A/165-165
color ignore_color, /scaffold_model//D/172-172
select donor5_site, donor5_site | /fused_model//D/175-175
color donor_color, donor5_site
#cutpoint scaffold->donor END
color donor_color, /donor5_model//A/166-176
select donor5_site, donor5_site | /fused_model//D/176-186
color donor_color, donor5_site
color ignore_color, /scaffold_model//D/173-183
#cutpoint donor->scaffold BEGIN
color donor_color, /donor5_model//A/177-182
color ignore_color, /scaffold_model//D/184-189
select donor5_site, donor5_site | /fused_model//D/187-192
color donor_color, donor5_site
color overlap_color, /scaffold_model//D/190-193
color overlap_color, /donor5_model//A/183-186
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor5_site
rotate y, 90
#ray
png fused_site_donor5_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor5_model.png, dpi=300
disable all
enable donor5_model
#ray
png donor5_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor5_model.png, dpi=300
disable all
enable fused_model
orient donor5_site
rotate y, 90
#ray
png fused_site_donor5_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor5_model.png, dpi=300
disable all
enable donor5_model
#ray
png donor5_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor5_model.png, dpi=300
color scaffold_color, /fused_model//D/197-256
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
pseudoatom peptide_atom_cal1r, selection=/fused_model//D/37/CA
pseudoatom peptide_atom_c1r, selection=/fused_model//D/37/C
pseudoatom peptide_atom_n1r, selection=/fused_model//D/38/N
pseudoatom peptide_atom_car1r, selection=/fused_model//D/38/CA
group atoms1r, peptide_atom_cal1r peptide_atom_c1r peptide_atom_n1r peptide_atom_car1r
group peptide_atoms, atoms1r
show sticks, /fused_model//D/37/CA /fused_model//D/37/C
show sticks, /fused_model//D/38/N /fused_model//D/38/CA
color carbon, /fused_model//D/37/C
color nitrogen, /fused_model//D/38/N
dihedral angle1r, peptide_atom_cal1r, peptide_atom_c1r, peptide_atom_n1r, peptide_atom_car1r
distance dist1r, peptide_atom_c1r, peptide_atom_n1r
group peptide_good, angle1r
group peptide_good, dist1r
pseudoatom peptide_atom_cal2l, selection=/fused_model//D/51/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//D/51/C
pseudoatom peptide_atom_n2l, selection=/fused_model//D/52/N
pseudoatom peptide_atom_car2l, selection=/fused_model//D/52/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//D/51/CA /fused_model//D/51/C
show sticks, /fused_model//D/52/N /fused_model//D/52/CA
color carbon, /fused_model//D/51/C
color nitrogen, /fused_model//D/52/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal2r, selection=/fused_model//D/66/CA
pseudoatom peptide_atom_c2r, selection=/fused_model//D/66/C
pseudoatom peptide_atom_n2r, selection=/fused_model//D/67/N
pseudoatom peptide_atom_car2r, selection=/fused_model//D/67/CA
group atoms2r, peptide_atom_cal2r peptide_atom_c2r peptide_atom_n2r peptide_atom_car2r
group peptide_atoms, atoms2r
show sticks, /fused_model//D/66/CA /fused_model//D/66/C
show sticks, /fused_model//D/67/N /fused_model//D/67/CA
color carbon, /fused_model//D/66/C
color nitrogen, /fused_model//D/67/N
dihedral angle2r, peptide_atom_cal2r, peptide_atom_c2r, peptide_atom_n2r, peptide_atom_car2r
distance dist2r, peptide_atom_c2r, peptide_atom_n2r
group peptide_bad, angle2r
group peptide_bad, dist2r
pseudoatom peptide_atom_cal3l, selection=/fused_model//D/86/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//D/86/C
pseudoatom peptide_atom_n3l, selection=/fused_model//D/87/N
pseudoatom peptide_atom_car3l, selection=/fused_model//D/87/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//D/86/CA /fused_model//D/86/C
show sticks, /fused_model//D/87/N /fused_model//D/87/CA
color carbon, /fused_model//D/86/C
color nitrogen, /fused_model//D/87/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//D/98/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//D/98/C
pseudoatom peptide_atom_n3r, selection=/fused_model//D/99/N
pseudoatom peptide_atom_car3r, selection=/fused_model//D/99/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//D/98/CA /fused_model//D/98/C
show sticks, /fused_model//D/99/N /fused_model//D/99/CA
color carbon, /fused_model//D/98/C
color nitrogen, /fused_model//D/99/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
pseudoatom peptide_atom_cal4l, selection=/fused_model//D/137/CA
pseudoatom peptide_atom_c4l, selection=/fused_model//D/137/C
pseudoatom peptide_atom_n4l, selection=/fused_model//D/138/N
pseudoatom peptide_atom_car4l, selection=/fused_model//D/138/CA
group atoms4l, peptide_atom_cal4l peptide_atom_c4l peptide_atom_n4l peptide_atom_car4l
group peptide_atoms, atoms4l
show sticks, /fused_model//D/137/CA /fused_model//D/137/C
show sticks, /fused_model//D/138/N /fused_model//D/138/CA
color carbon, /fused_model//D/137/C
color nitrogen, /fused_model//D/138/N
dihedral angle4l, peptide_atom_cal4l, peptide_atom_c4l, peptide_atom_n4l, peptide_atom_car4l
distance dist4l, peptide_atom_c4l, peptide_atom_n4l
group peptide_bad, angle4l
group peptide_bad, dist4l
pseudoatom peptide_atom_cal4r, selection=/fused_model//D/160/CA
pseudoatom peptide_atom_c4r, selection=/fused_model//D/160/C
pseudoatom peptide_atom_n4r, selection=/fused_model//D/161/N
pseudoatom peptide_atom_car4r, selection=/fused_model//D/161/CA
group atoms4r, peptide_atom_cal4r peptide_atom_c4r peptide_atom_n4r peptide_atom_car4r
group peptide_atoms, atoms4r
show sticks, /fused_model//D/160/CA /fused_model//D/160/C
show sticks, /fused_model//D/161/N /fused_model//D/161/CA
color carbon, /fused_model//D/160/C
color nitrogen, /fused_model//D/161/N
dihedral angle4r, peptide_atom_cal4r, peptide_atom_c4r, peptide_atom_n4r, peptide_atom_car4r
distance dist4r, peptide_atom_c4r, peptide_atom_n4r
group peptide_bad, angle4r
group peptide_bad, dist4r
pseudoatom peptide_atom_cal5l, selection=/fused_model//D/174/CA
pseudoatom peptide_atom_c5l, selection=/fused_model//D/174/C
pseudoatom peptide_atom_n5l, selection=/fused_model//D/175/N
pseudoatom peptide_atom_car5l, selection=/fused_model//D/175/CA
group atoms5l, peptide_atom_cal5l peptide_atom_c5l peptide_atom_n5l peptide_atom_car5l
group peptide_atoms, atoms5l
show sticks, /fused_model//D/174/CA /fused_model//D/174/C
show sticks, /fused_model//D/175/N /fused_model//D/175/CA
color carbon, /fused_model//D/174/C
color nitrogen, /fused_model//D/175/N
dihedral angle5l, peptide_atom_cal5l, peptide_atom_c5l, peptide_atom_n5l, peptide_atom_car5l
distance dist5l, peptide_atom_c5l, peptide_atom_n5l
group peptide_good, angle5l
group peptide_good, dist5l
pseudoatom peptide_atom_cal5r, selection=/fused_model//D/192/CA
pseudoatom peptide_atom_c5r, selection=/fused_model//D/192/C
pseudoatom peptide_atom_n5r, selection=/fused_model//D/193/N
pseudoatom peptide_atom_car5r, selection=/fused_model//D/193/CA
group atoms5r, peptide_atom_cal5r peptide_atom_c5r peptide_atom_n5r peptide_atom_car5r
group peptide_atoms, atoms5r
show sticks, /fused_model//D/192/CA /fused_model//D/192/C
show sticks, /fused_model//D/193/N /fused_model//D/193/CA
color carbon, /fused_model//D/192/C
color nitrogen, /fused_model//D/193/N
dihedral angle5r, peptide_atom_cal5r, peptide_atom_c5r, peptide_atom_n5r, peptide_atom_car5r
distance dist5r, peptide_atom_c5r, peptide_atom_n5r
group peptide_good, angle5r
group peptide_good, dist5r
disable peptide_good
disable peptide_bad
disable peptide_atoms
