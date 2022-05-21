#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1thf_1qo2/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1thf_1qo2/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1thf_1qo2/donor1.pdb, donor1_model
load example/example_files/output/biol/1thf_1qo2/donor2.pdb, donor2_model
load example/example_files/output/biol/1thf_1qo2/donor3.pdb, donor3_model
hide all
cd example/example_files/output/biol/1thf_1qo2

#END HEADER
show cartoon, /fused_model//A
show cartoon, /scaffold_model//A
color scaffold_color, /fused_model//A
color scaffold_color, /scaffold_model//A
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//D/
color ignore_color, /donor1_model//D/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/3-6
color overlap_color, /donor1_model//D/6-9
color donor_color, /donor1_model//D/10-12
color ignore_color, /scaffold_model//A/7-9
select donor1_site, donor1_site | /fused_model//A/7-9
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//D/13-31
select donor1_site, donor1_site | /fused_model//A/10-28
color donor_color, donor1_site
color ignore_color, /scaffold_model//A/10-31
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//D/32-37
color ignore_color, /scaffold_model//A/32-37
select donor1_site, donor1_site | /fused_model//A/29-34
color donor_color, donor1_site
color overlap_color, /scaffold_model//A/38-42
color overlap_color, /donor1_model//D/38-42
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
show cartoon, /donor2_model//D/
color ignore_color, /donor2_model//D/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/122-125
color overlap_color, /donor2_model//D/125-128
color donor_color, /donor2_model//D/129-131
color ignore_color, /scaffold_model//A/126-128
select donor2_site, donor2_site | /fused_model//A/123-125
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//D/132-152
select donor2_site, donor2_site | /fused_model//A/126-146
color donor_color, donor2_site
color ignore_color, /scaffold_model//A/129-145
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//D/153-155
color ignore_color, /scaffold_model//A/146-148
select donor2_site, donor2_site | /fused_model//A/147-149
color donor_color, donor2_site
color overlap_color, /scaffold_model//A/149-155
color overlap_color, /donor2_model//D/156-162
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
show cartoon, /donor3_model//D/
color ignore_color, /donor3_model//D/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/160-163
color overlap_color, /donor3_model//D/167-170
color donor_color, /donor3_model//D/171-172
color ignore_color, /scaffold_model//A/164-165
select donor3_site, donor3_site | /fused_model//A/165-166
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//D/173-183
select donor3_site, donor3_site | /fused_model//A/167-177
color donor_color, donor3_site
color ignore_color, /scaffold_model//A/166-176
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//D/184-189
color ignore_color, /scaffold_model//A/177-182
select donor3_site, donor3_site | /fused_model//A/178-183
color donor_color, donor3_site
color overlap_color, /scaffold_model//A/183-186
color overlap_color, /donor3_model//D/190-193
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
color scaffold_color, /fused_model//A/188-242
color scaffold_color, /scaffold_model//A/187-241
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//A/6/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//A/6/C
pseudoatom peptide_atom_n1l, selection=/fused_model//A/7/N
pseudoatom peptide_atom_car1l, selection=/fused_model//A/7/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//A/6/CA /fused_model//A/6/C
show sticks, /fused_model//A/7/N /fused_model//A/7/CA
color carbon, /fused_model//A/6/C
color nitrogen, /fused_model//A/7/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal1r, selection=/fused_model//A/34/CA
pseudoatom peptide_atom_c1r, selection=/fused_model//A/34/C
pseudoatom peptide_atom_n1r, selection=/fused_model//A/35/N
pseudoatom peptide_atom_car1r, selection=/fused_model//A/35/CA
group atoms1r, peptide_atom_cal1r peptide_atom_c1r peptide_atom_n1r peptide_atom_car1r
group peptide_atoms, atoms1r
show sticks, /fused_model//A/34/CA /fused_model//A/34/C
show sticks, /fused_model//A/35/N /fused_model//A/35/CA
color carbon, /fused_model//A/34/C
color nitrogen, /fused_model//A/35/N
dihedral angle1r, peptide_atom_cal1r, peptide_atom_c1r, peptide_atom_n1r, peptide_atom_car1r
distance dist1r, peptide_atom_c1r, peptide_atom_n1r
group peptide_bad, angle1r
group peptide_bad, dist1r
pseudoatom peptide_atom_cal2l, selection=/fused_model//A/122/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//A/122/C
pseudoatom peptide_atom_n2l, selection=/fused_model//A/123/N
pseudoatom peptide_atom_car2l, selection=/fused_model//A/123/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//A/122/CA /fused_model//A/122/C
show sticks, /fused_model//A/123/N /fused_model//A/123/CA
color carbon, /fused_model//A/122/C
color nitrogen, /fused_model//A/123/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal2r, selection=/fused_model//A/149/CA
pseudoatom peptide_atom_c2r, selection=/fused_model//A/149/C
pseudoatom peptide_atom_n2r, selection=/fused_model//A/150/N
pseudoatom peptide_atom_car2r, selection=/fused_model//A/150/CA
group atoms2r, peptide_atom_cal2r peptide_atom_c2r peptide_atom_n2r peptide_atom_car2r
group peptide_atoms, atoms2r
show sticks, /fused_model//A/149/CA /fused_model//A/149/C
show sticks, /fused_model//A/150/N /fused_model//A/150/CA
color carbon, /fused_model//A/149/C
color nitrogen, /fused_model//A/150/N
dihedral angle2r, peptide_atom_cal2r, peptide_atom_c2r, peptide_atom_n2r, peptide_atom_car2r
distance dist2r, peptide_atom_c2r, peptide_atom_n2r
group peptide_bad, angle2r
group peptide_bad, dist2r
pseudoatom peptide_atom_cal3l, selection=/fused_model//A/164/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//A/164/C
pseudoatom peptide_atom_n3l, selection=/fused_model//A/165/N
pseudoatom peptide_atom_car3l, selection=/fused_model//A/165/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//A/164/CA /fused_model//A/164/C
show sticks, /fused_model//A/165/N /fused_model//A/165/CA
color carbon, /fused_model//A/164/C
color nitrogen, /fused_model//A/165/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//A/183/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//A/183/C
pseudoatom peptide_atom_n3r, selection=/fused_model//A/184/N
pseudoatom peptide_atom_car3r, selection=/fused_model//A/184/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//A/183/CA /fused_model//A/183/C
show sticks, /fused_model//A/184/N /fused_model//A/184/CA
color carbon, /fused_model//A/183/C
color nitrogen, /fused_model//A/184/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
disable peptide_good
disable peptide_bad
disable peptide_atoms
