#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1a53_1thf/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1a53_1thf/combined_0_fused.pdb, fused_model
load example/example_files/output/biol/1a53_1thf/donor1.pdb, donor1_model
load example/example_files/output/biol/1a53_1thf/donor2.pdb, donor2_model
load example/example_files/output/biol/1a53_1thf/donor3.pdb, donor3_model
hide all
cd example/example_files/output/biol/1a53_1thf

#END HEADER
show cartoon, /fused_model//D
show cartoon, /scaffold_model//D
color scaffold_color, /fused_model//D
color scaffold_color, /scaffold_model//D
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//D/6-6
color overlap_color, /donor1_model//A/47-47
color donor_color, /donor1_model//A/48-51
color ignore_color, /scaffold_model//D/7-10
select donor1_site, donor1_site | /fused_model//D/7-10
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/52-64
select donor1_site, donor1_site | /fused_model//D/11-23
color donor_color, donor1_site
color ignore_color, /scaffold_model//D/11-31
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/65-66
color ignore_color, /scaffold_model//D/32-33
select donor1_site, donor1_site | /fused_model//D/24-25
color donor_color, donor1_site
color overlap_color, /scaffold_model//D/34-41
color overlap_color, /donor1_model//A/67-74
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
color overlap_color, /scaffold_model//D/128-130
color overlap_color, /donor2_model//A/156-158
color donor_color, /donor2_model//A/159-159
color ignore_color, /scaffold_model//D/131-131
select donor2_site, donor2_site | /fused_model//D/123-123
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/160-164
select donor2_site, donor2_site | /fused_model//D/124-128
color donor_color, donor2_site
color ignore_color, /scaffold_model//D/132-152
#cutpoint donor->scaffold BEGIN
color donor_color, /donor2_model//A/165-170
color ignore_color, /scaffold_model//D/153-158
select donor2_site, donor2_site | /fused_model//D/129-134
color donor_color, donor2_site
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
color overlap_color, /donor3_model//A/175-176
color donor_color, /donor3_model//A/177-179
color ignore_color, /scaffold_model//D/169-171
select donor3_site, donor3_site | /fused_model//D/145-147
color donor_color, donor3_site
#cutpoint scaffold->donor END
color donor_color, /donor3_model//A/180-193
select donor3_site, donor3_site | /fused_model//D/148-161
color donor_color, donor3_site
color ignore_color, /scaffold_model//D/172-183
#cutpoint donor->scaffold BEGIN
color donor_color, /donor3_model//A/194-194
color ignore_color, /scaffold_model//D/184-184
select donor3_site, donor3_site | /fused_model//D/162-162
color donor_color, donor3_site
color overlap_color, /scaffold_model//D/185-189
color overlap_color, /donor3_model//A/195-199
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
color scaffold_color, /fused_model//D/168-231
color scaffold_color, /scaffold_model//D/190-253
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//D/6/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//D/6/C
pseudoatom peptide_atom_n1l, selection=/fused_model//D/7/N
pseudoatom peptide_atom_car1l, selection=/fused_model//D/7/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//D/6/CA /fused_model//D/6/C
show sticks, /fused_model//D/7/N /fused_model//D/7/CA
color carbon, /fused_model//D/6/C
color nitrogen, /fused_model//D/7/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal1r, selection=/fused_model//D/25/CA
pseudoatom peptide_atom_c1r, selection=/fused_model//D/25/C
pseudoatom peptide_atom_n1r, selection=/fused_model//D/26/N
pseudoatom peptide_atom_car1r, selection=/fused_model//D/26/CA
group atoms1r, peptide_atom_cal1r peptide_atom_c1r peptide_atom_n1r peptide_atom_car1r
group peptide_atoms, atoms1r
show sticks, /fused_model//D/25/CA /fused_model//D/25/C
show sticks, /fused_model//D/26/N /fused_model//D/26/CA
color carbon, /fused_model//D/25/C
color nitrogen, /fused_model//D/26/N
dihedral angle1r, peptide_atom_cal1r, peptide_atom_c1r, peptide_atom_n1r, peptide_atom_car1r
distance dist1r, peptide_atom_c1r, peptide_atom_n1r
group peptide_bad, angle1r
group peptide_bad, dist1r
pseudoatom peptide_atom_cal2l, selection=/fused_model//D/122/CA
pseudoatom peptide_atom_c2l, selection=/fused_model//D/122/C
pseudoatom peptide_atom_n2l, selection=/fused_model//D/123/N
pseudoatom peptide_atom_car2l, selection=/fused_model//D/123/CA
group atoms2l, peptide_atom_cal2l peptide_atom_c2l peptide_atom_n2l peptide_atom_car2l
group peptide_atoms, atoms2l
show sticks, /fused_model//D/122/CA /fused_model//D/122/C
show sticks, /fused_model//D/123/N /fused_model//D/123/CA
color carbon, /fused_model//D/122/C
color nitrogen, /fused_model//D/123/N
dihedral angle2l, peptide_atom_cal2l, peptide_atom_c2l, peptide_atom_n2l, peptide_atom_car2l
distance dist2l, peptide_atom_c2l, peptide_atom_n2l
group peptide_bad, angle2l
group peptide_bad, dist2l
pseudoatom peptide_atom_cal3l, selection=/fused_model//D/144/CA
pseudoatom peptide_atom_c3l, selection=/fused_model//D/144/C
pseudoatom peptide_atom_n3l, selection=/fused_model//D/145/N
pseudoatom peptide_atom_car3l, selection=/fused_model//D/145/CA
group atoms3l, peptide_atom_cal3l peptide_atom_c3l peptide_atom_n3l peptide_atom_car3l
group peptide_atoms, atoms3l
show sticks, /fused_model//D/144/CA /fused_model//D/144/C
show sticks, /fused_model//D/145/N /fused_model//D/145/CA
color carbon, /fused_model//D/144/C
color nitrogen, /fused_model//D/145/N
dihedral angle3l, peptide_atom_cal3l, peptide_atom_c3l, peptide_atom_n3l, peptide_atom_car3l
distance dist3l, peptide_atom_c3l, peptide_atom_n3l
group peptide_bad, angle3l
group peptide_bad, dist3l
pseudoatom peptide_atom_cal3r, selection=/fused_model//D/162/CA
pseudoatom peptide_atom_c3r, selection=/fused_model//D/162/C
pseudoatom peptide_atom_n3r, selection=/fused_model//D/163/N
pseudoatom peptide_atom_car3r, selection=/fused_model//D/163/CA
group atoms3r, peptide_atom_cal3r peptide_atom_c3r peptide_atom_n3r peptide_atom_car3r
group peptide_atoms, atoms3r
show sticks, /fused_model//D/162/CA /fused_model//D/162/C
show sticks, /fused_model//D/163/N /fused_model//D/163/CA
color carbon, /fused_model//D/162/C
color nitrogen, /fused_model//D/163/N
dihedral angle3r, peptide_atom_cal3r, peptide_atom_c3r, peptide_atom_n3r, peptide_atom_car3r
distance dist3r, peptide_atom_c3r, peptide_atom_n3r
group peptide_bad, angle3r
group peptide_bad, dist3r
disable peptide_good
disable peptide_bad
disable peptide_atoms
