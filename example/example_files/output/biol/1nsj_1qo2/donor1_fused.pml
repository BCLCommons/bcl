#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1qo2/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1qo2/donor1_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1qo2/donor1.pdb, donor1_model
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
color scaffold_color, /fused_model//A/23-225
color scaffold_color, /scaffold_model//A/39-241
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
disable peptide_good
disable peptide_bad
disable peptide_atoms
