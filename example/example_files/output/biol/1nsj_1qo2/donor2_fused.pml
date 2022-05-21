#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1qo2/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1qo2/donor2_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1qo2/donor2.pdb, donor1_model
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
color overlap_color, /scaffold_model//A/98-98
color overlap_color, /donor1_model//A/79-79
color donor_color, /donor1_model//A/80-82
color ignore_color, /scaffold_model//A/99-101
select donor1_site, donor1_site | /fused_model//A/99-101
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/83-88
select donor1_site, donor1_site | /fused_model//A/102-107
color donor_color, donor1_site
color ignore_color, /scaffold_model//A/102-109
#cutpoint donor->scaffold BEGIN
color donor_color, /donor1_model//A/89-92
color ignore_color, /scaffold_model//A/110-113
select donor1_site, donor1_site | /fused_model//A/108-111
color donor_color, donor1_site
color overlap_color, /scaffold_model//A/114-117
color overlap_color, /donor1_model//A/93-96
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
color scaffold_color, /fused_model//A/116-239
color scaffold_color, /scaffold_model//A/118-241
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
pseudoatom peptide_atom_cal1l, selection=/fused_model//A/98/CA
pseudoatom peptide_atom_c1l, selection=/fused_model//A/98/C
pseudoatom peptide_atom_n1l, selection=/fused_model//A/99/N
pseudoatom peptide_atom_car1l, selection=/fused_model//A/99/CA
group atoms1l, peptide_atom_cal1l peptide_atom_c1l peptide_atom_n1l peptide_atom_car1l
group peptide_atoms, atoms1l
show sticks, /fused_model//A/98/CA /fused_model//A/98/C
show sticks, /fused_model//A/99/N /fused_model//A/99/CA
color carbon, /fused_model//A/98/C
color nitrogen, /fused_model//A/99/N
dihedral angle1l, peptide_atom_cal1l, peptide_atom_c1l, peptide_atom_n1l, peptide_atom_car1l
distance dist1l, peptide_atom_c1l, peptide_atom_n1l
group peptide_bad, angle1l
group peptide_bad, dist1l
pseudoatom peptide_atom_cal1r, selection=/fused_model//A/111/CA
pseudoatom peptide_atom_c1r, selection=/fused_model//A/111/C
pseudoatom peptide_atom_n1r, selection=/fused_model//A/112/N
pseudoatom peptide_atom_car1r, selection=/fused_model//A/112/CA
group atoms1r, peptide_atom_cal1r peptide_atom_c1r peptide_atom_n1r peptide_atom_car1r
group peptide_atoms, atoms1r
show sticks, /fused_model//A/111/CA /fused_model//A/111/C
show sticks, /fused_model//A/112/N /fused_model//A/112/CA
color carbon, /fused_model//A/111/C
color nitrogen, /fused_model//A/112/N
dihedral angle1r, peptide_atom_cal1r, peptide_atom_c1r, peptide_atom_n1r, peptide_atom_car1r
distance dist1r, peptide_atom_c1r, peptide_atom_n1r
group peptide_good, angle1r
group peptide_good, dist1r
disable peptide_good
disable peptide_bad
disable peptide_atoms
