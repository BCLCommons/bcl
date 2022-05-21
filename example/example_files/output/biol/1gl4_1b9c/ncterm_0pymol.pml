#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1gl4_1b9c/scaffold.pdb, scaffold_model

#END HEADER
load example/example_files/output/biol/1gl4_1b9c/ncterm_0fused.pdb, fused_model
load example/example_files/output/biol/1gl4_1b9c/nterm_donor.pdb, donor_nterm_model
load example/example_files/output/biol/1gl4_1b9c/cterm_donor.pdb, donor_cterm_model
hide all
cd example/example_files/output/biol/1gl4_1b9c
show cartoon, /fused_model//A
show cartoon, /scaffold_model//A
color scaffold_color, /fused_model//A
color scaffold_color, /scaffold_model//A
select donor_nterm_site, donor_nterm_model & scaffold_model
show cartoon, /donor_nterm_model//A/
color ignore_color, /donor_nterm_model//A/
color donor_color, /donor_nterm_model//A/1-10
color ignore_color, /scaffold_model//A/1-59
#cutpoint donor->scaffold BEGIN
color donor_color, /donor_nterm_model//A/11-12
color ignore_color, /scaffold_model//A/60-61
select donor_nterm_site, donor_nterm_site | /fused_model//A/11-12
color donor_color, donor_nterm_site
color overlap_color, /scaffold_model//A/62-62
color overlap_color, /donor_nterm_model//A/13-13
#cutpoint donor->scaffold END
disable all
enable fused_model
orient donor_nterm_site
rotate y, 90
#ray
png fused_site_donor_nterm_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor_nterm_model.png, dpi=300
disable all
enable donor_nterm_model
#ray
png donor_nterm_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor_nterm_model.png, dpi=300
select donor_cterm_site, donor_cterm_model & scaffold_model
show cartoon, /donor_cterm_model//A/
color ignore_color, /donor_cterm_model//A/
#cutpoint scaffold->donor BEGIN
color overlap_color, /scaffold_model//A/266-266
color overlap_color, /donor_cterm_model//A/266-266
color donor_color, /donor_cterm_model//A/267-276
color ignore_color, /scaffold_model//A/267-276
select donor_cterm_site, donor_cterm_site | /fused_model//A/218-227
color donor_color, donor_cterm_site
#cutpoint scaffold->donor END
disable all
enable fused_model
orient donor_cterm_site
rotate y, 90
#ray
png fused_site_donor_cterm_model.png, dpi=300
orient fused_model
#ray
png fused_model_donor_cterm_model.png, dpi=300
disable all
enable donor_cterm_model
#ray
png donor_cterm_model.png, dpi=300
disable all
enable scaffold_model
#ray
png scaffold_model_donor_cterm_model.png, dpi=300
pseudoatom peptide_atom_cal0r, selection=/fused_model//A/12/CA
pseudoatom peptide_atom_c0r, selection=/fused_model//A/12/C
pseudoatom peptide_atom_n0r, selection=/fused_model//A/13/N
pseudoatom peptide_atom_car0r, selection=/fused_model//A/13/CA
group atoms0r, peptide_atom_cal0r peptide_atom_c0r peptide_atom_n0r peptide_atom_car0r
group peptide_atoms, atoms0r
show sticks, /fused_model//A/12/CA /fused_model//A/12/C
show sticks, /fused_model//A/13/N /fused_model//A/13/CA
color carbon, /fused_model//A/12/C
color nitrogen, /fused_model//A/13/N
dihedral angle0r, peptide_atom_cal0r, peptide_atom_c0r, peptide_atom_n0r, peptide_atom_car0r
distance dist0r, peptide_atom_c0r, peptide_atom_n0r
group peptide_bad, angle0r
group peptide_bad, dist0r
pseudoatom peptide_atom_cal0l, selection=/fused_model//A/217/CA
pseudoatom peptide_atom_c0l, selection=/fused_model//A/217/C
pseudoatom peptide_atom_n0l, selection=/fused_model//A/218/N
pseudoatom peptide_atom_car0l, selection=/fused_model//A/218/CA
group atoms0l, peptide_atom_cal0l peptide_atom_c0l peptide_atom_n0l peptide_atom_car0l
group peptide_atoms, atoms0l
show sticks, /fused_model//A/217/CA /fused_model//A/217/C
show sticks, /fused_model//A/218/N /fused_model//A/218/CA
color carbon, /fused_model//A/217/C
color nitrogen, /fused_model//A/218/N
dihedral angle0l, peptide_atom_cal0l, peptide_atom_c0l, peptide_atom_n0l, peptide_atom_car0l
distance dist0l, peptide_atom_c0l, peptide_atom_n0l
group peptide_good, angle0l
group peptide_good, dist0l
enable all
disable peptide_good
disable peptide_bad
disable peptide_atoms
