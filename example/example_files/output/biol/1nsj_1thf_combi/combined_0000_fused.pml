#BEGIN HEADER
bg_color white
set_color scaffold_color= [0.95 , 0.78 , 0.00]
set_color overlap_color=  [0.00 , 0.53 , 0.22]
set_color donor_color=    [0.02 , 0.50 , 0.72]
set_color ignore_color=   [1.00 , 0.00 , 0.00]
load example/example_files/output/biol/1nsj_1thf_combi/scaffold.pdb, scaffold_model
load example/example_files/output/biol/1nsj_1thf_combi/combined_0000_fused.pdb, fused_model
load example/example_files/output/biol/1nsj_1thf_combi/donor1.pdb, donor1_model
load example/example_files/output/biol/1nsj_1thf_combi/donor2.pdb, donor2_model
hide all
cd example/example_files/output/biol/1nsj_1thf_combi

#END HEADER
show cartoon, /fused_model//D
show cartoon, /scaffold_model//D
color scaffold_color, /fused_model//D
color scaffold_color, /scaffold_model//D
select donor1_site, donor1_model & scaffold_model
show cartoon, /donor1_model//A/
color ignore_color, /donor1_model//A/
#cutpoint scaffold->donor BEGIN
color donor_color, /donor1_model//A/3-6
color ignore_color, /scaffold_model//D/7-10
select donor1_site, donor1_site | /fused_model//D/7-10
color donor_color, donor1_site
#cutpoint scaffold->donor END
color donor_color, /donor1_model//A/7-14
select donor1_site, donor1_site | /fused_model//D/11-18
color donor_color, donor1_site
color ignore_color, /scaffold_model//D/11-31
#cutpoint donor->scaffold BEGIN
color overlap_color, /scaffold_model//D/32-38
color overlap_color, /donor1_model//A/15-21
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
select donor2_site, donor2_model & scaffold_model
show cartoon, /donor2_model//A/
color ignore_color, /donor2_model//A/
#cutpoint scaffold->donor BEGIN
color donor_color, /donor2_model//A/79-82
color ignore_color, /scaffold_model//D/99-102
select donor2_site, donor2_site | /fused_model//D/86-89
color donor_color, donor2_site
#cutpoint scaffold->donor END
color donor_color, /donor2_model//A/83-87
select donor2_site, donor2_site | /fused_model//D/90-94
color donor_color, donor2_site
color ignore_color, /scaffold_model//D/103-109
#cutpoint donor->scaffold BEGIN
color overlap_color, /scaffold_model//D/110-118
color overlap_color, /donor2_model//A/88-96
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
color scaffold_color, /scaffold_model//D/119-253
group sites, donor*_site
disable sites
group donor_models, donor*_model
disable donor_models
#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' -delay 0 fused_site.png -delay 100 label:'site' animation.gif
enable all
disable peptide_good
disable peptide_bad
disable peptide_atoms
