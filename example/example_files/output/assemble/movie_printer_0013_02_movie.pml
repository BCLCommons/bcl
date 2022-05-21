python
from pymol import cmd

def color_bcl_model( object_name, min_resi, max_resi):
  """
AUTHOR
  Nils Woetzel
USAGE
  color_bcl_model( object_name, min_resi, max_resi)

  This function colors any protein by the minimum and maximum residue id, according to the rainbow palette.
  This is helpful is multiple structures have to be colored, but the structures have only different residues
  present, but each residue in different structures have to have the same color.
  """

  # process arguments
  try:
    min_resi = int( min_resi)
    max_resi = int( max_resi)
  except ValueError:
    cmd.spectrum( "count", "rainbow", object_name, "byres")
    return
  object_name = str( object_name)

  # set the b factor of each residue to the resi
  for residue in range( min_resi, max_resi):
    cmd.alter( object_name + " and resi " + str( residue), "b=" + str( residue))

  # color according to b factor, which is the actual residue id
  cmd.spectrum( "b", "rainbow", object_name, min_resi, max_resi)
cmd.extend( "color_bcl_model", color_bcl_model)

python end
load example/example_files/output/assemble/movie_printer_0013_02_00123_improved.pdb, model_current
hide everything, model_current
cartoon automatic, model_current
show cartoon, model_current
color_bcl_model( "model_current", min_seq_id, max_seq_id)
zoom model_start
cmd.move( 'x', camera_translate_x)
label label0, "aaclash          0.00 "
label label1, "aadist          -2.38 "
label label2, "aaneigh       -332.16 "
label label3, "aaneigh_ent      0.00 "
label label4, "loop           -33.95 "
label label5, "loop_closure     0.00 "
label label6, "rgyr            73.29 "
label label7, "co_score       -52.57 "
label label8, "sseclash         0.00 "
label label9, "ssepack_fr    -174.96 "
label label10, "strand_fr    -1285.55 "
label label11, "sum          -1808.28 "
set label_color, improved, legend_improved
wait_counter = "%04d" % counter
cmd.png("example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
set label_color, black, legend_improved
delete model_current
load example/example_files/output/assemble/movie_printer_0013_02_final.pdb, model_current
hide lines, model_current
cartoon automatic, model_current
show cartoon, model_current
color_bcl_model( "model_current", min_seq_id, max_seq_id)
zoom model_start
cmd.move( 'x', camera_translate_x)
label label0, "aaclash          0.00 "
label label1, "aadist          -2.38 "
label label2, "aaneigh       -332.16 "
label label3, "aaneigh_ent      0.00 "
label label4, "loop           -33.95 "
label label5, "loop_closure     0.00 "
label label6, "rgyr            73.29 "
label label7, "co_score       -52.57 "
label label8, "sseclash         0.00 "
label label9, "ssepack_fr    -174.96 "
label label10, "strand_fr    -1285.55 "
label label11, "sum          -1808.28 "
wait_counter = "%04d" % counter
cmd.png("example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
