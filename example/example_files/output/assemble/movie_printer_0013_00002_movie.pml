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
load example/example_files/output/assemble/movie_printer_0013_00002_00123_improved.pdb, model_current
hide everything, model_current
cartoon automatic, model_current
show cartoon, model_current
color_bcl_model( "model_current", min_seq_id, max_seq_id)
zoom model_current
rotate x, 90, model_current
wait_counter = "%04d" % counter
cmd.png("example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity west -extent 720x480 example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -fill green -draw 'text 600,405 \"improved\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -draw 'text 600,420 \"accepted\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -draw 'text 600,435 \"rejected\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -draw 'text 600,450 \"skipped\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -draw 'text 10,470 \"\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
cmd.system("ln -s example/example_files/output/assemble/movie_printer_" + wait_counter + ".png " + "example/example_files/output/assemble/movie_printer_" + "%04d" % counter + ".png")
counter += 1
delete model_current
load example/example_files/output/assemble/movie_printer_0013_00002_final.pdb, model_current
hide lines, model_current
cartoon automatic, model_current
show cartoon, model_current
color_bcl_model( "model_current", min_seq_id, max_seq_id)
zoom model_current
rotate x, 90, model_current
wait_counter = "%04d" % counter
cmd.png("example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity west -extent 720x480 example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
cmd.system("convert example/example_files/output/assemble/movie_printer_" + wait_counter + ".png -gravity northwest -pointsize 12 -font courier -draw 'text 10,470 \"\"' example/example_files/output/assemble/movie_printer_" + wait_counter + ".png")
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
