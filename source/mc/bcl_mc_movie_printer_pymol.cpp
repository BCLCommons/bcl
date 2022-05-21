// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "mc/bcl_mc_movie_printer_pymol.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb.h"
#include "score/bcl_score_protein_model.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
  //////////
  // data //
  //////////

    const size_t MoviePrinterPymol::s_WaitTimes[] =
    {
      5, //e_Improved
      2, //e_Accepted
      1, //e_Skipped
      1, //e_Rejected
      0  //s_NumberStepStatus
    };

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoviePrinterPymol::s_Instance( GetObjectInstances().AddInstance( new MoviePrinterPymol()));

    //! pymol function to color by residue id
    const std::string MoviePrinterPymol::s_ColorFunction
    (
      "from pymol import cmd\n\n"
      "def color_bcl_model( object_name, min_resi, max_resi):\n"
      "  \"\"\"\n"
      "AUTHOR\n"
      "  Nils Woetzel\n"
      "USAGE\n"
      "  color_bcl_model( object_name, min_resi, max_resi)\n\n"
      "  This function colors any protein by the minimum and maximum residue id, according to the rainbow palette.\n"
      "  This is helpful is multiple structures have to be colored, but the structures have only different residues\n"
      "  present, but each residue in different structures have to have the same color.\n"
      "  \"\"\"\n\n"
      "  # process arguments\n"
      "  try:\n"
      "    min_resi = int( min_resi)\n"
      "    max_resi = int( max_resi)\n"
      "  except ValueError:\n"
      "    cmd.spectrum( \"count\", \"rainbow\", object_name, \"byres\")\n"
      "    return\n"
      "  object_name = str( object_name)\n\n"
      "  # set the b factor of each residue to the resi\n"
      "  for residue in range( min_resi, max_resi):\n"
      "    cmd.alter( object_name + \" and resi \" + str( residue), \"b=\" + str( residue))\n\n"
      "  # color according to b factor, which is the actual residue id\n"
      "  cmd.spectrum( \"b\", \"rainbow\", object_name, min_resi, max_resi)\n"
      "cmd.extend( \"color_bcl_model\", color_bcl_model)\n\n"
    );

    const size_t MoviePrinterPymol::s_TableWidth( 240);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoviePrinterPymol::MoviePrinterPymol() :
      m_RayTracing( true),
      m_Width( 720),
      m_Height( 480)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoviePrinterInterface
    MoviePrinterPymol *MoviePrinterPymol::Clone() const
    {
      return new MoviePrinterPymol( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoviePrinterPymol::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default extension for movie script file
    //! @return extension for script files
    const std::string &MoviePrinterPymol::GetScriptFileExtension() const
    {
      // extension
      static const std::string s_extension( "pml");

      // return
      return s_extension;
    }

    //! @brief initialize with all required information
    //! @param PREFIX prefix as absolute path
    //! @param COL_NAMES the column names to be considered
    //! @param ROW_NAMES the row name to be considered
    //! @param WIDTH width in pixel
    //! @param HEIGHT height in pixel
    //! @param ENABLE_RAY true or false
    void MoviePrinterPymol::Initialize
    (
      const std::string &PREFIX,
      const storage::Vector< std::string> &COL_NAMES,
      const storage::Vector< std::string> &ROW_NAMES,
      const size_t WIDTH, const size_t HEIGHT,
      const bool ENABLE_RAY
    )
    {
      // check for absolute path of prefix
      if( !PREFIX.empty() && !io::File::IsAbsolutePath( PREFIX))
      {
        BCL_MessageVrb
        (
          PREFIX + " does not appear to be an absolute path, but it is recommended to use absolute path!"
        )
      }

      m_FirstScript.clear();
      m_MovieScript.clear();
      m_Prefix = PREFIX;
      m_HasStartFrame = false;

      m_LabelColNames = COL_NAMES;
      m_LabelRowNames = ROW_NAMES;

      // find the length of the widest row name
      size_t label_width( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_LabelRowNames.Begin()), itr_end( m_LabelRowNames.End());
        itr != itr_end;
        ++itr
      )
      {
        if
        (
          itr->size() > label_width &&
          *itr != score::ProteinModel::GetTypeName( score::ProteinModel::e_Structure) &&
          *itr != score::ProteinModel::GetTypeName( score::ProteinModel::e_Sequence)
        )
        {
          label_width = itr->size();
        }
      }

      m_LabelFormats.clear();
      m_LabelFormats.push_back( util::Format().W( label_width).L());
      m_LabelFormats.push_back( util::Format().W( 9).R().FFP( 2));

      BCL_Assert( WIDTH > s_TableWidth, "Movie width must be larger than " + util::Format()( s_TableWidth) + " pixels.");
      m_Width = WIDTH;
      m_Height = HEIGHT;

      m_RayTracing = ENABLE_RAY;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset
    void MoviePrinterPymol::Reset( const std::string &PREFIX)
    {
      m_Prefix = PREFIX;
      m_MovieScript.clear();
    }

    //! @brief add start frame
    //! @param FILE_NAME the filename of the starting frame
    //! @param KEEP_FRAME_UP is start frame permanent (like a density map)
    void MoviePrinterPymol::SetStartFrame
    (
      const std::string &FILE_NAME,
      const bool KEEP_FRAME_UP
    )
    {
      // check for absolute path of FILE_NAME
      if( !FILE_NAME.empty() && !io::File::IsAbsolutePath( FILE_NAME))
      {
        BCL_MessageVrb
        (
          FILE_NAME + " does not appear to be a absolute path, but it is recommended to use absolute path!"
        )
      }

      // empty the first script
      m_FirstScript.clear();

      // add color function
      m_FirstScript += "#sequence length, so that incomplete models can be colored rainbow relative to full sequence\n";
      m_FirstScript += "min_seq_id=\"unknown\"\n";
      m_FirstScript += "max_seq_id=\"unknown\"\n";

      // add counter
      m_FirstScript += "counter=0\n";

      // define window size
      m_FirstScript += "viewport " + util::Format()( m_Width - s_TableWidth) + ", " + util::Format()( m_Height) + '\n';

      // define background color
      m_FirstScript += "bg_color white\n";

      // load the given file
      if( !FILE_NAME.empty())
      {
        m_HasStartFrame = true;
        // if mrc map
        if( io::File::GetLastExtension( FILE_NAME) == "mrc" || io::File::GetLastExtension( FILE_NAME) == "ccp4")
        {
          m_FirstScript += "set ray_shadows, off\n";
          m_FirstScript += "unset normalize_ccp4_maps\n";
          m_FirstScript += "'mrc' : CCP4 map\n";
          // load volume
          m_FirstScript += "load " + FILE_NAME + ", model_start\n";
          m_FirstScript += "isomesh surface2, model_start, 1\n";
          m_FirstScript += "set mesh_width, 0.5\n";
        }

        // if pdb file
        else if( io::File::GetLastExtension( FILE_NAME) == pdb::GetDefaultFileExtension())
        {
          m_FirstScript += "load " + FILE_NAME + ", model_start\n"; // load structure
          m_FirstScript += "hide everything, model_start\n";             // hides atoms and bonds
          m_FirstScript += "cartoon automatic, model_start\n";      // shows it as ribbon
          m_FirstScript += "show cartoon, model_start\n";           // shows it as ribbon
          m_FirstScript += "color_bcl_model( \"model_start\", min_seq_id, max_seq_id)\n"; // rainbow the model
          m_FirstScript += "set cartoon_transparency, 0.75, model_start\n"; // transparent model
        }

        // rotate
        m_FirstScript += "rotate x, 90, model_start\n";

        // make frame
        if( m_RayTracing)
        {
          m_FirstScript += "ray " + util::Format()( m_Width - s_TableWidth) + ", " + util::Format()( m_Height) + '\n';
        }

        const std::string png_name( m_Prefix + "\" + \"%04d\" % counter + \".png");
        m_FirstScript += "cmd.png(\"" + png_name + "\")\n";

        // unload given file if frame should nto stay up
        if( !KEEP_FRAME_UP)
        {
          // unload structure
          m_FirstScript += "set ray_shadows, on\n";
          m_FirstScript += "hide cartoon, model_start";
        }

        m_FirstScript += TableToLabel( storage::Table< double>(), png_name);
        m_FirstScript += WriteWaterMark( png_name);
        m_FirstScript += "counter += 1\n";
      }
    }

    //! @brief add frame to movie
    //! @param FILE_NAME the filename of the pdb the model was written to
    //! @param STEP_STATUS the mc step status
    //! @param LABEL add label to frame
    void MoviePrinterPymol::AddFrame
    (
      const std::string &FILE_NAME,
      const opti::StepStatus &STEP_STATUS,
      const storage::Table< double> &LABEL
    )
    {
      // check for absolute path of FILE_NAME
      if( !FILE_NAME.empty() && !io::File::IsAbsolutePath( FILE_NAME))
      {
        BCL_MessageVrb
        (
          FILE_NAME + " does not appear to be a absolute path, but it is recommended to use absolute path!"
        );
      }

      // create script for that frame
      std::string frame_script;

      // add to movie script
      frame_script += "load " + FILE_NAME + ", model_current\n"; // opens file
      frame_script += "hide everything, model_current\n";        // hides atoms and bonds
      frame_script += "cartoon automatic, model_current\n";      // shows it as ribbon
      frame_script += "show cartoon, model_current\n";           // shows it as ribbon
      frame_script += "color_bcl_model( \"model_current\", min_seq_id, max_seq_id)\n"; // colors the model
      if( m_HasStartFrame)
      {
        frame_script += "zoom model_start\n";           // zoom to start model
      }
      else
      {
        frame_script += "zoom model_current\n";           // zoom to current model
      }

      // label
      const std::string png_name( m_Prefix + "\" + wait_counter + \".png");

      // rotate
      frame_script += "rotate x, 90, model_current\n";

      // make frame
      if( m_RayTracing)
      {
        frame_script += "ray " + util::Format()( m_Width - s_TableWidth) + ", " + util::Format()( m_Height) + '\n';
      }

      frame_script += "wait_counter = \"%04d\" % counter\n";
      frame_script += "cmd.png(\"" + png_name + "\")\n";
      frame_script += TableToLabel( LABEL, png_name);
      frame_script += WriteLegend( STEP_STATUS, png_name);
      frame_script += WriteWaterMark( png_name);
      frame_script += "counter += 1\n";

      for( size_t count( 1); count < s_WaitTimes[ STEP_STATUS]; ++count)
      {
        frame_script += "cmd.system(\"ln -s " + m_Prefix + "\" + wait_counter + \".png \" + \"" + m_Prefix + "\" + \"%04d\" % counter + \".png\")\n";
        frame_script += "counter += 1\n";
      }

      // delete model
      frame_script += "delete model_current\n";

      // add frame script to MovieScript
      m_MovieScript += frame_script;
    }

    //! @brief add final frame to movie
    //! @param FILE_NAME the filename of the pdb the model was written to
    //! @param LABEL add label to frame
    void MoviePrinterPymol::AddFinalFrame
    (
      const std::string &FILE_NAME,
      const storage::Table< double> &LABEL
    )
    {
      // check for absolute path of FILE_NAME
      if( !FILE_NAME.empty() && !io::File::IsAbsolutePath( FILE_NAME))
      {
        BCL_MessageVrb
        (
          FILE_NAME + " does not appear to be a absolute path, but it is recommended to use absolute path!"
        );
      }

      // create script for that frame
      std::string frame_script;

      // add to movie script
      frame_script += "load " + FILE_NAME + ", model_current\n"; // opens file
      frame_script += "hide lines, model_current\n";             // hides atoms and bonds
      frame_script += "cartoon automatic, model_current\n";      // shows it as ribbon
      frame_script += "show cartoon, model_current\n";           // shows it as ribbon
      frame_script += "color_bcl_model( \"model_current\", min_seq_id, max_seq_id)\n"; // colors the model
      if( m_HasStartFrame)
      {
        frame_script += "zoom model_start\n";           // zoom to start model
      }
      else
      {
        frame_script += "zoom model_current\n";           // zoom to current model
      }

      // label
      const std::string png_name( m_Prefix + "\" + wait_counter + \".png");

      // rotate
      frame_script += "rotate x, 90, model_current\n";

      // make frame
      if( m_RayTracing)
      {
        frame_script += "ray " + util::Format()( m_Width - s_TableWidth) + ", " + util::Format()( m_Height) + '\n';
      }

      frame_script += "wait_counter = \"%04d\" % counter\n";
      frame_script += "cmd.png(\"" + png_name + "\")\n";
      frame_script += TableToLabel( LABEL, png_name);
      frame_script += WriteWaterMark( png_name);
      frame_script += "counter += 1\n";

      for( size_t count( 1); count < GetFinalWaitTime(); ++count)
      {
        frame_script += "cmd.system(\"ln -s " + m_Prefix + "\" + wait_counter + \".png \" + \"" + m_Prefix + "\" + \"%04d\" % counter + \".png\")\n";
        frame_script += "counter += 1\n";
      }

      // add frame script to MovieScript
      m_MovieScript += frame_script;
    }

    //! @brief set the watermark
    //! @param WATERMARK the watermark
    void MoviePrinterPymol::SetWaterMark( const std::string &WATERMARK)
    {
      m_WaterMark = WATERMARK;
    }

    //! @brief create ffmpeg command line
    //! @return commandline to be called with ffmpeg to generate the movie from the generated pngs
    std::string MoviePrinterPymol::FFMpegCommandLine() const
    {
      const std::string command_line
      (
        "ffmpeg" // command
        " -r 25" // frame rate
        " -i " + m_Prefix + "%04d.png" // pattern for files
        " -y" // overwriting existing file
        " -f mp4" // format
        " -s " + util::Format()( m_Width) + "x" + util::Format()( m_Height) + // size
        " -qscale 1"       // quality
        " " + m_Prefix + "movie.mp4" // output name
      );

      // end
      return command_line;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the movie script to a file
    //! @param OSTREAM
    std::ostream &MoviePrinterPymol::WriteScript
    (
      std::ostream &OSTREAM
    ) const
    {
      // write coloring script
      OSTREAM << "python\n" << s_ColorFunction << "python end\n";

      // write first portion and all frames
      OSTREAM << m_FirstScript;

      // write all frames
      OSTREAM << m_MovieScript;

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoviePrinterPymol::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoviePrinterPymol::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert table to label
    //! @param TABLE
    //! @param PNG_NAME
    //! @return string
    std::string MoviePrinterPymol::TableToLabel
    (
      const storage::Table< double> &TABLE,
      const std::string &PNG_NAME
    ) const
    {
      // write table to string
      std::stringstream tmp;
      TABLE.WriteSubTable( tmp, m_LabelColNames, m_LabelRowNames, m_LabelFormats, false);

      // add space for labels
      std::string label_script
      (
        "cmd.system(\"convert " + PNG_NAME + " -gravity west -extent " +
        util::Format()( m_Width) + "x" + util::Format()( m_Height) + " " + PNG_NAME + "\")\n"
      );

      // split label to modify
      const storage::Vector< std::string> labels( util::SplitString( tmp.str(), "\n"));

      size_t y_pos( 50);
      // iterate over all labels
      for
      (
        storage::Vector< std::string>::const_iterator itr( labels.Begin()), itr_end( labels.End());
        itr != itr_end; ++itr
      )
      {
        // check for undefined to replace with empty string
        std::string this_label( *itr);
        const size_t nan_pos( this_label.find( "nan"));
        if( nan_pos != std::string::npos)
        {
          this_label.replace( nan_pos, 3, "   ");
        }

        // set the string of the label
        label_script += "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
                        "-draw 'text " + util::Format()( m_Width - s_TableWidth) + "," + util::Format()( y_pos) +
                        " \\\"" + this_label + "\\\"' " + PNG_NAME + "\")\n";

        // move down
        y_pos += 15;
      }

      // end
      return label_script;
    }

    //! @brief Write the step status legend
    //! @param STEP_STATUS given step status
    //! @param PNG_NAME
    //! @return the step status legend
    std::string MoviePrinterPymol::WriteLegend( const opti::StepStatus &STEP_STATUS, const std::string &PNG_NAME) const
    {
      const size_t legend_height( 75);
      std::string legend;

      // build string
      legend += "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
             ( STEP_STATUS == opti::e_Improved ? "-fill green " : "") + "-draw 'text " + util::Format()( m_Width - s_TableWidth / 2) +
             "," + util::Format()( m_Height - legend_height) + " \\\"improved\\\"' " + PNG_NAME + "\")\n";
      legend += "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
             ( STEP_STATUS == opti::e_Accepted ? "-fill yellow " : "") + "-draw 'text " + util::Format()( m_Width - s_TableWidth / 2) +
             "," + util::Format()( m_Height - legend_height + 15) + " \\\"accepted\\\"' " + PNG_NAME + "\")\n";
      legend += "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
             ( STEP_STATUS == opti::e_Rejected ? "-fill red " : "") + "-draw 'text " + util::Format()( m_Width - s_TableWidth / 2) +
             "," + util::Format()( m_Height - legend_height + 30) + " \\\"rejected\\\"' " + PNG_NAME + "\")\n";
      legend += "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
             "-draw 'text " + util::Format()( m_Width - s_TableWidth / 2) +
             "," + util::Format()( m_Height - legend_height + 45) + " \\\"skipped\\\"' " + PNG_NAME + "\")\n";

      // end
      return legend;
    }

    //! @brief Write the watermark
    //! @param PNG_NAME
    //! @return the watermark
    std::string MoviePrinterPymol::WriteWaterMark( const std::string &PNG_NAME)
    {
      return "cmd.system(\"convert " + PNG_NAME + " -gravity northwest -pointsize 12 -font courier " +
          "-draw 'text 10," + util::Format()( m_Height - 10) + " \\\"" + m_WaterMark + "\\\"' " + PNG_NAME + "\")\n";
    }

  } // namespace mc
} // namespace bcl
