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
#include "mc/bcl_mc.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_movie_printer_chimera.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
  //////////
  // data //
  //////////

    const size_t MoviePrinterChimera::s_WaitTimes[] =
    {
      5, //e_Improved
      1, //e_Accepted
      1, //e_Skipped
      1  //e_Rejected
    };

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoviePrinterChimera::s_Instance( GetObjectInstances().AddInstance( new MoviePrinterChimera()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoviePrinterChimera::MoviePrinterChimera() :
      m_RayTracing( false),
      m_Width( 720),
      m_Height( 480)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoviePrinterInterface
    MoviePrinterChimera *MoviePrinterChimera::Clone() const
    {
      return new MoviePrinterChimera( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoviePrinterChimera::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default extension for movie script file
    //! @return extension for script files
    const std::string &MoviePrinterChimera::GetScriptFileExtension() const
    {
      // extension
      static const std::string s_extension( "cmd");

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
    void MoviePrinterChimera::Initialize
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

      m_Path       = io::File::SplitToPathAndFileName( PREFIX).First();
      m_FilePrefix = io::File::SplitToPathAndFileName( PREFIX).Second();

      m_LabelColNames = COL_NAMES;
      m_LabelRowNames = ROW_NAMES;

      m_Width = WIDTH;
      m_Height = HEIGHT;

      m_RayTracing = ENABLE_RAY;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset
    void MoviePrinterChimera::Reset( const std::string &PREFIX)
    {
      m_FilePrefix = PREFIX;
      m_FirstScript.clear();
      m_MovieScript.clear();
    }

    //! @brief add start frame
    //! @param FILE_NAME the filename of the starting frame
    //! @param KEEP_FRAME_UP is start frame permanent (like a density map)
    void MoviePrinterChimera::SetStartFrame
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

      // define window size
      m_FirstScript += "windowsize " + util::Format()( m_Width) + " " + util::Format()( m_Height) + "\n";

      // define background color
      m_FirstScript += "set bg_color white\n";

      // define colors for each step kind
      m_FirstScript += "colordef " + GetStepStatusName( opti::e_Accepted) + " 0.0 0.0 1.0 1.0\n"; // blue
      m_FirstScript += "colordef " + GetStepStatusName( opti::e_Rejected) + " 1.0 0.0 0.0 1.0\n"; // red
      m_FirstScript += "colordef " + GetStepStatusName( opti::e_Improved) + " 0.0 1.0 0.0 1.0\n"; // green
      m_FirstScript += "colordef greytrans 0.5 0.5 0.5 0.5\n"; // grey and transparent

      // labels
      m_FirstScript += "2dlabels create 1\n";

      // load the given file
      if( !FILE_NAME.empty())
      {
        // load structure
        m_FirstScript += "open 2 " + FILE_NAME + "\n";

        // if mrc map
        if( io::File::GetLastExtension( FILE_NAME) == "mrc")
        {
          m_FirstScript += "volume #2 style surface\n";
          m_FirstScript += "volume #2 transparency 0.5\n";
          m_FirstScript += "volume #2 step 1\n";
        }

        // if pdb file
        if( io::File::GetLastExtension( FILE_NAME) == pdb::GetDefaultFileExtension())
        {
          m_FirstScript += "~show #2\n"; // hides atoms and bonds
          m_FirstScript += "color greytrans #2\n"; // colors the model
          m_FirstScript += "ribbon #2\n"; // shows it as ribbon
          m_FirstScript += "ribrepr supersmooth #2\n"; // supersmooth ribbon
        }
      }

      // legend
      m_FirstScript += "2dlabels create 3 xpos 0.05 ypos 0.25 size 16 color " + opti::GetStepStatusName( opti::e_Accepted) + " text \"" + GetStepStatusName( opti::e_Accepted) + "\"\n";
      m_FirstScript += "2dlabels create 4 xpos 0.05 ypos 0.20 size 16 color " + opti::GetStepStatusName( opti::e_Rejected) + " text \"" + GetStepStatusName( opti::e_Rejected) + "\"\n";
      m_FirstScript += "2dlabels create 5 xpos 0.05 ypos 0.15 size 16 color " + opti::GetStepStatusName( opti::e_Improved) + " text \"" + GetStepStatusName( opti::e_Improved) + "\"\n";

      // make frame
      m_FirstScript += "movie record raytrace " + std::string( m_RayTracing ? "true" : "false") + " fformat png pattern " + m_FilePrefix + "* directory " + m_Path + "\n";
      m_FirstScript += "movie stop\n";

      // unload given file if frame should nto stay up
      if( !FILE_NAME.empty() && !KEEP_FRAME_UP)
      {
        // unload structure
        m_FirstScript += "close 2";
      }
    }

    //! @brief add frame to movie
    //! @param FILE_NAME the filename of the pdb the model was written to
    //! @param STEP_STATUS the mc step status
    //! @param LABEL add label to frame
    void MoviePrinterChimera::AddFrame
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
      frame_script += "open 1 " + FILE_NAME + "\n"; // opens file
      frame_script += "~show #1\n"; // hides atoms and bonds
      frame_script += "color " + opti::GetStepStatusName( STEP_STATUS) + " #1\n"; // colors the model
      frame_script += "ribbon #1\n"; // shows it as ribbon
      frame_script += "ribrepr supersmooth #1\n"; // supersmooth ribbon

      // label
      frame_script += "2dlabels change 1 color black xpos 0.05 ypos 0.95 size 16 text \"";
      frame_script += TableToLabel( LABEL) + "\"\n";

      // record frame
      frame_script += "movie record raytrace " + std::string( m_RayTracing ? "true" : "false") + " fformat png pattern " + m_FilePrefix + "* directory " + m_Path + "\n";
      frame_script += "wait " + util::Format()( s_WaitTimes[ STEP_STATUS]) + "\n";
      frame_script += "movie stop\n";
      frame_script += "close #1\n";

      // add frame script to MovieScript
      m_MovieScript += frame_script;
    }

    //! @brief add final frame to movie
    //! @param FILE_NAME the filename of the pdb the model was written to
    //! @param LABEL add label to frame
    void MoviePrinterChimera::AddFinalFrame
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
      frame_script += "open 1 " + FILE_NAME + "\n"; // opens file
      frame_script += "~show #1\n"; // hides atoms and bonds
      frame_script += "color " + opti::GetStepStatusName( opti::e_Improved) + " #1\n"; // colors the model
      frame_script += "ribbon #1\n"; // shows it as ribbon
      frame_script += "ribrepr supersmooth #1\n"; // supersmooth ribbon

      // label
      frame_script += "2dlabels change 1 color black xpos 0.05 ypos 0.95 size 16 text \"";
      frame_script += TableToLabel( LABEL) + "\"\n";

      // record frame
      frame_script += "movie record raytrace " + std::string( m_RayTracing ? "true" : "false") + " fformat png pattern " + m_FilePrefix + "* directory " + m_Path + "\n";
      frame_script += "wait " + util::Format()( GetFinalWaitTime()) + "\n";
      frame_script += "movie stop\n";

      // add frame script to MovieScript
      m_MovieScript += frame_script;
    }

    //! @brief set the watermark
    //! @param WATERMARK the watermark
    void MoviePrinterChimera::SetWaterMark( const std::string &WATERMARK)
    {
      // create label
      m_FirstScript += "2dlabels create 2\n";

      // label
      m_FirstScript += "2dlabels change 2 color grey xpos 0.05 ypos 0.05 size 16 text \"" + WATERMARK + "\"\n";
    }

    //! @brief create ffmpeg command line
    //! @return commandline to be called with ffmpeg to generate the movie from the generated pngs
    std::string MoviePrinterChimera::FFMpegCommandLine() const
    {
      const std::string command_line
      (
        "ffmpeg" // command
        " -r 25" // frame rate
        " -i " + m_Path + m_FilePrefix + "%05d.png" // pattern for files
        " -y" // overwriting existing file
        " -f mp4" // format
        " -s 720x480"    // size
        " -qscale 1"       // quality
        " " + m_Path + m_FilePrefix + "movie.mp4" // output name
      );

      // end
      return command_line;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the movie script to a file
    //! @param OSTREAM
    std::ostream &MoviePrinterChimera::WriteScript
    (
      std::ostream &OSTREAM
    ) const
    {
      // write first portion and all frames
      OSTREAM << m_FirstScript;

      // write all frames
      OSTREAM << m_MovieScript;

      // write encoding of movie
      OSTREAM << "movie encode mformat mpeg output " + m_Path + m_FilePrefix + "movie.mpeg preset dvd resetMode keep\n";

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoviePrinterChimera::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoviePrinterChimera::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert table to label
    //! @param TABLE
    //! @return string
    std::string MoviePrinterChimera::TableToLabel( const storage::Table< double> &TABLE) const
    {
      // write table to string
      std::stringstream tmp;
      TABLE.WriteSubTable( tmp, m_LabelColNames, m_LabelRowNames);

      // copy label to modify
      std::string label( tmp.str());
      static const std::string s_line_break( "\\n");
      // iterate over label and replace line breaks with in text line breaks
      for( std::string::size_type itr( 0); itr != label.length(); ++itr)
      {
        // replace line breaks with actual line breaks
        if( label[ itr] == '\n')
        {
          label.replace( itr, 1, s_line_break);
          itr += s_line_break.length() - 1;
        }
      }

      // end
      return label;
    }

  } // namespace mc
} // namespace bcl
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
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "mc/bcl_mc_movie_printers.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! @brief get the default bcl watermark
    const std::string &MoviePrinterInterface::GetDefaultWatermark()
    {
      return GetCopyright();
    }

    //! @brief return command line flag for printing the script necessary to create a movie
    //! @return command line flag for printing the script necessary to create a movie
    util::ShPtr< command::FlagInterface> &MoviePrinterInterface::GetFlagMoviePrinter()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "movie",
          "\twrite a movie script and change preferences for movie generation"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( GetParameterMoviePrefix());
        flag->PushBack( GetParameterMoviePrinterType());
        flag->PushBack( GetParameterMovieWidth());
        flag->PushBack( GetParameterMovieHeight());
        flag->PushBack( GetParameterRayTrace());
      }

      return s_flag;
    }

    //! @brief return command line parameter for specifying the prefix for where movie files are generated
    //! @return command line parameter for specifying the prefix for where movie files are generated
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMoviePrefix()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_prefix",
          "prefix for where movie files are generated",
          ""
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the type of movie printer that should be used
    //! @return command line parameter for specifying the type of movie printer that should be used
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMoviePrinterType()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movieprinter",
          "choice of movie printer for different programs",
          command::ParameterCheckEnumerate< MoviePrinters>(),
          GetMoviePrinters().e_Pymol
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the width resolution of the movie
    //! @return command line parameter for specifying the width resolution of the movie
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMovieWidth()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_width",
          "width of movie",
          command::ParameterCheckRanged< size_t>( 480, 1920),
          "720"
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the height resolution of the movie
    //! @return command line parameter for specifying the height resolution of the movie
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMovieHeight()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_height",
          "height of movie",
          command::ParameterCheckRanged< size_t>( 320, 1200),
          "480"
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying whether or not all movie frames should be ray traced
    //! @return command line parameter for specifying whether or not all movie frames should be ray traced
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterRayTrace()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "ray",
          "raytrace all frames - 0=no, 1=yes",
          command::ParameterCheckRanged< size_t>( 0, 1),
          "0"
        )
      );

      return s_parameter;
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_movie_printers.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_movie_printer_chimera.h"
#include "mc/bcl_mc_movie_printer_pymol.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Default constructor
    MoviePrinters::MoviePrinters() :
      e_Chimera( AddEnum( "Chimera", util::ShPtr< MoviePrinterInterface>( new MoviePrinterChimera()))),
      e_Pymol( AddEnum( "Pymol", util::ShPtr< MoviePrinterInterface>( new MoviePrinterPymol())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoviePrinters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get enumerated list of MoviePrinters
    const MoviePrinters &GetMoviePrinters()
    {
      return MoviePrinters::GetEnums();
    }

  } // namespace mc

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< mc::MoviePrinterInterface>, mc::MoviePrinters>;

  } // namespace util
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_add.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopAdd::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopAdd)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopAdd::MutateLoopAdd( const std::string &LOOP_LIBRARY_FILENAME) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopAdd
    MutateLoopAdd *MutateLoopAdd::Clone() const
    {
      return new MutateLoopAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopAdd::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopAdd::GetAlias() const
    {
      static const std::string s_name( "MutateLoopAdd");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopAdd::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Constructs a missing loop region in a protein model using a template library.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopAdd::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopAdd::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the missing loops for construction
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // find a suitable loop template in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( loop));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        BCL_MessageVrb( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_add_resize.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopAddResize::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopAddResize())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    //! @param MIN_SIZES minimum lengths of helices and strands after resize
    MutateLoopAddResize::MutateLoopAddResize
    (
      const std::string &LOOP_LIBRARY_FILENAME,
      const storage::VectorND< 2, size_t> MIN_SIZES
    ) :
      m_LoopLocator(),
      m_MinSizes( MIN_SIZES),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopAddResize
    MutateLoopAddResize *MutateLoopAddResize::Clone() const
    {
      return new MutateLoopAddResize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopAddResize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopAddResize::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief returns the default minimum lengths of helices and strands after resize
    //! @return the default minimum lengths of helices and strands after resize
    storage::VectorND< 2, size_t> MutateLoopAddResize::GetDefaultMinSizes()
    {
      static storage::VectorND< 2, size_t> s_min_sizes( 5, 3);
      return s_min_sizes;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopAddResize::GetAlias() const
    {
      static const std::string s_name( "MutateLoopAddResize");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopAddResize::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Constructs a missing loop region in a protein model using a template library while randomly resizing the anchor SSEs."
      );
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopAddResize::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @detail this mutate randomly selects a missing loop in a protein model. the anchor SSEs of the selected loop are
    //! randomly resized and the resulting loop is constructed using a suitable conformation from a template library.
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopAddResize::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the missing loops for construction
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // get the anchor SSEs of this loop
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > anchors( GetAnchors( MODEL, loop));

      // randomly resize the anchor SSEs
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > resized_sses( ResizeSSEs( anchors));
      const util::ShPtr< assemble::SSE> &sp_n_anchor( resized_sses.First());
      const util::ShPtr< assemble::SSE> &sp_c_anchor( resized_sses.Second());

      // compute the new loop for the resized anchors
      const util::ShPtr< fold::LoopParameters> sp_new_loop
      (
        fold::LoopParameters::Create( *sp_n_anchor->GetLastMember(), *sp_c_anchor->GetFirstMember())
      );

      // replace the anchor SSEs with the resized ones
      result_model->ReplaceResize( sp_n_anchor);
      result_model->ReplaceResize( sp_c_anchor);

      // find a suitable loop template in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( loop));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        // return start model if no suitable templates were found
        BCL_MessageVrb( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the anchor SSEs of a given loop
    //! @param MODEL protein model in which to find the anchor SSEs
    //! @param LOOP loop for which to return the anchor SSEs
    //! @return anchor SSEs of the given loop
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateLoopAddResize::GetAnchors
    (
      const assemble::ProteinModel &MODEL,
      const fold::LoopParameters &LOOP
    )
    {
      // get sequence and chain IDs of the anchor residues
      const storage::Vector< int> &anchor_ids( LOOP.GetAnchors());
      const char &chain_id( LOOP.GetChainID());

      // find the corresponding SSEs
      const assemble::Chain &chain( *MODEL.GetChain( chain_id));
      const util::SiPtrVector< const biol::AABase> residues( chain.GetSequence()->GetMembers());
      const biol::AABase &n_res( *residues( anchor_ids( 0) - 1));
      const biol::AABase &c_res( *residues( anchor_ids( 1) - 1));
      const util::ShPtr< assemble::SSE> sp_n_anchor( new assemble::SSE( *MODEL.GetSSE( n_res)));
      const util::ShPtr< assemble::SSE> sp_c_anchor( new assemble::SSE( *MODEL.GetSSE( c_res)));
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > anchors( sp_n_anchor, sp_c_anchor);

      return anchors;
    }

    //! @brief randomly resizes the given SSEs
    //! @param SSES the SSEs that shall be resized
    //! @return resized SSEs
    storage::VectorND< 2, util::ShPtr< assemble::SSE> > MutateLoopAddResize::ResizeSSEs
    (
      const storage::VectorND< 2, util::ShPtr< assemble::SSE> > &SSES
    )
    {
      // determine the sequence length of the loop and the length of the anchors
      const assemble::SSE &n_anchor( *SSES.First());
      const assemble::SSE &c_anchor( *SSES.Second());

      // determine which anchors to resize: 0 for both, 1 for n-terminal, 2 for c-terminal
      const size_t resize_targets( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)));

      // determine the resize ranges
      const size_t resize_n
      (
        n_anchor.GetSize() > 6 ? random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)) : 0
      );
      const size_t resize_c
      (
        c_anchor.GetSize() > 6 ? random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, 2)) : 0
      );

      // resize the SSEs
      storage::VectorND< 2, util::ShPtr< assemble::SSE> > resized_sses;
      // if( resize_targets == 0 || resize_targets == 2) // resize n-terminal SSE
      if( false)
      {
        resized_sses.First() = ShrinkSSE( n_anchor, resize_n, biol::AASequenceFlexibility::e_CTerminal);
      }
      else
      {
        resized_sses.First() = SSES.First();
      }
      if( resize_targets == 1 || resize_targets == 2) // resize c-terminal SSE
      {
        resized_sses.Second() = ShrinkSSE( c_anchor, resize_c, biol::AASequenceFlexibility::e_NTerminal);
      }
      else
      {
        resized_sses.Second() = SSES.Second();
      }

      return resized_sses;
    }

    //! @brief shrinks the given SSE by the given length
    //! @param SSE SSE to be shrunk
    //! @param LENGTH length by which the SSE will get shrunk
    //! @param DIRECTION terminus at which to shrink the SSE
    //! @return the shrunk SSE
    util::ShPtr< assemble::SSE> MutateLoopAddResize::ShrinkSSE
    (
      const assemble::SSE &SSE, size_t LENGTH, const biol::AASequenceFlexibility::SequenceDirection &DIRECTION
    )
    {
      // check if the SSE is longer than the size reduction
      const size_t new_length( SSE.GetSize() - LENGTH);
      BCL_Assert( new_length > 0, "Operation would result in SSE of negative length.");

      // create a new shrunk SSE
      util::ShPtr< assemble::SSE> sp_shrunk_sse;
      if( DIRECTION == biol::AASequenceFlexibility::e_NTerminal)
      {
        sp_shrunk_sse = util::ShPtr< assemble::SSE>
          (
            new assemble::SSE( SSE.SubSequence( LENGTH, SSE.GetLength()), SSE.GetType())
          );
      }
      else if( DIRECTION == biol::AASequenceFlexibility::e_CTerminal)
      {
        sp_shrunk_sse = util::ShPtr< assemble::SSE>
          (
            new assemble::SSE( SSE.SubSequence( 0, new_length), SSE.GetType())
          );
      }
      else
      {
        BCL_MessageCrt( "Direction " + util::Format()( DIRECTION) + " is unsupported. Returning original SSE.");
        sp_shrunk_sse = util::CloneToShPtr( SSE);
      }

      return sp_shrunk_sse;
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_fragment_add.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopFragmentAdd::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopFragmentAdd)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopFragmentAdd::MutateLoopFragmentAdd( const std::string &LOOP_LIBRARY_FILENAME) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopFragmentAdd
    MutateLoopFragmentAdd *MutateLoopFragmentAdd::Clone() const
    {
      return new MutateLoopFragmentAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopFragmentAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopFragmentAdd::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopFragmentAdd::GetAlias() const
    {
      static const std::string s_name( "MutateLoopFragmentAdd");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopFragmentAdd::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Constructs a fragment of a missing loop in a protein model.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopFragmentAdd::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopFragmentAdd::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      typedef fold::LocatorMissingCoordinates::Span Span;
      const storage::Vector< Span> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the missing loops for construction
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const Span &loop( loops( loop_index));
      const size_t seq_distance( loop.Second() - loop.First() + 1);
      if( seq_distance < 3)
      {
        BCL_MessageVrb( "No missing loops with at least length 3 in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // find a suitable loop template in the library
      const math::Range< size_t> fragment_range( 1, seq_distance - 1);
      const size_t fragment_length( random::GetGlobalRandom().SizeT( fragment_range));
      storage::Vector< int> new_anchors;
      new_anchors.PushBack( loop.First() - 1);
      new_anchors.PushBack( loop.First() + fragment_length);
      const fold::LoopParameters loop_fragment( new_anchors, loop.Third());
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( fragment_length));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        BCL_MessageTop( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop_fragment, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_fragment_replace.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopFragmentReplace::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopFragmentReplace)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopFragmentReplace::MutateLoopFragmentReplace( const std::string &LOOP_LIBRARY_FILENAME) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopFragmentReplace
    MutateLoopFragmentReplace *MutateLoopFragmentReplace::Clone() const
    {
      return new MutateLoopFragmentReplace( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopFragmentReplace::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopFragmentReplace::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopFragmentReplace::GetAlias() const
    {
      static const std::string s_name( "MutateLoopFragmentReplace");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopFragmentReplace::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Replaces a fragment of a loop in a protein model.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopFragmentReplace::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopFragmentReplace::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));
      util::ShPtrVector< fold::LoopParameters> loops_filtered;
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        const fold::LoopParameters &loop_params( **loop_it);
        const char chain_id( loop_params.GetChainID());
        const storage::Pair< int, int> res_ids( loop_params.GetAnchors()( 0) + 1, loop_params.GetAnchors()( 1) - 1);
        const size_t length( res_ids.Second() - res_ids.First());
        const biol::AASequence loop( MODEL.GetChain( chain_id)->GetSequence()->SubSequence( res_ids.First() - 1, length));
        if( IsPartiallyDefined( loop))
        {
          loops_filtered.PushBack( *loop_it);
        }
      }

      // randomly select one fragment for reconstruction
      if( loops_filtered.IsEmpty())
      {
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      const math::Range< size_t> range( 0, loops_filtered.GetSize() - 1);
      const size_t loop_index( random::GetGlobalRandom().SizeT( range));
      const fold::LoopParameters &selected_loop( *loops_filtered( loop_index));
      const char chain_id( selected_loop.GetChainID());
      const storage::Pair< int, int> res_ids( selected_loop.GetAnchors()( 0), selected_loop.GetAnchors()( 1));
      const biol::AASequence loop_sequence
      (
        MODEL.GetChain( chain_id)->GetSequence()->SubSequence( res_ids.First(), res_ids.Second() - res_ids.First())
      );
      const util::ShPtr< fold::LoopParameters> selected_fragment( SelectFragment( loop_sequence));
      if( !selected_fragment.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      const size_t frag_length( selected_fragment->GetSequenceDistance() + 2);
      const math::Range< size_t> frag_range( 1, frag_length);
      const size_t selected_frag_length( random::GetGlobalRandom().SizeT( frag_range));
      storage::Vector< int> rebuild_ids;
      rebuild_ids.PushBack( selected_fragment->GetAnchors()( 1) - selected_frag_length);
      rebuild_ids.PushBack( selected_fragment->GetAnchors()( 1) + 1);
      const fold::LoopParameters fragment_rebuild( rebuild_ids, selected_fragment->GetChainID());

      // find a suitable template for the selected fragment in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( selected_frag_length));
      if( templates.IsEmpty())
      {
        BCL_MessageTop( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }
      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, fragment_rebuild, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Returns whether the provided sequence is partially defined
    //! @param SEQUENCE sequence to be evaluated
    //! @return true, if the provided sequence is partially defined
    bool MutateLoopFragmentReplace::IsPartiallyDefined( const biol::AASequence &SEQUENCE) const
    {
      for( auto res_it( SEQUENCE.Begin()), res_it_end( SEQUENCE.End()); res_it != res_it_end; ++res_it)
      {
        if( ( **res_it).HasDefinedCoordinates())
        {
          for( auto res_int_it( res_it + 1); res_int_it != res_it_end; ++res_int_it)
          {
            if( ( **res_int_it).HasDefinedCoordinates())
            {
              return true;
            }
          }
        }
      }
      return false;
    }

    //! @brief Randomly select a fragment for reconstruction
    //! @param SEQUENCE sequence to select from
    //! @return randomly selected fragment
    util::ShPtr< fold::LoopParameters> MutateLoopFragmentReplace::SelectFragment( const biol::AASequence &SEQUENCE) const
    {
      const size_t n_id( ( **SEQUENCE.Begin()).GetSeqID());
      for( auto res_it( SEQUENCE.Begin() + 1), res_it_end( SEQUENCE.End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &current_res( **res_it);
        if( !current_res.HasDefinedCoordinates())
        {
          storage::Vector< int> anchors;
          anchors.PushBack( n_id);
          anchors.PushBack( current_res.GetSeqID() - 1);
          const util::ShPtr< fold::LoopParameters> sp_frag( new fold::LoopParameters( anchors, current_res.GetChainID()));
          return sp_frag;
        }
      }
      return util::ShPtr< fold::LoopParameters>();
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_remove.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_complete.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopRemove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopRemove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateLoopRemove::MutateLoopRemove() :
      m_LoopLocator( true, false) // locate defined loops and ignore terminal loops
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopRemove
    MutateLoopRemove *MutateLoopRemove::Clone() const
    {
      return new MutateLoopRemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopRemove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopRemove::GetAlias() const
    {
      static const std::string s_name( "MutateLoopRemove");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopRemove::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Removes a loop from a protein model.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopRemove::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLocator = fold::LocatorLoop();
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopRemove::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // clone the original model
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find defined loop regions in the model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are now defined loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No defined loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the loop regions
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // find the loop in the protein model
      const util::ShPtr< assemble::Chain> &sp_chain( result_model->GetChain( loop.GetChainID()));
      const util::SiPtrVector< const biol::AABase> res_tmp( sp_chain->GetAminoAcids());
      const biol::AABase &first_loop_res( *res_tmp( loop.GetAnchors()( 0)));
      const assemble::SSE &selected_loop( *result_model->GetSSE( first_loop_res));
      const util::SiPtrVector< const biol::AABase> residues( selected_loop.GetMembers());

      // create a copy of the selected loop with undefined coordinates
      util::ShPtrVector< biol::AABase> new_seq;
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        const util::ShPtr< biol::AABase> sp_new_res( new biol::AAComplete( ( **res_it).GetData()));
        new_seq.PushBack( sp_new_res);
      }
      const biol::AASequence new_sequence( new_seq, selected_loop.GetFirstAA()->GetChainID());
      const util::ShPtr< assemble::SSE> sp_new_loop( new assemble::SSE( new_sequence, selected_loop.GetType()));

      // replace the selected loop in the protein model with the new undefined one
      result_model->Replace( sp_new_loop);
      util::ShPtr< assemble::SSE> sp_loop( new assemble::SSE( selected_loop, biol::GetSSTypes().COIL));

      // return the result
      const math::MutateResult< assemble::ProteinModel> mutate_result( result_model, *this);
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutate_loop_replace.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protein_geometry.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> MutateLoopReplace::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateLoopReplace())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param LOOP_LIBRARY_FILENAME path of the loop template library
    MutateLoopReplace::MutateLoopReplace
    (
      const std::string &LOOP_LIBRARY_FILENAME
    ) :
      m_LoopLocator(),
      m_LoopLibraryFilename( LOOP_LIBRARY_FILENAME)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief copy constructor
    //! @return pointer to a new MutateLoopReplace
    MutateLoopReplace *MutateLoopReplace::Clone() const
    {
      return new MutateLoopReplace( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &MutateLoopReplace::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default path of the loop template library
    //! @return the default path of the loop template library
    const std::string &MutateLoopReplace::GetDefaultLibraryPath()
    {
      static std::string s_default_path( "histogram/loop_library");
      return s_default_path;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopReplace::GetAlias() const
    {
      static const std::string s_name( "MutateLoopReplace");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopReplace::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Replaces a loop in a protein model with a conformation from the template library.");
      serializer.AddInitializer
      (
        "loop library",
        "path to the loop template library",
        io::Serialization::GetAgentInputFilename( &m_LoopLibraryFilename)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateLoopReplace::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_LoopLibrary = fold::LoopLibrary::CreateLoopLibrary( m_LoopLibraryFilename);
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief applies a mutation to the given protein model
    //! @param MODEL protein model to which to apply the mutation
    //! @return result mutating the given protein model
    math::MutateResult< assemble::ProteinModel> MutateLoopReplace::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // model resulting from applying the mutate
      util::ShPtr< assemble::ProteinModel> result_model( MODEL.Clone());

      // find missing loop regions in the protein model
      const util::ShPtrVector< fold::LoopParameters> loops( m_LoopLocator.Locate( MODEL));

      // if there are no missing loops in the protein model return the original model
      if( loops.IsEmpty())
      {
        BCL_MessageVrb( "No missing loops in the protein model. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // randomly select one of the existing loops for replacement
      const size_t loop_index( random::GetGlobalRandom().SizeT( math::Range< size_t>( 0, loops.GetSize() - 1)));
      const fold::LoopParameters &loop( *loops( loop_index));

      // find a suitable loop template in the library
      const util::ShPtrVector< fold::LoopParameters> templates( m_LoopLibrary->FindTemplates( loop));

      // loop region can't be constructed if there is no suitable template
      if( templates.IsEmpty())
      {
        BCL_MessageVrb( "No suitable loop template found. Returning original model.");
        return math::MutateResult< assemble::ProteinModel>( result_model, *this);
      }

      // select a random template from the found templates
      const util::ShPtrVector< fold::LoopParameters>::const_iterator selected_template_it
      (
        random::GetGlobalRandom().Iterator( templates.Begin(), templates.End(), templates.GetSize())
      );
      const fold::LoopParameters &selected_template( **selected_template_it);

      // fit the loop to the selected template and add it to the model
      result_model = fold::ProteinGeometry::FitToTemplate( MODEL, loop, selected_template);

      // return the result
      return math::MutateResult< assemble::ProteinModel>( result_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_mutates.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_mutate_loop_add.h"
#include "mc/bcl_mc_mutate_loop_add_resize.h"
#include "mc/bcl_mc_mutate_loop_remove.h"
#include "mc/bcl_mc_mutate_loop_replace.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> Mutates::s_Instance
    (
      GetObjectInstances().AddInstance( new Mutates())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mutates::Mutates()
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new Mutates
    Mutates *Mutates::Clone() const
    {
      return new Mutates( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @returns single instance of this class
    //! @return single instance of this class
    Mutates &Mutates::GetInstance()
    {
      static Mutates s_single_instance;
      return s_single_instance;
    }

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Mutates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the mutates and adds them to the enumerator
    void Mutates::Initialize()
    {
      // initialize the default mutates
      fold::DefaultMutates::GetInstance().InitializeMutates();

      // initialize mutates for loop hashing
      InitializeLoopHashMutates();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &Mutates::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &Mutates::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize mutates required for loop hashing
    void Mutates::InitializeLoopHashMutates()
    {
      // get the enumerated mutates
      fold::Mutates &mutates( fold::GetMutates());

      // add a mutate for adding loops
      e_MutateLoopHashAdd = mutates.AddMutate( MutateLoopAdd());

      // add a mutate for replacing loops
      e_MutateLoopHashReplace = mutates.AddMutate( MutateLoopReplace());

      // add a mutate for resizing loops
      e_MutateLoopHashResize = mutates.AddMutate( MutateLoopAddResize());

      // add a mutate for removing loops
      e_MutateLoopHashRemove = mutates.AddMutate( MutateLoopRemove());
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_optimization_ccd.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_handler_locator_loop_domain.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> OptimizationCCD::s_Instance
    (
      util::Enumerated< opti::OptimizationInterface< assemble::ProteinModel> >::AddInstance( new OptimizationCCD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OptimizationCCD::OptimizationCCD()
    {
    }

    //! @brief construct from members
    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param MUTATES mutates to sample models
    //! @param CRITERION termination criterion
    //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
    OptimizationCCD::OptimizationCCD
    (
      const score::ProteinModelScoreSum &SCORE_FUNCTION,
      const math::MutateInterface< assemble::ProteinModel> &MUTATES,
      const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
      const Metropolis< double> &METROPOLIS
    ) :
      OptimizationMCM( SCORE_FUNCTION, MUTATES, CRITERION, METROPOLIS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new OptimizationCCD
    OptimizationCCD *OptimizationCCD::Clone() const
    {
      return new OptimizationCCD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &OptimizationCCD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &OptimizationCCD::GetAlias() const
    {
      static const std::string s_alias( "CCDOptimizer");
      return s_alias;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs pre-processing of the argument before the optimization is performed
    //! @detail adds initial coordinates for missing loop regions and adds a loop domain locator to the protein model data
    //! @param ARGUMENT data to be pre-processed
    //! @return pre-processed data
    void OptimizationCCD::PreProcess( assemble::ProteinModel &ARGUMENT) const
    {
      const fold::ProtocolLoopClose &protocol( fold::ProtocolLoopClose::GetInstance());

      // complete the model with coordinates in missing places
      protocol.AddLoopCoordinates( ARGUMENT);

      // split the coils, where they are not peptide bonded
      protocol.SplitCoilsAtNonPetideBond( ARGUMENT);

      // add backbone hydrogens
      ARGUMENT = *protocol.AddNitrogenHydrogens( ARGUMENT);

      // add a loop locator to the protein model data
      const util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_locator
      (
        fold::HandlerLocatorLoopDomain().CreateBidirectionalLocatorsForInteriorCoil( ARGUMENT)
      );
      util::ShPtr< assemble::ProteinModelData> sp_model_data( ARGUMENT.GetProteinModelData());
      if( !sp_model_data->GetData( assemble::ProteinModelData::e_LoopDomainLocators).IsDefined())
      {
        sp_model_data->Insert( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
      else
      {
        sp_model_data->Replace( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace mc
} // namespace bcl
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

#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_mutate_membrane_chain_move.h"
#include "fold/bcl_fold_mutate_protein_model_chain_move.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_optimization_docking.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_quality_docking.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "quality/bcl_quality_rmsd.h"
namespace bcl
{
    namespace mc
    {
  //////////
  // data //
  //////////

        //! single instance of this class
        const util::SiPtr< const util::ObjectInterface> OptimizationDocking::s_instance
        (
          GetObjectInstances().AddInstance( new OptimizationDocking())
        );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

        //! @brief default constructor
        OptimizationDocking::OptimizationDocking()
        {
        }

        //! @brief constructor that takes in initializers
        //! @param SCORE_FUNCTION scoring function for ranking docked models
        //! @param MUTATES Monte Carlo moves for sampling conformational space
        //! @param CRITERION termination criterion
        //! @param METROPOLIS Metropolis parameters
        //! @param MEASURES what measures to use to evaluate quality of models
        OptimizationDocking::OptimizationDocking
        (
          const score::ProteinModelScoreSum &SCORE_FUNCTION,
          const math::MutateInterface< assemble::ProteinModel> &MUTATES,
          const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
          const Metropolis< double> &METROPOLIS,
          const storage::Set< quality::Measure> &MEASURES
        ) :
          m_ScoreFunction( SCORE_FUNCTION),
          m_Mutates( MUTATES),
          m_Criterion( CRITERION),
          m_Metropolis( METROPOLIS),
          m_QualityMeasures( MEASURES)
        {
        }

        //! @brief clone function
        //! @return a pointer to a new OptimizationMCMDocking object
        OptimizationDocking *OptimizationDocking::Clone() const
        {
            return new OptimizationDocking( *this);
        }

  /////////////////
  // data access //
  /////////////////

        //! @brief gets the name of this class
        //! @return the name of this class
        const std::string &OptimizationDocking::GetClassIdentifier() const
        {
            return GetStaticClassName( *this);
        }

        //! @brief gets the name of the object when used in a dynamic context
        //! @return the name of the object when used in a dynamic context
        const std::string &OptimizationDocking::GetAlias() const
        {
            static const std::string s_alias( "DockingOptimizer");
            return s_alias;
        }

        //! @brief get the scoring function
        //! @return the scoring function
        const score::ProteinModelScoreSum &OptimizationDocking::GetScoreFunction() const
        {
          return m_ScoreFunction;
        }

        //! @brief get quality measures
        //! @return quality measures
        const storage::Set< quality::Measure> &OptimizationDocking::GetQualityMeasures() const
        {
          return m_QualityMeasures;
        }

        //! @brief gets parameters for data members that are set up from data labels
        //! @return parameters for data members that are set up from data labels
        io::Serializer OptimizationDocking::GetSerializer() const
        {
            // gets the serializer from the base class
            io::Serializer serializer;

            // sets class description
            serializer.SetClassDescription( "Monte Carlo sampling for protein-protein docking");

            // add initializers
            serializer.AddInitializer
            (
              "score function",
              "scoring function to be used along with Monte Carlo sampling",
              io::Serialization::GetAgent( &m_ScoreFunction)
            );
            serializer.AddInitializer
            (
              "mutates",
              "mutates to sample the protein models",
              io::Serialization::GetAgent( &m_Mutates)
            );
            serializer.AddInitializer
            (
              "termination criterion",
              "criterion when the optimization will be terminated",
              io::Serialization::GetAgent( &m_Criterion)
            );
            serializer.AddInitializer
            (
              "metropolis",
              "metropolis criterion to decide which mutates are accepted",
              io::Serialization::GetAgent( &m_Metropolis)
            );
            serializer.AddInitializer
            (
              "quality measures",
              "measures for evaluating how close the model is to the native structure",
              io::Serialization::GetAgent( &m_QualityMeasures),
              "(RMSD)"
            );

            return serializer;
        }

  ////////////////
  // operations //
  ////////////////

        //! @brief optimizes the provided protein model
        //! @param MODEL protein model to be optimized
        void OptimizationDocking::Optimize( assemble::ProteinModel &MODEL) const
        {
          // create the approximator
          Approximator< assemble::ProteinModel, double> approximator
          (
            m_ScoreFunction, *m_Mutates, m_Metropolis, *m_Criterion, MODEL
          );

          // optimize the given protein model
          approximator.Approximate();

          // get the result of the optimization
          MODEL = *util::ShPtr< assemble::ProteinModel>( approximator.GetTracker().GetBest()->First().HardCopy());
        }
  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_optimization_mcm.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_approximator.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> OptimizationMCM::s_Instance
    (
      util::Enumerated< opti::OptimizationInterface< assemble::ProteinModel> >::AddInstance( new OptimizationMCM())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OptimizationMCM::OptimizationMCM() :
      m_ScoreFunction(),
      m_Mutates(),
      m_Criterion(),
      m_Metropolis()
    {
    }

    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param MUTATES mutates to sample models
    //! @param CRITERION termination criterion
    //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
    OptimizationMCM::OptimizationMCM
    (
      const score::ProteinModelScoreSum &SCORE_FUNCTION,
      const math::MutateInterface< assemble::ProteinModel> &MUTATES,
      const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
      const Metropolis< double> &METROPOLIS
    ) :
      m_ScoreFunction( SCORE_FUNCTION),
      m_Mutates( MUTATES),
      m_Criterion( CRITERION),
      m_Metropolis( METROPOLIS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new OptimizationMCM
    OptimizationMCM *OptimizationMCM::Clone() const
    {
      return new OptimizationMCM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &OptimizationMCM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &OptimizationMCM::GetAlias() const
    {
      static const std::string s_alias( "MCMOptimizer");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer OptimizationMCM::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
      // serializer.Merge( opti::OptimizationInterface< assemble::ProteinModel>::GetSerializer());
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the sampled protein models",
        io::Serialization::GetAgent( &m_ScoreFunction)
      );
      serializer.AddInitializer
      (
        "mutates",
        "mutates to sample the protein models",
        io::Serialization::GetAgent( &m_Mutates)
      );
      serializer.AddInitializer
      (
        "termination criterion",
        "criterion when the optimization will be terminated",
        io::Serialization::GetAgent( &m_Criterion)
      );
      serializer.AddInitializer
      (
        "metropolis",
        "metropolis criterion to decide which mutates are accepted",
        io::Serialization::GetAgent( &m_Metropolis)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief optimizes the provided protein model
    //! @param MODEL protein model to be optimized
    void OptimizationMCM::Optimize( assemble::ProteinModel &MODEL) const
    {
      // create the approximator
      Approximator< assemble::ProteinModel, double> approximator
      (
        m_ScoreFunction, *m_Mutates, m_Metropolis, *m_Criterion, MODEL
      );

      // optimize the given protein model
      approximator.Approximate();

      // get the result of the optimization
      MODEL = *approximator.GetTracker().GetBest()->First().HardCopy();
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the default temperature control for Metropolis
    //! @return default temperature control for Metropolis
    util::ShPtr< TemperatureInterface> OptimizationMCM::GetDefaultTemperature()
    {
      return util::ShPtr< TemperatureInterface>( new TemperatureDefault);
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_stage.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutates.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_printer_combined.h"
#include "mc/bcl_mc_printer_with_criterion.h"
#include "opti/bcl_opti_criterion_all.h"
#include "opti/bcl_opti_criterion_phase.h"
#include "opti/bcl_opti_criterion_result_changed.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Stage::s_Instance
    (
      GetObjectInstances().AddInstance( new Stage())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Stage::Stage() :
      m_Name( ""),
      m_Number( util::GetUndefined< size_t>()),
      m_FoldProtocols(),
      m_ScoreProtocols(),
      m_MutateProtocols(),
      m_MaxNumberIterations( util::GetUndefined< size_t>()),
      m_MaxNumberUnimprovedIterations( util::GetUndefined< size_t>()),
      m_ScoreFunction(),
      m_ScoreWeightSet(),
      m_MutateTree(),
      m_Mutate(),
      m_Temperature(),
      m_NumberOfScoresDropped( 0),
      m_ModifyStartModel( true),
      m_PrintStartModel( false),
      m_PrintIterationModels( false),
      m_PrintEndModel( false),
      m_PrintTrackerHistory( false),
      m_PoolPostfix( ""),
      m_RoundNumber(util::GetUndefined< size_t>()),
      m_Prefix( ""),
      m_Path(),
      m_QualityMeasures()
    {
    }

    //! @brief constructor from members
    //! @param FOLD_PROTOCOLS Fold protocols
    //! @param SCORE_PROTOCOLS Score protocols
    //! @param MUTATE_PROTOCOLS Mutate protocols
    Stage::Stage
    (
      const storage::List< fold::Protocol> &FOLD_PROTOCOLS,
      const storage::List< fold::Protocol> &SCORE_PROTOCOLS,
      const storage::List< fold::Protocol> &MUTATE_PROTOCOLS
    ) :
      m_Name( ""),
      m_Number( util::GetUndefined< size_t>()),
      m_FoldProtocols( FOLD_PROTOCOLS),
      m_ScoreProtocols( SCORE_PROTOCOLS),
      m_MutateProtocols( MUTATE_PROTOCOLS),
      m_MaxNumberIterations( util::GetUndefined< size_t>()),
      m_MaxNumberUnimprovedIterations( util::GetUndefined< size_t>()),
      m_ScoreFunction(),
      m_ScoreWeightSet(),
      m_MutateTree(),
      m_Mutate(),
      m_Temperature(),
      m_NumberOfScoresDropped( 0),
      m_ModifyStartModel( true),
      m_PrintStartModel( false),
      m_PrintIterationModels( false),
      m_PrintEndModel( false),
      m_PrintTrackerHistory( false),
      m_PoolPostfix( ""),
      m_RoundNumber(util::GetUndefined< size_t>()),
      m_Prefix( ""),
      m_Path(),
      m_QualityMeasures()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Stage
    Stage *Stage::Clone() const
    {
      return new Stage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Stage::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the list of fold protocols used
    //! @return the list of fold protocols used
    const storage::List< fold::Protocol> &Stage::GetFoldProtocols() const
    {
      return m_FoldProtocols;
    }

    //! @brief sets the list of fold protocols used
    //! @param FOLD_PROTOCOLS list of fold protocols
    void Stage::SetFoldProtocols( const storage::List< fold::Protocol> &FOLD_PROTOCOLS)
    {
      m_FoldProtocols = FOLD_PROTOCOLS;
    }

    //! @brief returns the list of score protocols used
    //! @return the list of score protocols used
    const storage::List< fold::Protocol> &Stage::GetScoreProtocols() const
    {
      return m_ScoreProtocols;
    }

    //! @brief sets the list of score protocols used
    //! @param SCORE_PROTOCOLS list of score protocols
    void Stage::SetScoreProtocols( const storage::List< fold::Protocol> &SCORE_PROTOCOLS)
    {
      m_ScoreProtocols = SCORE_PROTOCOLS;
    }

    //! @brief returns the list of mutate protocols used
    //! @return the list of mutate protocols used
    const storage::List< fold::Protocol> &Stage::GetMutateProtocols() const
    {
      return m_MutateProtocols;
    }

    //! @brief sets the list of mutate protocols used
    //! @param MUTATE_PROTOCOLS list of mutate protocols
    void Stage::SetMutateProtocols( const storage::List< fold::Protocol> &MUTATE_PROTOCOLS)
    {
      m_MutateProtocols = MUTATE_PROTOCOLS;
    }

    //! @brief returns score weightset
    //! @return scoring weightset
    const util::ShPtr< fold::ScoreWeightSet> &Stage::GetScoreWeightSet() const
    {
      return m_ScoreWeightSet;
    }

    //! @brief returns score weightset
    //! @return scoring weightset
    void Stage::SetScoreWeightSet( const util::ShPtr< fold::ScoreWeightSet> &SCORE_WEIGHT_SET)
    {
      m_ScoreWeightSet = SCORE_WEIGHT_SET;
    }

    //! @brief returns score function
    //! @return score function
    const util::ShPtr< score::ProteinModelScoreSum> &Stage::GetScoreFunction() const
    {
      return m_ScoreFunction;
    }

    //! @brief returns score function
    //! @return score function
    void Stage::SetScoreFunction( const util::ShPtr< score::ProteinModelScoreSum> &SCORE_FUNCTION)
    {
      m_ScoreFunction = SCORE_FUNCTION;
    }

    //! @brief get mutate tree
    //! @return mutate tree
    const util::ShPtr< fold::MutateTree> &Stage::GetMutateTree() const
    {
      return m_MutateTree;
    }

    //! @brief get mutate tree
    //! @return mutate tree
    void Stage::SetMutateTree( const util::ShPtr< fold::MutateTree> &MUTATE_TREE)
    {
      m_MutateTree = MUTATE_TREE;
    }

    //! @brief get mutate functions for protein models
    //! @return mutate functions for protein models
    const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &Stage::GetMutate() const
    {
      return m_Mutate;
    }

    //! @brief get mutate functions for protein models
    //! @return mutate functions for protein models
    void Stage::SetMutate( const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &MUTATE)
    {
      m_Mutate = MUTATE;
    }

    //! @brief return the temperature
    //! @return the temperature
    const util::ShPtr< TemperatureInterface> &Stage::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief sets the temperature
    //! @param SP_TEMPERATURE new temperature to be usedif( is_last_stage)
    void Stage::SetTemperature( const util::ShPtr< TemperatureInterface> &SP_TEMPERATURE)
    {
      m_Temperature = SP_TEMPERATURE;
    }

    //! @brief sets the boolean indicating if the start model should be modified or not
    //! @param MODIFY_START_MODEL indicates true if the start model should be modified - false if not
    void Stage::SetModifyStartModel( const bool MODIFY_START_MODEL)
    {
      m_ModifyStartModel = MODIFY_START_MODEL;
    }

    //! @brief sets the boolean indicating if the start model should be printed or not
    //! @param MODIFY_START_MODEL indicates true if the start model should be printed - false if not
    void Stage::SetPrintStartModel( const bool PRINT_START_MODEL)
    {
      m_PrintStartModel = PRINT_START_MODEL;
    }

    //! @brief sets the boolean indicating if every model should be printed or not
    //! @param PRINT_MODEL indicates true if the all models should be printed - false if not
    void Stage::SetPrintIterationModels( const bool PRINT_MODEL)
    {
      m_PrintIterationModels = PRINT_MODEL;
    }

    //! @brief sets the boolean indicating if the end model should be printed or not
    //! @param MODIFY_END_MODEL indicates true if the end model should be printed - false if not
    void Stage::SetPrintEndModel( const bool PRINT_END_MODEL)
    {
      m_PrintEndModel = PRINT_END_MODEL;
    }

    //! @brief gets the filename postfix of the pool for use with this stage
    //! @return string which is the filename postfix of the pool for use with this stage
    const std::string &Stage::GetPoolPostfix() const
    {
      return m_PoolPostfix;
    }

    //! @brief sets the filename postfix of the pool for use with this stage
    //! @param POOL_POSTFIX the filename postfix of the pool for use with this stage
    void Stage::SetPoolPostfix( const std::string &POOL_POSTFIX)
    {
      m_PoolPostfix = POOL_POSTFIX;
    }
    //! @brief return quality measures to be calculated
    //! @return quality measures to be calculated
    void Stage::SetQualityMeasures( const storage::Set< quality::Measure> QUALITY_MEASURES)
    {
       m_QualityMeasures = QUALITY_MEASURES;
    }
    //! @brief gives the prefix object
    //! @return the prefix that is prepended to output files
    void Stage::SetPrefix( const std::string &PREFIX)
    {
       m_Prefix = PREFIX;
    }

    //! @brief sets the boolean indicating if the tracker history should be printed
    //! @param TRACKER_HISTORY true if the tracker history should be printer
    void Stage::SetPrintTrackerHistory( const bool &TRACKER_HISTORY)
    {
      m_PrintTrackerHistory = TRACKER_HISTORY;
    }

    //! @brief sets the path for the output
    //! @param PATH path for the output
    void Stage::SetPath( const std::string &PATH)
    {
       m_Prefix = PATH;
    }

    //! @brief initiates the approximation and returns a shared pointer to the approximation result
    //! @param SP_PROTEIN_MODEL shared pointer to the starting model for the approximation
    //! @return shared pointer to the approximation result
    util::ShPtr< assemble::ProteinModel> Stage::Approximate( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const
    {
      // if a new sse pool is needed for this stage
      if( !m_PoolPostfix.empty())
      {
        InitializePool( SP_PROTEIN_MODEL);
      }

      // modify the start model
      if( m_ModifyStartModel)
      {
        BCL_MessageStd( "modifying start model for stage " + util::Format()( m_Number + 1));
        ModifyStartModel( *SP_PROTEIN_MODEL);
      }

    /////////////////////////////////////
    // setting up printers and storage //
    /////////////////////////////////////

      // create the combined printer
      PrinterCombined< assemble::ProteinModel, double> printer_combine;

      // create the criterion for the protein model printer
      storage::Set< opti::PhaseEnum> phases;
      if( m_PrintStartModel)
      {
        phases.Insert( opti::e_Start);
      }
      if( m_PrintIterationModels)
      {
        phases.Insert( opti::e_Iteration);
      }
      if( m_PrintEndModel)
      {
        phases.Insert( opti::e_End);
      }

      util::ShPtr< opti::CriterionInterface< assemble::ProteinModel, double> > sp_print_criterion
      (
        new opti::CriterionPhase< assemble::ProteinModel, double>( phases)
      );
      if( m_PrintIterationModels)
      {
        opti::CriterionAll< assemble::ProteinModel, double> comb;
        comb.InsertCriteria( *sp_print_criterion);
        comb.InsertCriteria( opti::CriterionResultChanged< assemble::ProteinModel, double>());
        sp_print_criterion = util::ShPtr< opti::CriterionInterface< assemble::ProteinModel, double> >( comb.Clone());
      }
      // create the protein model printer and add it to the combined printer
      util::ShPtr< PrintInterface< assemble::ProteinModel, double> > sp_printer
      (
        new assemble::PrinterProteinModel
        (
          fold::GetSetup().GetPrefix() + ( !m_IsLastStage ? m_Name : ""),
          fold::GetSetup().GetStorage(),
          fold::GetSetup().GetSuperimposeMeasure()
        )
      );
      util::ShPtr< PrintInterface< assemble::ProteinModel, double> > sp_model_printer
      (
        new PrinterWithCriterion< assemble::ProteinModel, double>( sp_printer, sp_print_criterion)
      );
      printer_combine.Insert( sp_model_printer);
      printer_combine.Initialize( m_RoundNumber, m_Number);
      printer_combine.SetPrefix( m_Prefix);

      // create the pdb factory
      util::ShPtr< pdb::Factory> factory( new pdb::Factory());
      ModifyFactory( factory);

      // modify the factory so that it adds the scoring table to the pdb when printed
      util::ShPtr< assemble::ProteinStorageFile> sp_storage( fold::GetSetup().GetStorage());
      sp_storage->SetFactory( factory);

    ////////////////////////////////////////////
    // preparing and performing approximation //
    ////////////////////////////////////////////

      // create the metropolis criterion
      const double metropolis_min_change( 0.0001);
      const util::ShPtr< Metropolis< double> > sp_metropolis
      (
        new Metropolis< double>( m_Temperature, true, metropolis_min_change)
      );

      // create the termination criteria for the approximation
      opti::CriterionCombine< assemble::ProteinModel, double> criterion_combine;
      ModifyCriterion( criterion_combine);

      storage::Vector< std::string> dropped_functions( m_NumberOfScoresDropped);
      storage::Vector< double> dropped_values( m_NumberOfScoresDropped);
      util::ShPtr< score::ProteinModelScoreSum> sum( m_ScoreFunction);
      static storage::Set< std::string> s_prev_dropped_functions;
      static bool s_score_dropped_continuity( false);
      if( m_NumberOfScoresDropped)
      {
        storage::Vector< std::string> functions;
        if( s_score_dropped_continuity)
        {
          storage::Vector< std::string> prev_dropped_functions;
          functions = m_ScoreFunction->GetFunctionSchemes();
          functions.PopBack();
          functions.Shuffle();
          size_t n_not_dropped( 0);
          for( auto itr_a( functions.Begin()), itr_place( functions.Begin()), itr_end( functions.End()); itr_a != itr_end; ++itr_a)
          {
            if( !s_prev_dropped_functions.Contains( *itr_a))
            {
              if( itr_a != itr_place)
              {
                *itr_place = *itr_a;
              }
              ++n_not_dropped;
              ++itr_place;
            }
            else
            {
              prev_dropped_functions.PushBack( *itr_a);
            }
          }
          functions.Resize( n_not_dropped);
          functions.InsertElements( 0, prev_dropped_functions);
          s_prev_dropped_functions.Reset();
          s_prev_dropped_functions.InsertElements( functions.Begin(), functions.Begin() + m_NumberOfScoresDropped);
        }
        else
        {
          functions = m_ScoreFunction->GetFunctionSchemes();
          functions.PopBack();
          functions.Shuffle();
        }
        for( size_t i( 0); i < m_NumberOfScoresDropped; ++i)
        {
          dropped_functions( i) = functions( i);
          dropped_values( i) = sum->GetTerm( functions( i)).First();
          sum->SetCoefficient( functions( i), 0.0);
        }
      }

      // create the approximator
      Approximator< assemble::ProteinModel, double> approximator
      (
        *m_ScoreFunction,
        *m_Mutate,
        *sp_metropolis,
        criterion_combine,
        *SP_PROTEIN_MODEL
      );

      // add the combined printer to the approximator
      approximator.SetPrinter( printer_combine);

      // start the approximation
      approximator.Approximate();

      // get the approximation result
      util::ShPtr< storage::Pair< assemble::ProteinModel, double> > sp_final_model_and_score_pair
      (
        approximator.GetTracker().GetBest()
      );
      util::ShPtr< assemble::ProteinModel> sp_model_return
      (
        util::ShPtr< assemble::ProteinModel>( approximator.GetTracker().GetBest()->First().HardCopy())
      );

      // output information regarding the final model of this stage
      BCL_MessageCrt( "#SSEs: " + util::Format()( sp_model_return->GetNumberSSEs()));

      // if this stage is the last stage print the final model
      if( m_IsLastStage)
      {
        // initialize with round number
        sp_printer->Initialize( m_RoundNumber);

        // set the prefix
        sp_printer->SetPrefix( m_Prefix);

        // print the final structure
        sp_printer->Print( approximator.GetTracker());

        s_prev_dropped_functions.Reset();
      }

      if( m_NumberOfScoresDropped)
      {
        for( size_t i( 0); i < m_NumberOfScoresDropped; ++i)
        {
          sum->SetCoefficient( dropped_functions( i), dropped_values( i));
        }
      }

      return sp_model_return;
    }

    //! @brief informs the stage that it is the last stage in the approximation progress
    void Stage::SetIsLastStage()
    {
      m_IsLastStage = true;
    }

    //! @brief reset certain members using the round number and stage number passed
    //! @param ROUND_NUMBER round of optimization
    //! @param STAGE_NUMBER stage of optimization
    void Stage::Stage::Reset( const size_t ROUND_NUMBER, const size_t STAGE_NUMBER)
    {
      // reset the temperature
      m_Temperature->Reset();

      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // call reset
        ( **protocol_itr)->Reset();
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize all score
    void Stage::InitializeScores()
    {
      // iterate over the list
      for
      (
        storage::List< fold::Protocol>::iterator itr( m_ScoreProtocols.Begin()), itr_end( m_ScoreProtocols.End());
        itr != itr_end; ++itr
      )
      {
        // initialize scores for this protocol
        ( **itr)->InitializeScores();
      }
    }

    //! @brief initialize all mutates
    void Stage::InitializeMutates()
    {
      // iterate over the list
      for
      (
        storage::List< fold::Protocol>::iterator itr( m_MutateProtocols.Begin()), itr_end( m_MutateProtocols.End());
        itr != itr_end; ++itr
      )
      {
        // initialize mutates for this protocol
        ( **itr)->InitializeMutates();
      }
    }

    void Stage::InitializePool( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const
    {
      // initialize empty pool
      util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool());

      // initialize map to hold the min pool lengths
      storage::Map< biol::SSType, size_t> min_pool_sse_lengths( assemble::SSEPool::GetCommandLineMinSSELengths());

      // open pool file
      io::IFStream read;
      const std::string pool_file
      (
        assemble::SSEPool::GetFlagPoolPrefix()->GetFirstParameter()->GetValue() + "." + GetPoolPostfix()
      );
      BCL_MessageStd( "Reading pool from file for stage setting " + pool_file);
      io::File::MustOpenIFStream( read, pool_file);

      // read pool
      sp_pool->ReadSSEPool
      (
        read,
        *SP_PROTEIN_MODEL,
        min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
        min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
      );

      // if separate pool flag was provided and pool is not overlapping
      if( fold::DefaultFlags::GetFlagPoolSeparate()->GetFlag() && !sp_pool->IsOverlapping())
      {
        // separate pools
        BCL_MessageStd( "separating adjoining SSEs in the pool");
        const bool success
        (
          sp_pool->Separate
          (
            assemble::SSEPool::GetCommandLineMinSSELengths(),
            fold::DefaultFlags::GetFlagPoolSeparate()->GetFirstParameter()->GetNumericalValue< size_t>()
          )
        );
        // make sure to prune it to remove loops
        sp_pool->Prune( assemble::SSEPool::GetCommandLineMinSSELengths());

        BCL_MessageStd( "separating success: " + util::Format()( success));
      }

      BCL_MessageStd( "pool set to ");
      sp_pool->WriteSSEPool( util::GetLogger());

      SP_PROTEIN_MODEL->SetSSEPoolData( sp_pool);
    }

    //! @brief modify the starting model
    //! @param START_MODEL starting model to be used
    void Stage::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the start model
        ( **protocol_itr)->ModifyStartModel( START_MODEL);
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    void Stage::ModifyCriterion( opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the terminate
        ( **protocol_itr)->ModifyCriterion( CRITERION, *this);
      }
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    void Stage::ModifyPrinter( PrinterCombined< assemble::ProteinModel, double> &PRINTER) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the printer
        ( **protocol_itr)->ModifyPrinter( PRINTER, *this);
      }
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to modify
    void Stage::ModifyFactory( util::ShPtr< pdb::Factory> &FACTORY) const
    {
      // reset factory printers
      FACTORY->ResetPrinters();

      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the factory
        ( **protocol_itr)->ModifyFactory( FACTORY, *this);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Stage::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Stage::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_temperature_accepted.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureAccepted::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureAccepted())
    );

    //! @brief return command line flag for setting temperature adjustment based on accepted steps ratio
    //! @return command line flag for setting temperature adjustment based on accepted steps ratio
    util::ShPtr< command::FlagInterface> &TemperatureAccepted::GetFlagTemperature()
    {
      // initialize static instance of the flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mc_temperature_fraction",
          "\tmodify the start and end fractions of accepted steps for temperature adjustments for monte carlo minimization"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( GetParameterStartFraction());
        flag->PushBack( GetParameterEndFraction());
        flag->PushBack( GetParameterStartTemperature());
        flag->PushBack( GetParameterUpdateInterval());
      }

      // end
      return s_flag;
    }

    //! @brief return command line parameter for setting the start fraction
    //! @return command line parameter for setting the start fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterStartFraction()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "start_ratio",
          "\tthis a fraction of accepted steps at the start of the algorithm",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.5"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterEndFraction()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "end_ratio",
          "\tthis a fraction of accepted steps at the end of the algorithm",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.2"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterStartTemperature()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "start_temperature", "\tthis the starting temperature ", "1"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterUpdateInterval()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "update_interval", "\tnumber of steps between each update", "10"
        )
      );

      // end
      return s_parameter;
    }

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureAccepted::TemperatureAccepted() :
      m_Temperature(),
      m_StartFraction(),
      m_EndFraction(),
      m_StartTemperature(),
      m_Delta(),
      m_UpdateInterval()
    {
    }

    //! @brief constructor from number iterations, the other paramaters are taken from the commandline arguments
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureAccepted::TemperatureAccepted
    (
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_Temperature( 0, GetParameterStartTemperature()->GetNumericalValue< double>()),
      m_StartFraction( GetParameterStartFraction()->GetNumericalValue< double>()),
      m_EndFraction( GetParameterEndFraction()->GetNumericalValue< double>()),
      m_StartTemperature( GetParameterStartTemperature()->GetNumericalValue< double>()),
      m_Delta( 0.0),
      m_UpdateInterval( GetParameterUpdateInterval()->GetNumericalValue< size_t>())
    {
    }

    //! @brief constructor from starting and ending temperature
    //! @param START_FRACTION starting fraction
    //! @param END_FRACTION ending fraction
    //! @param NUMBER_OF_ITERATIONS number of iterations
    //! @param START_TEMPERATURE starting temperature
    //! @param UPDATE_INTERVAL interval length between each update
    TemperatureAccepted::TemperatureAccepted
    (
      const double START_FRACTION,
      const double END_FRACTION,
      const size_t NUMBER_OF_ITERATIONS,
      const double START_TEMPERATURE,
      const size_t UPDATE_INTERVAL
    ) :
      m_Temperature( 0, START_TEMPERATURE),
      m_StartFraction( START_FRACTION),
      m_EndFraction( END_FRACTION),
      m_StartTemperature( START_TEMPERATURE),
      m_Delta( 0.0),
      m_UpdateInterval( UPDATE_INTERVAL)
    {
      //check that start fraction is in interval (0, 1], end fraction in interval [0, 1) and that start > end
      if( m_StartFraction > 0.0 && m_StartFraction <= 1.0 && m_EndFraction > 0.0 && m_EndFraction <= 1.0)
      {
        // calculate the delta as in the linear case
        // the /2 is required since not the correct fraction but the overall fraction of accepted steps is monitored
        m_Delta = ( m_EndFraction - m_StartFraction) / NUMBER_OF_ITERATIONS;
      }
    }

    //! @brief virtual copy constructor
    TemperatureAccepted *TemperatureAccepted::Clone() const
    {
      return new TemperatureAccepted( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureAccepted::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureAccepted::GetAlias() const
    {
      static const std::string s_name( "AcceptedTemperatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureAccepted::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Temperature control that adjusts the temperature to achieve a given fraction of accepted versus"
        "rejected mutations."
      );
      serializer.AddInitializer
      (
        "start fraction",
        "start fraction of accepted versus rejected mutations",
        io::Serialization::GetAgent( &m_StartFraction),
        "0.5"
      );
      serializer.AddInitializer
      (
        "end fraction",
        "end fraction of accepted versus rejected mutations",
        io::Serialization::GetAgent( &m_EndFraction),
        "0.2"
      );
      serializer.AddInitializer
      (
        "start temperature",
        "temperature at the beginning of the simulation",
        io::Serialization::GetAgent( &m_StartTemperature),
        "1"
      );
      serializer.AddInitializer
      (
        "delta",
        "magnitude of temperature change",
        io::Serialization::GetAgent( &m_Delta),
        "0.00015"
      );
      serializer.AddInitializer
      (
        "update interval",
        "length of the interval between updates",
        io::Serialization::GetAgent( &m_UpdateInterval),
        "10"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset temperature
    void TemperatureAccepted::Reset()
    {
      // reset members
      m_Temperature = storage::Pair< size_t, double>( 0, m_StartTemperature);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return current temperature
    //! @param TRACKER the current tracker
    //! @return current temperature
    double TemperatureAccepted::GetTemperature( const opti::TrackerBase &TRACKER) const
    {
      // if the iteration number does not match the counter
      if( TRACKER.GetIteration() - m_Temperature.First() >= m_UpdateInterval)
      {
        // update temperature
        UpdateTemperature( TRACKER);
      }

      // return temperature
      return m_Temperature.Second();
    }

    //! @brief calculate new temperature
    //! @param TRACKER the base tracker
    void TemperatureAccepted::UpdateTemperature( const opti::TrackerBase &TRACKER) const
    {
      // update only every m_ConstantIntervalLength' step
      if( TRACKER.GetIteration() - m_Temperature.First() < m_UpdateInterval)
      {
        return;
      }

      // update the iteration number
      m_Temperature.First() = TRACKER.GetIteration();

      // calculate target fraction
      const double target_fraction( m_StartFraction + m_Temperature.First() * m_Delta);

      const size_t points_to_consider
      (
        std::min( size_t( std::max( 1.0 / std::max( target_fraction, 1.0e-6), 10.0)), m_UpdateInterval)
      );
      const size_t nth_point_desired
      (
        std::max( size_t( target_fraction * points_to_consider), size_t( 1)) - 1
      );
      // copy the vector of deltas under consideration
      storage::Vector< double> deltas_copy( m_PreviousDeltas.End() - points_to_consider, m_PreviousDeltas.End());
      std::nth_element( deltas_copy.Begin(), deltas_copy.Begin() + nth_point_desired, deltas_copy.End());
      static const double ln2( std::log( 2.0));
      m_Temperature.Second() = 0.7 * m_Temperature.Second() + 0.3 * deltas_copy( nth_point_desired) / ln2;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureAccepted::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Temperature, ISTREAM);
      io::Serialize::Read( m_StartFraction, ISTREAM);
      io::Serialize::Read( m_EndFraction, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);
      io::Serialize::Read( m_Delta, ISTREAM);
      io::Serialize::Read( m_UpdateInterval, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureAccepted::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Temperature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StartFraction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EndFraction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Delta, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UpdateInterval, OSTREAM, INDENT);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_temperature_default.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureDefault::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureDefault())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureDefault::TemperatureDefault() :
      m_Temperature( util::GetUndefinedDouble())
    {
    }

    //! @brief constructor from a starting and ending temperature
    //! @param TEMPERATURE temperature
    TemperatureDefault::TemperatureDefault( const double TEMPERATURE) :
      m_Temperature( TEMPERATURE)
    {
    }

    //! @brief virtual copy constructor
    TemperatureDefault *TemperatureDefault::Clone() const
    {
      return new TemperatureDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureDefault::GetAlias() const
    {
      static const std::string s_alias( "DefaultTemperatureControl");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureDefault::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Keeps the temperature constant at the defined value.");
      serializer.AddInitializer
      (
        "temperature",
        "temperature is kept constant at this value",
        io::Serialization::GetAgent( &m_Temperature)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset temperature
    void TemperatureDefault::Reset()
    {
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureDefault::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Temperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Temperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_temperature_exponential.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureExponential::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureExponential())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureExponential::TemperatureExponential() :
      m_LastIterationNumber( 0),
      m_Scale( 0),
      m_StartTemperature( 0)
    {
    }

    //! @brief constructor from a starting and ending temperature and total number of iterations
    //! @param START_TEMPERATURE starting temperature
    //! @param END_TEMPERATURE ending temperature
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureExponential::TemperatureExponential
    (
      const double START_TEMPERATURE,
      const double END_TEMPERATURE,
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_LastIterationNumber( 0),
      m_Scale( 1.0),
      m_StartTemperature( START_TEMPERATURE)
    {
      // if the starting temperature is not 0 and the start is larger than the end
      if( START_TEMPERATURE != 0.0 && START_TEMPERATURE > END_TEMPERATURE)
      {
        // calculate the scale
        m_Scale = pow( END_TEMPERATURE / START_TEMPERATURE, 1.0 / double( NUMBER_OF_ITERATIONS));
      }
    }

    //! @brief virtual copy constructor
    TemperatureExponential *TemperatureExponential::Clone() const
    {
      return new TemperatureExponential( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureExponential::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureExponential::GetAlias() const
    {
      static const std::string s_name( "ExponentialTemparatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureExponential::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Temperature control that allows exponential change of the simulated temperature."
      );
      serializer.AddInitializer
      (
        "scaling factor for temperature adjustment",
        "scaling factor for the temperature adjustment",
        io::Serialization::GetAgent( &m_Scale)
      );
      serializer.AddInitializer
      (
        "starting temperature",
        "temperature at the start of the simulation",
        io::Serialization::GetAgent( &m_StartTemperature)
      );
      serializer.AddDataMember
      (
        "last iteration number",
        io::Serialization::GetAgent( &m_LastIterationNumber)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset this temperature
    void TemperatureExponential::Reset()
    {
      // reset members
      m_LastIterationNumber = 0;
    }

    //! @brief returns last calculated temperature without updating
    //! @return last calculated temperature without updating
    double TemperatureExponential::GetLastCalculatedTemperature() const
    {
      return m_StartTemperature * math::Pow( m_Scale, double( m_LastIterationNumber));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureExponential::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LastIterationNumber, ISTREAM);
      io::Serialize::Read( m_Scale, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureExponential::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LastIterationNumber, OSTREAM) << '\n';
      io::Serialize::Write( m_Scale, OSTREAM) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_temperature_linear.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureLinear::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureLinear())
    );

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureLinear::TemperatureLinear() :
      m_LastIterationNumber(),
      m_Delta(),
      m_StartTemperature()
    {
    }

    //! @brief constructor from a counter, starting and ending temperature and total number of iterations
    //! @param ITERATION_COUNTER iteration counter
    //! @param START_TEMPERATURE starting temperature
    //! @param END_TEMPERATURE ending temperature
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureLinear::TemperatureLinear
    (
      const double START_TEMPERATURE,
      const double END_TEMPERATURE,
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_LastIterationNumber( 0),
      m_Delta( ( END_TEMPERATURE - START_TEMPERATURE) / NUMBER_OF_ITERATIONS),
      m_StartTemperature( START_TEMPERATURE)
    {
    }

    //! @brief virtual copy constructor
    TemperatureLinear *TemperatureLinear::Clone() const
    {
      return new TemperatureLinear( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureLinear::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns last calculated temperature without updating
    //! @return last calculated temperature without updating
    double TemperatureLinear::GetLastCalculatedTemperature() const
    {
      return m_StartTemperature + m_LastIterationNumber * m_Delta;
    }

    //! @brief return current temperature
    //! @param TRACKER the current tracker
    //! @return current temperature
    double TemperatureLinear::GetTemperature( const opti::TrackerBase &TRACKER) const
    {
      // update the last iteration counter
      m_LastIterationNumber = TRACKER.GetIteration();

      // return temperature
      return GetLastCalculatedTemperature();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureLinear::GetAlias() const
    {
      static const std::string s_name( "LinearTemparatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureLinear::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Temperature control that allows linear adjustment of the temperature.");
      serializer.AddInitializer
      (
        "step size of temperature change",
        "delta by which the temperature is changed",
        io::Serialization::GetAgent( &m_Delta)
      );
      serializer.AddInitializer
      (
        "starting temperature",
        "temperature at the beginning of the simulation",
        io::Serialization::GetAgent( &m_StartTemperature)
      );
      serializer.AddDataMember
      (
        "last iteration number",
        io::Serialization::GetAgent( &m_LastIterationNumber)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset this temperature
    void TemperatureLinear::Reset()
    {
      // reset members
      m_LastIterationNumber = 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureLinear::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_LastIterationNumber, ISTREAM);
      io::Serialize::Read( m_Delta, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureLinear::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_LastIterationNumber, OSTREAM) << '\n';
      io::Serialize::Write( m_Delta, OSTREAM) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
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
#include "mc/bcl_mc_template_instantiations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    template class BCL_API PrintInterface< assemble::ProteinModel, double>;

  } // namespace mc
} // namespace bcl
