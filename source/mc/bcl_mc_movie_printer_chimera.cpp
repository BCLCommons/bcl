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
