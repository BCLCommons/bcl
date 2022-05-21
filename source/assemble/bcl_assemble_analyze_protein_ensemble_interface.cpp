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
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! @brief return command line flag for specifying the prefix prepended to each analysis' postfix
    //! @return command line flag for specifying the prefix prepended to each analysis' postfix
    const util::ShPtr< command::FlagInterface> &AnalyzeProteinEnsembleInterface::GetFlagOutFilePrefix()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "analysis_prefix", "\tthe prefix prepended to each analysis' postfix",
          command::Parameter( "prefix", "\tthe prefix prepended to each analysis' postfix")
        )
      );
      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief function to create the file that contains the analysis information
    //! @param OUT_FILE_PREFIX the string that will be prepended to the postfix to create the full filename
    //! @param ENSEMBLE the ensemble that will be analyzed
    void AnalyzeProteinEnsembleInterface::WriteAnalysisFile
    (
      const std::string &OUT_FILE_PREFIX, const ProteinEnsemble &ENSEMBLE
    ) const
    {
      // get the output file name
      const std::string out_filename( OUT_FILE_PREFIX + GetOutFilePostfix());

      // for writing the file
      io::OFStream write;

      // open the stream to the file
      io::File::MustOpenOFStream( write, out_filename);

      // write the analysis from the operator
      write << operator()( ENSEMBLE);

      // close the filestream
      io::File::CloseClearFStream( write);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
