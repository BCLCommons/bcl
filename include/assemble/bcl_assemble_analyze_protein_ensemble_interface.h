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

#ifndef BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_INTERFACE_H_
#define BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeProteinEnsembleInterface
    //! @brief interface for classes that can analyze a protein ensemble and output the analysis
    //! @details The analysis can be text, tables, scripts to create graphs, etc.
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date Jul 30, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeProteinEnsembleInterface :
      public util::FunctionInterface< ProteinEnsemble, std::string>,
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! @brief return command line flag for specifying the prefix prepended to each analysis' postfix
      //! @return command line flag for specifying the prefix prepended to each analysis' postfix
      static const util::ShPtr< command::FlagInterface> &GetFlagOutFilePrefix();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual AnalyzeProteinEnsembleInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      virtual const std::string &GetOutFilePostfix() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief function to create the file that contains the analysis information
      //! @param OUT_FILE_PREFIX the string that will be prepended to the postfix to create the full filename
      //! @param ENSEMBLE the ensemble that will be analyzed
      void WriteAnalysisFile( const std::string &OUT_FILE_PREFIX, const ProteinEnsemble &ENSEMBLE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class AnalyzeProteinEnsembleInterface

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_INTERFACE_H_
