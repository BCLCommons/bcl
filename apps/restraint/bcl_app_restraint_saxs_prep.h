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

#ifndef BCL_APP_RESTRAINT_SAXS_PREP_H_
#define BCL_APP_RESTRAINT_SAXS_PREP_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintSaxsPrep
    //! @brief RestraintSaxPreps is for preparing experimental SAXS data for analysis.  This provides the functionality
    //!        to simulate experimental error, reduce data points to avoid overfitting, and to provide input format for
    //!        folding
    //!
    //! @see @link example_app_restraint_saxs_prep.cpp @endlink
    //! @author putnamdk
    //! @date 07/25/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSaxsPrep :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag for experimental restraint filename
      util::ShPtr< command::FlagInterface> m_ExpRestraintFileFlag;

      //! flag for pdb file from which Shannon Sampling will be evaluated
      util::ShPtr< command::FlagInterface> m_PDBFilenameFlag;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputProteinModelFlag;

      //! flag to determine the norm factor with regula falsi (true) or pythagorean approximation (false)
      util::ShPtr< command::FlagInterface> m_UseAnalyticNormFactorFlag;

      //! Parameter for max dimension of the protein
      util::ShPtr< command::FlagInterface> m_DmaxFlag;

      //! Parameter for max dimension of the protein
      util::ShPtr< command::FlagInterface> m_GnomeFlag;

      //! Parameter for the number of iterations in Shannon Sampling
      util::ShPtr< command::FlagInterface> m_SamplingRoundsFlag;

      //! flag to simulate errors
      util::ShPtr< command::FlagInterface> m_SimulateErrorFlag;

      //! pass flag if you want to score the file
      util::ShPtr< command::FlagInterface> m_ScoreFileFlag;

      //! flag to simulate errors
      util::ShPtr< command::FlagInterface> m_ReadGnomeFitFlag;

      //! flag to set experimental error to provided value
      util::ShPtr< command::FlagInterface> m_SetErrorFlag;

      //! flag for Reducing DataSet using Shannon Sampling
      util::ShPtr< command::FlagInterface> m_ReduceDataFlag;

      //! flag for Reducing DataSet using min error Sampling
      util::ShPtr< command::FlagInterface> m_ReduceDataMinErrorFlag;

      //! flag to use experimental error
      util::ShPtr< command::FlagInterface> m_UseErrorFlag;

      //! flag for creating data file for analysis in R with BCL Data
      util::ShPtr< command::FlagInterface> m_BCLFileNameFlag;

      //! flag for creating data file for analysis in R with Crysol Data
      util::ShPtr< command::FlagInterface> m_FoxsFileNameFlag;

      //! flag for creating data file for analysis in R with Foxs Data
      util::ShPtr< command::FlagInterface> m_CrysolFileNameFlag;

      //! pass flag if you want to convert profiles to log scale for visualization
      util::ShPtr< command::FlagInterface> m_LogScaleFlag;

      //! pass flag if you want to convert profiles to log scale for visualization
      util::ShPtr< command::FlagInterface> m_DerivativeFlag;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSaxsPrep();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      RestraintSaxsPrep *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name

      const std::string &GetClassIdentifier() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      static const ApplicationType RestraintSaxsPrep_Instance;

    }; // RestraintSaxsPrep
  } // namespace app
} // namespace bcl

#endif // BCL_APP_RESTRAINT_SAXS_PREP_H_
