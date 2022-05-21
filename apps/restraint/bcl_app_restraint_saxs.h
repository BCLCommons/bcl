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

#ifndef BCL_APP_RESTRAINT_SAXS_H_
#define BCL_APP_RESTRAINT_SAXS_H_

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
    //! @class RestraintSaxs
    //! @brief RestraintSaxs is for simulating saxs data from a given protein model and corresponding
    //!        experimental data.
    //!
    //! @see @link example_app_restraint_saxs.cpp @endlink
    //! @author putnamdk, weinerbe, loweew, heinzes1
    //! @date 03/16/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSaxs :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! flag for pdb file from which the restraints will be simulated
      util::ShPtr< command::FlagInterface> m_PDBFilenameFlag;

      //! flag for experimental restraint filename
      util::ShPtr< command::FlagInterface> m_ExpRestraintFileFlag;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputFileFlag;

      //! flag for Solvent Accessible Surface Area file
      util::ShPtr< command::FlagInterface> m_SasaFileFlag;

      //! flag for printing RMSD score
      util::ShPtr< command::FlagInterface> m_PrintRMSD;

      //! flag for printing Stovgaard score
      util::ShPtr< command::FlagInterface> m_PrintStovgaard;

      //! flag for printing Stovgaard score
      util::ShPtr< command::FlagInterface> m_PrintProteinModelWithLoopsFileFlag;

      //! flag for printing Derivative
      util::ShPtr< command::FlagInterface> m_PrintDerivative;

      //! flag for Reducing DataSet
      util::ShPtr< command::FlagInterface> m_ReduceData;

      //! flag for using Simulated/Experimental Errors in Chi calculation
      util::ShPtr< command::FlagInterface> m_UseErrors;

      //! flag to approximate side chains
      util::ShPtr< command::FlagInterface> m_ApproximateSideChains;

      //! flag to approximate loops
      util::ShPtr< command::FlagInterface> m_ApproximateLoops;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSaxs();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      RestraintSaxs *Clone() const;

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

      static const ApplicationType RestraintSaxs_Instance;

    }; // RestraintSaxs
  } // namespace app
} // namespace bcl

#endif // BCL_APP_RESTRAINT_SAXS_H_
