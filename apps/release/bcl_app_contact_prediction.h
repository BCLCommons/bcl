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

#ifndef BCL_APP_CONTACT_PREDICTION_H_
#define BCL_APP_CONTACT_PREDICTION_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContactPrediction
    //! @brief Class predicts contacts for the given protein sequences
    //! @details Class is the application for predicting contacts for the given protein sequences
    //!
    //! @see @link example_app_contact_prediction.cpp @endlink
    //! @author karakam
    //! @date 03/24/2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContactPrediction :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! filename for input
      util::ShPtr< command::ParameterInterface> m_InputFilenameParam;

      //! pdb list - to be used if providing a list of pdbs
      util::ShPtr< command::FlagInterface> m_PDBListFlag;

      //! merged threshold - residue pairs with merged predictions below this threshold do not get printed out
      util::ShPtr< command::FlagInterface> m_ThresholdFlag;

      //!real contacts
      util::ShPtr< command::FlagInterface> m_RealContactsFlag;

      //! no prediction
      util::ShPtr< command::FlagInterface> m_NoPredictionFlag;

      //!output path
      util::ShPtr< command::FlagInterface> m_OutputPathFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ContactPrediction();

    public:

      //! @brief Clone function
      //! @return pointer to new ContactPrediction
      ContactPrediction *Clone() const
      {
        return new ContactPrediction( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      std::string GetBCLScopedName() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      // instantiate enumerator for ContactPrediction class
      static const ApplicationType ContactPrediction_Instance;

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief creates the contact map and real contact map( if requested) and outputs them
      //! @brief for the provided chain and pdb_tag
      //! @param SP_CHAIN ShPtr to the chain for which the contact map is going to be generated
      //! @param PDB_TAG pdb tag of the protein of interest
      void CreateContactMap
      (
        const util::ShPtr< assemble::Chain> &SP_CHAIN,
        const std::string &PDB_TAG
      ) const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

    }; // class ContactPrediction

  } // namespace app
} // namespace bcl

#endif // BCL_APP_CONTACT_PREDICTION_H_
