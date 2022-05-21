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

#ifndef BCL_APP_LOOP_TEMPLATE_H_
#define BCL_APP_LOOP_TEMPLATE_H_

// include the namespace header
#include "app/bcl_app.h"

// include other forward headers - sorted alphabetically
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_loop_parameters.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopTemplate
    //! @brief Creates a library of loop conformations from a given set of pdbs
    //!
    //! @see @link example_app_loop_template.cpp @endlink
    //! @author fischea
    //! @date Dec 16, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LoopTemplate :
      public Interface
    {

    ///////////
    // data //
    ///////////

    public:

      //! single instance of this class
      static const ApplicationType LoopTemplates_Instance;

    private:

      //! flag for a list of pdbs to create the library from
      util::ShPtr< command::FlagStatic> m_PDBList;

      //! flag for output name of the template library
      util::ShPtr< command::FlagStatic> m_TemplateLibraryFile;

      //! flag for output name of the statistics file
      util::ShPtr< command::FlagStatic> m_StatisticsFile;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      LoopTemplate();

      //! @brief clone function
      //! @return pointer to a new LoopTemplate
      LoopTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns a shared pointer to the command object
      //! @return shared pointer to the command object
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the main function of this application
      //! @return exit code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads members from a given input stream
      //! @param ISTREAM input stream to read members from
      //! @return input stream which members were read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into a given output stream
      //! @param OSTREAM output stream to write members into
      //! @param INDENT number of indentations
      //! @return output stream into which members were written
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! @brief writes out the results file
      //! @param PARAMETERS parameters of the gathered loop conformations
      //! @param FILE_NAME name of the output file
      static void WriteResults
      (
        const storage::Vector< fold::LoopParameters> &PARAMETERS, const std::string &FILE_NAME
      );

      //! @brief writes out statistics about the loops
      //! @param PARAMETERS parameters of the gathered loop conformations
      //! @param FILE_NAME name of the output file
      static void WriteStatistics
      (
        const storage::Vector< fold::LoopParameters> &PARAMETERS, const std::string &FILE_NAME
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief returns the loop conformations in the given protein model
      //! @detail loop conformations are parameterized through their dihedral angles
      //! @param MODEL protein model to get loop conformations from
      //! @return the loop conformations in the given protein model
      static storage::Vector< fold::LoopParameters> GetLoops( const assemble::ProteinModel &MODEL);

      //! @brief checks if the given loop parameters are valid
      //! @param LOOP loop parameters to be checked
      //! @return true, if the given loop parameters are valid
      static bool IsValid( const fold::LoopParameters &LOOP);

    }; // class LoopTemplate

  } // namespace app
} // namespace bcl

#endif // BCL_APP_LOOP_TEMPLATE_H_
