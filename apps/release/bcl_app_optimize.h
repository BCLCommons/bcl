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

#ifndef BCL_APP_OPTIMIZE_H_
#define BCL_APP_OPTIMIZE_H_

// include the namespace header
#include "app/bcl_app.h"

// include other forward headers - sorted alphabetically
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Optimize
    //! @brief Optimizes a set of protein models or predicts an ensemble of conformations.
    //!
    //! @see @link example_app_optimize.cpp @endlink
    //! @author fischea
    //! @date Oct 31, 2016
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Optimize :
      public Interface
    {

    ///////////
    // data //
    ///////////

    public:

      //! single instance of this class
      static const ApplicationType Optimizes_Instance;

    private:

      //! flag for a list of PDBs to be optimized
      util::ShPtr< command::FlagStatic> m_PDBList;

      //! flag for the prefix of the output files
      util::ShPtr< command::FlagStatic> m_OutputPrefix;

      //! flag for the path to the loop construction pipeline
      util::ShPtr< command::FlagStatic> m_OptimizerFileName;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      Optimize();

      //! @brief clone function
      //! @return pointer to a new Optimize
      Optimize *Clone() const;

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

    }; // class Optimize

  } // namespace app
} // namespace bcl

#endif // BCL_APP_OPTIMIZE_H_
