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

#ifndef BCL_APP_MODEL_TEST_H_
#define BCL_APP_MODEL_TEST_H_

// include header of this class
#include "app/bcl_app_apps.h"

// include forward headers
#include "descriptor/bcl_descriptor.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ModelTest
    //! @brief Application for testing a model (model::Interface e.g., NeuralNetwork, SVM). Testing a model::Interface
    //!        trained in the Approximator framework for different quality measurements through evaluation of
    //!        objective functions.
    //!
    //! @see @link example_app_model_test.cpp @endlink
    //! @author butkiem1
    //! @date 06/08/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ModelTest :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      // virtual screening or model evaluation
      util::ShPtr< command::FlagInterface> m_FlagOutputBase;

      // true to only output the average prediction, not the individual predictions
      util::ShPtr< command::FlagInterface> m_FlagAverageOnly;

      // dataset finder
      util::ShPtr< command::FlagInterface> m_FlagDataSetRetriever;

      //! flag for storing the final trained model
      util::ShPtr< command::FlagInterface> m_FlagModelInterfaceStorage;

      //! flag for the list of id labels
      util::ShPtr< command::FlagInterface> m_FlagIDCode;

    ///////////////////////////////////
    // construction and destruction  //
    ///////////////////////////////////

      //! default constructor
      ModelTest();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ModelTest *Clone() const
      {
        return new ModelTest( *this);
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

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////////
    // helper methods //
    ////////////////////

    private:

      //! @brief write out model predictions and experimental data
      //! @param DATA dataset with experimental data and basis of predictions
      //! @param PREDICTIONS predicted data
      //! @param FILENAME filename for written out predictions
      void WritePredictions
      (
        const descriptor::Dataset &DATA,
        const linal::MatrixConstInterface< float> &PREDICTIONS,
        const std::string &FILENAME
      ) const;

    public:

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

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

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType ModelTest_Instance;

    }; // ModelTest

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MODEL_TEST_H_
