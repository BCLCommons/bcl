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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "model/bcl_model_restricted_boltzmann_machine_layer.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_restricted_boltzmann_machine_layer.cpp
  //!
  //! @author mendenjl
  //! @date Aug 02, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRestrictedBoltzmannMachineLayer :
    public ExampleInterface
  {
  public:

    ExampleModelRestrictedBoltzmannMachineLayer *Clone() const
    {
      return new ExampleModelRestrictedBoltzmannMachineLayer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      static const size_t s_number_cols( 2);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a simple boltzmann network
      static const size_t s_number_hidden_neurons( 8);
      model::RestrictedBoltzmannMachineLayer boltz_23
      (
        s_number_cols,
        s_number_hidden_neurons,
        model::RestrictedBoltzmannMachineLayer::e_StochasticSigmoid
      );
      BCL_ExampleCheck( boltz_23.GetNumberInputNeurons(), s_number_cols);
      BCL_ExampleCheck( boltz_23.GetNumberHiddenNeurons(), s_number_hidden_neurons);
      BCL_ExampleCheck( boltz_23.GetVisibleBias().GetSize(), s_number_cols);
      BCL_ExampleCheck( boltz_23.GetHiddenBias().GetSize(), s_number_hidden_neurons);
      BCL_ExampleCheck( boltz_23.GetWeight().GetNumberRows(), s_number_cols);
      BCL_ExampleCheck( boltz_23.GetWeight().GetNumberCols(), s_number_hidden_neurons);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( boltz_23);
      BCL_MessageVrb( "read object");
      model::RestrictedBoltzmannMachineLayer boltz_23_read;
      ReadBCLObject( boltz_23_read);

      BCL_ExampleIndirectCheck( boltz_23_read.GetNumberInputNeurons(), s_number_cols, "I/O");
      BCL_ExampleIndirectCheck( boltz_23_read.GetNumberHiddenNeurons(), s_number_hidden_neurons, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  };
  //end ExampleModelRestrictedBoltzmannMachineLayer

  const ExampleClass::EnumType ExampleModelRestrictedBoltzmannMachineLayer::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRestrictedBoltzmannMachineLayer())
  );

} // namespace bcl
