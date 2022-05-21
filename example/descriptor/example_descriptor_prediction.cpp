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
#include "descriptor/bcl_descriptor_prediction.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_iterator.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_prediction.cpp
  //!
  //! @author mendenjl
  //! @date Feb 12, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorPrediction :
    public ExampleInterface
  {
  public:

    ExampleDescriptorPrediction *Clone() const
    {
      return new ExampleDescriptorPrediction( *this);
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

      // get the directory containing the models
      const std::string directory( AddExampleInputPathToFilename( e_Descriptor, ""));

      // prefix, to distinguish these models from others in that directory
      const std::string prefix( "predict");

      // create an implementation to do the prediction; in this case for 2 * AASeqID + 1
      util::Implementation< descriptor::Base< biol::AABase, float> > predicter
      (
        "Prediction(storage=File(directory=" + directory + ", prefix= " + prefix + "))"
      );

      // check sizes
      BCL_ExampleCheck( predicter->GetSizeOfFeatures(), 2);

      // load in a simple model
      assemble::ProteinModelWithCache mini_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "mini.pdb")),
        true
      );

      descriptor::Iterator< biol::AABase> itr( predicter->GetType());
      itr.SetObject( mini_model);
      predicter->SetObject( mini_model);

      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance
      (
        predicter->operator ()( itr),
        linal::MakeVector< float>( 4.18846, 3.04201),
        0.0001
      );
      BCL_ExampleCheckWithinTolerance
      (
        predicter->operator ()( ++itr),
        linal::MakeVector< float>( 5.00283, 5.00434),
        0.0001
      );

      // try changing the dimensions
      itr =
        descriptor::Iterator< biol::AABase>( descriptor::Type( 2, false, descriptor::Type::e_Symmetric), mini_model);
      predicter->SetDimension( 2);
      BCL_ExampleIndirectCheckWithinTolerance
      (
        predicter->operator ()( itr),
        linal::MakeVector< float>( 4.18846, 3.04201, 5.00283, 5.00434),
        0.0001,
        "SetDimenion compatibility with prediction"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorPrediction

  const ExampleClass::EnumType ExampleDescriptorPrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorPrediction())
  );

} // namespace bcl
