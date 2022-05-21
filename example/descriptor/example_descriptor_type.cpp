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
#include "descriptor/bcl_descriptor_type.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_type.cpp
  //! @details this example demonstrates how different Type will simulate a density map from a list of atoms
  //!
  //! @author mendenjl
  //! @date Dec 04, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorType :
    public ExampleInterface
  {
  public:

    ExampleDescriptorType *Clone() const
    {
      return new ExampleDescriptorType( *this);
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

      // create a scalar, elementwise, and pairwise types
      const descriptor::Type scalar( 0, true, descriptor::Type::e_Symmetric);
      const descriptor::Type element( 1, true, descriptor::Type::e_Symmetric);
      const descriptor::Type pairwise_asymmetric_repeats( 2, true, descriptor::Type::e_Asymmetric);
      const descriptor::Type triplet_symmetric_no_repeats( 3, false, descriptor::Type::e_Symmetric);
      const descriptor::Type triplet_symmetric_repeats( 3, true, descriptor::Type::e_Symmetric);
      const descriptor::Type triplet_asymmetric_no_repeats( 3, false, descriptor::Type::e_Asymmetric);

      // check some of the enums
      BCL_ExampleCheck( scalar.GetDimension(), 0);
      BCL_ExampleCheck( element.GetDimension(), 1);
      BCL_ExampleCheck( pairwise_asymmetric_repeats.GetDimension(), 2);
      BCL_ExampleCheck( triplet_symmetric_no_repeats.GetDimension(), 3);
      BCL_ExampleCheck( scalar.GetSymmetry(), descriptor::Type::e_Symmetric);
      BCL_ExampleCheck( element.GetSymmetry(), descriptor::Type::e_Symmetric);
      BCL_ExampleCheck( pairwise_asymmetric_repeats.GetSymmetry(), descriptor::Type::e_Asymmetric);

      // try out GetNumberFeatures
      BCL_ExampleCheck( scalar.GetNumberFeatures( 100), 1);
      BCL_ExampleCheck( element.GetNumberFeatures( 100), 100);
      BCL_ExampleCheck( pairwise_asymmetric_repeats.GetNumberFeatures( 100), 100 * 100);
      BCL_ExampleCheck( triplet_symmetric_no_repeats.GetNumberFeatures( 100), 100 * 99 * 98 / 6);
      BCL_ExampleCheck( triplet_symmetric_repeats.GetNumberFeatures( 100), 102 * 101 * 100 / 6);
      BCL_ExampleCheck( triplet_asymmetric_no_repeats.GetNumberFeatures( 100), 100 * 99 * 98);

      // try out generalize
      descriptor::Type pairwise_asymmetric_repeats_generalized( pairwise_asymmetric_repeats);
      pairwise_asymmetric_repeats_generalized.GeneralizeToHandle( triplet_symmetric_no_repeats);
      BCL_ExampleIndirectCheck
      (
        pairwise_asymmetric_repeats_generalized.GetDimension(),
        3,
        "GeneralizeToHandle"
      );
      BCL_ExampleIndirectCheck
      (
        pairwise_asymmetric_repeats_generalized.GetSymmetry(),
        descriptor::Type::e_Asymmetric,
        "GeneralizeToHandle"
      );
      BCL_ExampleIndirectCheck
      (
        pairwise_asymmetric_repeats_generalized.ConsiderRepeatedObjects(),
        true,
        "GeneralizeToHandle"
      );

      storage::Vector< size_t> one_two_three( storage::Vector< size_t>::Create( 1, 2, 3));

      BCL_ExampleCheck( descriptor::Type( 3, false, descriptor::Type::e_Symmetric).GetPosition( one_two_three, 5), 6);
      BCL_ExampleCheck( descriptor::Type( 3, false, descriptor::Type::e_Asymmetric).GetPosition( one_two_three, 5), 16);
      BCL_ExampleCheck( descriptor::Type( 3, true,  descriptor::Type::e_Asymmetric).GetPosition( one_two_three, 5), 38);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleDescriptorType

  const ExampleClass::EnumType ExampleDescriptorType::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorType())
  );

} // namespace bcl
