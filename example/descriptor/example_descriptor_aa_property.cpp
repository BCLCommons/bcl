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
#include "descriptor/bcl_descriptor_aa_property.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_property.cpp
  //!
  //! @author teixeipl
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAProperty :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAProperty *Clone() const
    {
      return new ExampleDescriptorAAProperty( *this);
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
      // form an AA and set the blast profile
      biol::AA amino_acid( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( amino_acid);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AAProperty property( biol::AATypeData::e_Volume);

      // clone
      util::ShPtr< descriptor::AAProperty> sp_property( property.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( property.GetSizeOfFeatures(), 1);

      // sequence to chain to protein model to protein model with cache
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( util::CloneToShPtr( aa_sequence)));
      assemble::ProteinModelWithCache protein_model_with_cache( assemble::ProteinModel( sp_chain), false);

      descriptor::Iterator< biol::AABase> itr( property.GetType());
      itr.SetObject( protein_model_with_cache);
      property.SetObject( protein_model_with_cache);

      linal::Vector< float> volume_result( 1);
      volume_result( 0) = amino_acid.GetType()->GetAAProperty( biol::AATypeData::e_Volume);

      // check the descriptions produced
      BCL_ExampleCheck( property( itr), volume_result);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAAProperty

  const ExampleClass::EnumType ExampleDescriptorAAProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAProperty())
  );
  
} // namespace bcl
