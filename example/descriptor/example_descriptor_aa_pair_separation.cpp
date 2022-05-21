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
#include "descriptor/bcl_descriptor_aa_pair_separation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_pair_separation.cpp
  //!
  //! @author teixeipl, mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAPairSeparation :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAPairSeparation *Clone() const
    {
      return new ExampleDescriptorAAPairSeparation( *this);
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
      // form an AA and set the aa_separation profile
      biol::AA amino_acid( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AAPairSeparation aa_separation;

      // clone
      util::ShPtr< descriptor::AAPairSeparation> sp_aa_separation( aa_separation.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( aa_separation.GetSizeOfFeatures(), 1);

    ////////////////
    // operations //
    ////////////////

      // form an AA and set the blast profile
      biol::AA amino_acid_a( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));
      biol::AA amino_acid_b( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().TYR, 2, 1)));
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( amino_acid_a);
      aa_sequence.PushBack( amino_acid_b);

      // sequence to chain to protein model to protein model with cache
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( util::CloneToShPtr( aa_sequence)));
      assemble::ProteinModelWithCache protein_model_with_cache( assemble::ProteinModel( sp_chain), false);

      descriptor::Iterator< biol::AABase> itr( aa_separation.GetType());
      itr.SetObject( protein_model_with_cache);
      aa_separation.SetObject( protein_model_with_cache);

      linal::Vector< float> separation_result( 1);
      separation_result( 0) = 0;

      // check the descriptions produced
      BCL_ExampleCheck( aa_separation( itr), separation_result);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAAPairSeparation

  const ExampleClass::EnumType ExampleDescriptorAAPairSeparation::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAPairSeparation())
  );

} // namespace bcl
