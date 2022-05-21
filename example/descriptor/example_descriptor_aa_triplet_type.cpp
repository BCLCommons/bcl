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
#include "descriptor/bcl_descriptor_aa_triplet_type.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_triplet_type.cpp
  //!
  //! @author mendenjl
  //! @date Apr 17, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAATripletType :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAATripletType *Clone() const
    {
      return new ExampleDescriptorAATripletType( *this);
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
      biol::AA alanine( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));
      biol::AA phenylalanine( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().PHE, 2, 2)));
      biol::AA arginine( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ARG, 3, 3)));
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( alanine);
      aa_sequence.PushBack( phenylalanine);
      aa_sequence.PushBack( arginine);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AATripletType aatriplettype;

      // clone
      util::ShPtr< descriptor::AATripletType> sp_blast( aatriplettype.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( aatriplettype.GetSizeOfFeatures(), 8000);

    ////////////////
    // operations //
    ////////////////

      // sequence to chain to protein model to protein model with cache
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( util::CloneToShPtr( aa_sequence)));
      assemble::ProteinModelWithCache protein_model_with_cache( assemble::ProteinModel( sp_chain), false);

      descriptor::Iterator< biol::AABase> itr( aatriplettype.GetType());
      itr.SetObject( protein_model_with_cache);
      aatriplettype.SetObject( protein_model_with_cache);

      // goto the middle AA
      ++itr;

      // check the descriptions produced
      BCL_ExampleCheckWithinAbsTolerance( aatriplettype( itr).Sum(), 1.0, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( aatriplettype( itr).Max(), 1.0, 0.001);
      BCL_ExampleCheckWithinAbsTolerance
      (
        aatriplettype( itr)( biol::GetAATypes().PHE.GetIndex() + biol::GetAATypes().ALA.GetIndex() * 20 + biol::GetAATypes().ARG.GetIndex() * 400),
        1.0,
        0.001
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAATripletType

  const ExampleClass::EnumType ExampleDescriptorAATripletType::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAATripletType())
  );

} // namespace bcl
