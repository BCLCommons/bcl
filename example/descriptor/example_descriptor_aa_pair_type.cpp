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
#include "descriptor/bcl_descriptor_aa_pair_type.h"

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
  //! @example example_descriptor_aa_pair_type.cpp
  //!
  //! @author mendenjl
  //! @date Oct 25, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAPairType :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAPairType *Clone() const
    {
      return new ExampleDescriptorAAPairType( *this);
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
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( alanine);
      aa_sequence.PushBack( phenylalanine);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AAPairType aapairtype;

      // clone
      util::ShPtr< descriptor::AAPairType> sp_blast( aapairtype.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( aapairtype.GetSizeOfFeatures(), 400);

    ////////////////
    // operations //
    ////////////////

      // sequence to chain to protein model to protein model with cache
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( util::CloneToShPtr( aa_sequence)));
      assemble::ProteinModelWithCache protein_model_with_cache( assemble::ProteinModel( sp_chain), false);

      descriptor::Iterator< biol::AABase> itr( aapairtype.GetType());
      itr.SetObject( protein_model_with_cache);
      aapairtype.SetObject( protein_model_with_cache);

      // check the descriptions produced
      BCL_ExampleCheckWithinAbsTolerance( aapairtype( itr).Sum(), 2.0, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( aapairtype( itr)( 13), 1.05, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAAPairType

  const ExampleClass::EnumType ExampleDescriptorAAPairType::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAPairType())
  );

} // namespace bcl
