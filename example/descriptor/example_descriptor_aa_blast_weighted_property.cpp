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
#include "descriptor/bcl_descriptor_aa_blast_weighted_property.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_blast_weighted_property.cpp
  //!
  //! @author mendenjl
  //! @date Mar 07, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAABlastWeightedProperty :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAABlastWeightedProperty *Clone() const
    {
      return new ExampleDescriptorAABlastWeightedProperty( *this);
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
      util::ShPtr< biol::AAData> aa_data_ptr( new biol::AAData( biol::GetAATypes().ALA, 1, 1));
      biol::AA amino_acid( aa_data_ptr);
      linal::Vector< double> aa_profile( biol::AATypes::s_NumberStandardAATypes, double( 1.0));

      // set the blast profile to be all 1's
      amino_acid.SetBlastProfile( biol::BlastProfile( aa_profile));
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( amino_acid);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AABlastWeightedProperty property( biol::AATypeData::e_Hydrophobicity);

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

      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance( property( itr), linal::MakeVector< float>( 0.0241), 0.001);

      // now change the blast profile to something more typical for a hydrophobic residue

      // a blast profile for a very hydrophobic residue
      storage::Vector< double> hydrophobic_vector
      (
        util::SplitStringToNumerical< double>( "-3  -5  -5  -4  4 -5  -5  -5  -5  2 3 -5  2 5 -3  -2  -2  5 0 -1", " ")
      );
      protein_model_with_cache.ResetCache();
      itr.SetObject( protein_model_with_cache);
      property.SetObject( protein_model_with_cache);
      ( *aa_sequence.Begin())->SetBlastProfile
      (
        linal::Vector< double>( hydrophobic_vector.Begin(), hydrophobic_vector.End())
      );

      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance( property( itr), linal::MakeVector< float>( 3.348), 0.001);

      // a blast profile for a very hydrophilic residue
      storage::Vector< double> hydrophilic_vector
      (
        util::SplitStringToNumerical< double>( "-4 -1 2 7 -6 -2 1 -1 -3 -6 -6 -3 -5 -6 -1 -1 -3 -6 -5 -5", " ")
      );
      ( *aa_sequence.Begin())->SetBlastProfile
      (
        linal::Vector< double>( hydrophilic_vector.Begin(), hydrophilic_vector.End())
      );
      protein_model_with_cache.ResetCache();
      itr.SetObject( protein_model_with_cache);
      property.SetObject( protein_model_with_cache);
      // note that the values for hydrophobicity (-2.546) is lower than the minimum hydrophobicity of any particular
      // amino acid (which is around -1.01).  This is because the blast profile has given us additional information
      // (e.g. this residue strongly favors hydrophilic types and strongly disfavors hydrophobic types) that suggest
      // suggest that the effective hydrophobicity is much lower
      BCL_ExampleCheckWithinTolerance( property( itr), linal::MakeVector< float>( -2.546), 0.001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAABlastWeightedProperty

  const ExampleClass::EnumType ExampleDescriptorAABlastWeightedProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAABlastWeightedProperty())
  );

} // namespace bcl
