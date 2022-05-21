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
#include "descriptor/bcl_descriptor_mutation_aa_property.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_mutations.h"
#include "biol/bcl_biol_protein_mutation_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_mutation_aa_property.cpp
  //!
  //! @author mendenjl
  //! @date Jan 22, 2019
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMutationAAProperty :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMutationAAProperty *Clone() const
    {
      return new ExampleDescriptorMutationAAProperty( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2ic8A.pdb"));

      // create protein model from pdb
      biol::ProteinMutationSet protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete),
        true,
        storage::Vector< biol::Mutation>::Create
        (
          biol::Mutation( 1, biol::GetAATypes().GLU, biol::GetAATypes().ALA),
          biol::Mutation( 2, biol::GetAATypes().ARG, biol::GetAATypes().ALA),
          biol::Mutation( 3, biol::GetAATypes().ALA, biol::GetAATypes().VAL)
        )
      );

      // create an iterator on the protein model
      descriptor::Iterator< biol::Mutation> itr_mut
      (
        descriptor::Type( 1, false, descriptor::Type::e_Symmetric),
        protein_model
      );

      // goto the 1st mutation in the protein
      itr_mut.GotoPosition( 0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor with a property that will calculate the norm charge for each residue, considering side chains and
      // backbone
      descriptor::MutationAAProperty native_aa_type
      (
        true,
        util::Implementation< descriptor::Base< biol::AABase, float> >( "ReflectingWindow(size=1,MaxIndex( AAType),alignment=JufoCenter)")
      );
      descriptor::MutationAAProperty mutant_aa_type
      (
        false,
        util::Implementation< descriptor::Base< biol::AABase, float> >( "ReflectingWindow(size=1,MaxIndex( AAType),alignment=JufoCenter)")
      );
      native_aa_type.SetObject( protein_model);
      mutant_aa_type.SetObject( protein_model);

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( native_aa_type.GetClassIdentifier(), GetStaticClassName( native_aa_type));

      // test the GetLength function
      BCL_ExampleCheck( native_aa_type.GetSizeOfFeatures(), 3);
      BCL_ExampleCheck( native_aa_type.GetType().GetDimension(), 1);
      BCL_ExampleCheck( native_aa_type.GetAlias(), "NativeAADescriptor");

    ///////////////
    // operators //
    ///////////////

      // test the operator
      BCL_MessageStd
      (
        "aa type native/mutant of " + itr_mut( 0)->ToString() + " = "
        + util::Format()( native_aa_type( itr_mut)) + " "
        + util::Format()( mutant_aa_type( itr_mut))
      );

      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 1.0, 6.0, 1.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 1.0, 0.0, 1.0), 0.0001);

      ++itr_mut;
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 6.0, 1.0, 0.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 6.0, 0.0, 0.0), 0.0001);

      ++itr_mut;
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 1.0, 0.0, 7.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 1.0, 19.0, 7.0), 0.0001);

      itr_mut.GotoPosition( 0);
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 1.0, 0.0, 1.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 1.0, 6.0, 1.0), 0.0001);

      ++itr_mut;
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 6.0, 0.0, 0.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 6.0, 1.0, 0.0), 0.0001);

      ++itr_mut;
      BCL_ExampleCheckWithinTolerance( mutant_aa_type( itr_mut), linal::Vector3D( 1.0, 19.0, 7.0), 0.0001);
      BCL_ExampleCheckWithinTolerance( native_aa_type( itr_mut), linal::Vector3D( 1.0, 0.0, 7.0), 0.0001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMutationAAProperty

  const ExampleClass::EnumType ExampleDescriptorMutationAAProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMutationAAProperty())
  );

} // namespace bcl
