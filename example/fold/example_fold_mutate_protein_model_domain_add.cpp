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
#include "fold/bcl_fold_mutate_protein_model_domain_add.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_all_possible_domains.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_domain_add.cpp
  //! @details checks the mutate and writes out two PDBs
  //!
  //! @author weinerbe
  //! @date Nov 8, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelDomainAdd :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelDomainAdd *Clone() const
    {
      return new ExampleFoldMutateProteinModelDomainAdd( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test the default constructor
      fold::MutateProteinModelDomainAdd def_construct;

      // test constructor from a domain
      assemble::CollectorAllPossibleDomains domain_collector( math::Range< size_t>( 3, 3), true);
      fold::MutateProteinModelDomainAdd domain_construct( domain_collector);

      // test constructor from a domain and bool to use a fold template
      fold::MutateProteinModelDomainAdd domain_template_construct( domain_collector, true);

      // test Clone
      util::ShPtr< fold::MutateProteinModelDomainAdd> clone_construct( domain_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      const std::string expected_identifier( "bcl::fold::MutateProteinModelDomainAdd");
      BCL_ExampleCheck( clone_construct->GetClassIdentifier(), expected_identifier);

      // test GetScheme
      BCL_ExampleCheck( domain_construct.GetScheme(), expected_identifier);

    ///////////////
    // operators //
    ///////////////

      // test operator() and write results to file
      const size_t nr_sses( 5);
      assemble::ProteinModel mutated_model_a( *( domain_construct( protein_model).GetArgument()));
      Proteins::WriteModelToPDB
      (
        mutated_model_a, AddExampleOutputPathToFilename( domain_construct, "1ubi_mutate_domain_add_a.pdb")
      );
      BCL_ExampleCheck( mutated_model_a.GetNumberSSEs(), nr_sses);

      assemble::ProteinModel mutated_model_b( *( domain_template_construct( protein_model).GetArgument()));
      BCL_ExampleCheck( mutated_model_b.GetNumberSSEs(), nr_sses);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( domain_construct);

      // read the object back in
      fold::MutateProteinModelDomainAdd read_construct;
      ReadBCLObject( read_construct);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelDomainAdd

  const ExampleClass::EnumType ExampleFoldMutateProteinModelDomainAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelDomainAdd())
  );

} // namespace bcl
