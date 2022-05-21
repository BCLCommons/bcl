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
#include "assemble/bcl_assemble_collector_common_aa.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_common_aa.cpp
  //!
  //! @author alexanns
  //! @date Apr 22, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorCommonAA :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorCommonAA *Clone() const
    {
      return new ExampleAssembleCollectorCommonAA( *this);
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
      // create string "pdb_filename_a" which has path for example pdb file
      const std::string pdb_filename_a( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create string "pdb_filename_b" which has path for example pdb file
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "1ubi_unpaired_strand.pdb"));

      // create ProteinModel "protein_model_a" from "pdb_a"
      BCL_MessageStd( "building model_a from pdb_a chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model_a( Proteins::GetModel( pdb_filename_a, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // create ProteinModel "protein_model_a" from "pdb_a"
      BCL_MessageStd( "building model_b from pdb_a chains and sse information");
      assemble::ProteinModel protein_model_b( Proteins::GetModel( pdb_filename_b, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // default constructor
      assemble::CollectorCommonAA collector;

      // test Collect function
      BCL_MessageStd( "test Collect function");

      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > common_aas;

      for( int i = 0; i < 50; ++i)
      {
        common_aas =
        (
          collector.Collect( storage::VectorND< 2, assemble::ProteinModel>( protein_model_a, protein_model_b))
        );
      }

      // make sure that the correct number of amino acids are collected
      const size_t expected_number_collected_amino_acids( 42);
      BCL_Example_Check
      (
        common_aas.First().GetSize() == expected_number_collected_amino_acids
        && common_aas.Second().GetSize() == expected_number_collected_amino_acids,
        "number of amino acids collected is " + util::Format()( common_aas.First().GetSize())
        + " and " + util::Format()( common_aas.Second().GetSize()) + " but should be "
        + util::Format()( expected_number_collected_amino_acids)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorCommonAA

  const ExampleClass::EnumType ExampleAssembleCollectorCommonAA::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorCommonAA())
  );

} // namespace bcl
