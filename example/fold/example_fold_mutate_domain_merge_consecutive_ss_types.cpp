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
#include "fold/bcl_fold_mutate_domain_merge_consecutive_ss_types.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_add_multiple.cpp
  //!
  //! @author bitterd
  //! @date Sep 11, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateDomainMergeConsecutiveSSTypes :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateDomainMergeConsecutiveSSTypes *Clone() const
    {
      return new ExampleFoldMutateDomainMergeConsecutiveSSTypes( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1gk8.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 3;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct the mutate fromt SSType coil
      fold::MutateDomainMergeConsecutiveSSTypes mutate( biol::GetSSTypes().COIL);
      protein_model.AddLoops( true, false);

      //! get a domain out of the ProteinModel
      const assemble::Domain domain( protein_model.GetSSEsAsDomain());
      util::SiPtrVector< const assemble::SSE> vector( domain.GetSSEs());

      std::string string_original;

      BCL_MessageStd( "Output SSEs of the original domain");
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator itr( vector.Begin()),
          itr_end( vector.End()); itr != itr_end; ++itr
      )
      {
        BCL_MessageStd( ( *itr)->GetIdentification());
        string_original += ( *itr)->GetIdentification() + '\n';
      }

      std::string compare_org_out( AddExampleInputPathToFilename( e_Biology, "FoldMutateMergeConsecutiveSSTtypes_original.data"));
      std::string new_org_out( AddExampleOutputPathToFilename( mutate, "FoldMutateMergeConsecutiveSSTtypes_original.data"));

      io::OFStream write_original;

      BCL_ExampleMustOpenOutputFile( write_original, new_org_out);

      write_original << string_original;

      io::File::CloseClearFStream( write_original);

      BCL_ExampleCheck
      (
        io::File::FilesMatch( new_org_out, compare_org_out), true
      );
    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "mutating the protein model");

      assemble::Domain mutated_domain( *mutate( domain).GetArgument());
      util::SiPtrVector< const assemble::SSE> mutated_vector( mutated_domain.GetSSEs());

      BCL_MessageStd( "Output SSEs of the mutated domain");

      std::string string_mutate;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator itr( mutated_vector.Begin()),
          itr_end( mutated_vector.End()); itr != itr_end; ++itr
      )
      {
        BCL_MessageStd( ( *itr)->GetIdentification());
        string_mutate += ( *itr)->GetIdentification() + '\n';
      }

      std::string compare_mutate_out( AddExampleInputPathToFilename( e_Biology, "FoldMutateMergeConsecutiveSSTtypes_mutate.data"));
      std::string new_mutate_out( AddExampleOutputPathToFilename( mutate, "FoldMutateMergeConsecutiveSSTtypes_mutate.data"));

      io::OFStream write_mutate;

      BCL_ExampleMustOpenOutputFile( write_mutate, new_mutate_out);

      write_mutate << string_mutate;

      io::File::CloseClearFStream( write_mutate);

      BCL_ExampleCheck( io::File::FilesMatch( compare_mutate_out, new_mutate_out), true);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateDomainMergeConsecutiveSSTypes

  const ExampleClass::EnumType ExampleFoldMutateDomainMergeConsecutiveSSTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateDomainMergeConsecutiveSSTypes())
  );

} // namespace bcl
