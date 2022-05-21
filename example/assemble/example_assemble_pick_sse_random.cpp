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
#include "assemble/bcl_assemble_pick_sse_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_pick_sse_random.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010 
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePickSSERandom :
    public ExampleInterface
  {
  public:

    ExampleAssemblePickSSERandom *Clone() const
    {
      return new ExampleAssemblePickSSERandom( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create collector
      assemble::CollectorSSE collector_sse;

      // test default constructor
      assemble::PickSSERandom pick_sse_random;

      // test Pick function
      BCL_MessageStd( "test Pick function");
      util::SiPtr< const assemble::SSE> sp_sse
      (
        pick_sse_random.Pick( collector_sse.Collect( protein_model))
      );

      BCL_MessageStd( "randomly choosen sse: " + sp_sse->GetIdentification());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssemblePickSSERandom

  const ExampleClass::EnumType ExampleAssemblePickSSERandom::s_Instance
  (
    GetExamples().AddEnum( ExampleAssemblePickSSERandom())
  );

} // namespace bcl
