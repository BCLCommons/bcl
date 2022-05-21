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
#include "fold/bcl_fold_add_parabolic_loops.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_add_parabolic_loops.cpp
  //! @details compare the calculated intensity with experimental intensity from CRYSOL for ubiquitin
  //!
  //! @author putnamdk, mendenjl
  //! @date Aug 28, 2012
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldAddParabolicLoops :
    public ExampleInterface
  {
  public:

    ExampleFoldAddParabolicLoops *Clone() const
    {
      return new ExampleFoldAddParabolicLoops( *this);
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

      fold::AddParabolicLoops added( false);
      assemble::ProteinModel model_with_loops( added( protein_model));
      pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
      io::OFStream write;
      const std::string out_filename( AddExampleOutputPathToFilename( model_with_loops, "1C1D_apxloops.pdb"));

      BCL_ExampleMustOpenOutputFile( write, out_filename);

      factory.WriteModelToPDB( model_with_loops, write);
      io::File::CloseClearFStream( write);
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldAddParabolicLoops

  const ExampleClass::EnumType ExampleFoldAddParabolicLoops::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldAddParabolicLoops())
  );

} // namespace bcl
