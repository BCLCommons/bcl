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
#include "biol/bcl_biol_dssp.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_dssp.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolDssp :
    public ExampleInterface
  {
  public:

    ExampleBiolDssp *Clone() const
    {
      return new ExampleBiolDssp( *this);
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
      // initialize
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      biol::DSSP dssp;

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const math::MutateResult< assemble::ProteinModel> result( dssp( model));

    //////////////////////
    // input and output //
    //////////////////////

      const std::string pdb_dssp_filename( AddExampleOutputPathToFilename( dssp, "1ubi_dssp.pdb"));
      Proteins::WriteModelToPDB( *result.GetArgument(), pdb_dssp_filename);

      const std::string out_filename( AddExampleOutputPathToFilename( dssp, "1ubi.dssp"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, out_filename);
      dssp.WriteToFile( write);
      io::File::CloseClearFStream( write);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolDssp

  const ExampleClass::EnumType ExampleBiolDssp::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolDssp())
  );

} // namespace bcl
