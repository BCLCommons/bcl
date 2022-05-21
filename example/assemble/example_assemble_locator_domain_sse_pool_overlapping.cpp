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
#include "assemble/bcl_assemble_locator_domain_sse_pool_overlapping.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_domain_sse_pool_overlapping.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Jun 15, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorDomainSSEPoolOverlapping :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorDomainSSEPoolOverlapping *Clone() const
    {
      return new ExampleAssembleLocatorDomainSSEPoolOverlapping( *this);
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
      // create protein model
      BCL_MessageStd( "building model");
      std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "4A2N.pdb"));
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create sse pool
      BCL_MessageStd( "read .pool");
      assemble::SSEPool sse_pool;
      io::IFStream read;
      std::string pool_filename( AddExampleInputPathToFilename( e_Biology, "4A2N.SSPredMC_JUFO9D_OCTOPUS.pool"));
      BCL_ExampleMustOpenInputFile( read, pool_filename);
      sse_pool.ReadSSEPool( read, protein_model, 3, 3);
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::LocatorDomainSSEPoolOverlapping def_constr;

      // constructor taking a pool
      assemble::LocatorDomainSSEPoolOverlapping param_constr( sse_pool);
      BCL_ExampleCheck( param_constr.GetPool().GetSize(), 12);

      // clone constructor
      util::ShPtr< assemble::LocatorDomainSSEPoolOverlapping> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetPool().GetSize(), 12);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< assemble::LocatorDomainSSEPoolOverlapping>(), clone_constr->GetClassIdentifier()
      );

      // GetPool
      BCL_ExampleCheck( param_constr.GetPool().GetSize(), 12);

    ////////////////
    // operations //
    ////////////////

      // Locate
      util::ShPtr< assemble::Domain> located_domain( param_constr.Locate( protein_model));
      BCL_ExampleCheck( located_domain->GetNumberSSEs(), 6);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( *clone_constr);
      assemble::LocatorDomainSSEPoolOverlapping read_locator;
      ReadBCLObject( read_locator);
      BCL_ExampleCheck( read_locator.GetPool().GetSize(), 12);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorDomainSSEPoolOverlapping

  const ExampleClass::EnumType ExampleAssembleLocatorDomainSSEPoolOverlapping::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorDomainSSEPoolOverlapping())
  );
  
} // namespace bcl
