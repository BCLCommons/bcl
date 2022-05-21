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
#include "fold/bcl_fold_handler_locator_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_handler_locator_loop_domain.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldHandlerLocatorLoopDomain :
    public ExampleInterface
  {
  public:

    ExampleFoldHandlerLocatorLoopDomain *Clone() const
    {
      return new ExampleFoldHandlerLocatorLoopDomain( *this);
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

      {
        assemble::ProteinModel model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"))
        );

        model.AddLoops( true, false);

        fold::HandlerLocatorLoopDomain handler;

        const std::string pdb_name( AddExampleOutputPathToFilename( handler, "add_loops_model.pdb"));
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, pdb_name);

        BCL_MessageDbg( "writing to file " + pdb_name);
        pdb::Factory().WriteModelToPDB( model, write);
        io::File::CloseClearFStream( write);

        const std::string domain_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_domain_a.txt"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, domain_filename);

        BCL_MessageDbg( "creating domain locator");
        fold::LocatorLoopDomain loop_domain_locator( handler.HandleRead( read, model));
        io::File::CloseClearFStream( read);

        // check that the correct number of loop segment locators have been read in
        BCL_MessageDbg( "loop segments are " + util::Format()( loop_domain_locator.GetLoopSegments()));
        BCL_ExampleCheck
        (
          loop_domain_locator.GetLoopSegments().GetSize(), 3
        );
      }

      {
        BCL_MessageDbg( "creating domain locator with full model");
        assemble::ProteinModel model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
        );

        model.AddLoops( true, false);

        // false means don't randomize the pseudo residue phi and psi if it can be calculated from model
        fold::HandlerLocatorLoopDomain handler( false);

        const std::string domain_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_domain_a.txt"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, domain_filename);

        BCL_MessageDbg( "creating domain locator");
        fold::LocatorLoopDomain loop_domain_locator( handler.HandleRead( read, model));
        io::File::CloseClearFStream( read);

        // check that the correct number of loop segment locators have been read in
        BCL_MessageDbg( "loop segments are " + util::Format()( loop_domain_locator.GetLoopSegments()));
        BCL_ExampleCheck
        (
          loop_domain_locator.GetLoopSegments().GetSize(), 3
        );
      }

      // test when randomizing phi psi of pseudo residue even if it can be calculated from the model
      {
        BCL_MessageDbg( "randomizing phi psi of pseudo resi even if it can be calc. from model");
        assemble::ProteinModel model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
        );

        model.AddLoops( true, false);

        // default true parameter value - randomize the pseudo residue phi and psi even if it can be calc from model
        fold::HandlerLocatorLoopDomain handler;

        const std::string domain_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_domain_a.txt"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, domain_filename);

        BCL_MessageDbg( "creating domain locator");
        fold::LocatorLoopDomain loop_domain_locator( handler.HandleRead( read, model));
        io::File::CloseClearFStream( read);

        // check that the correct number of loop segment locators have been read in
        BCL_ExampleCheck
        (
          loop_domain_locator.GetLoopSegments().GetSize(), 3
        );
      }

      // test HandleReadMultipleFunction
      {
        BCL_MessageDbg( "test HandleReadMultipleFunction");
        assemble::ProteinModel model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"))
        );

        model.AddLoops( true, false);

        fold::HandlerLocatorLoopDomain handler;

        const std::string pdb_name( AddExampleOutputPathToFilename( handler, "add_loops_model.pdb"));
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, pdb_name);

        BCL_MessageDbg( "writing to file " + pdb_name);
        pdb::Factory().WriteModelToPDB( model, write);
        io::File::CloseClearFStream( write);

        const std::string domain_filename( AddExampleInputPathToFilename( e_Fold, "handler_loop_domain_multiple_2.txt"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, domain_filename);

        BCL_MessageDbg( "creating domain locator");
        util::ShPtrList< fold::LocatorLoopDomain> loop_domain_locator( handler.HandleReadMultiple( read, model));
        io::File::CloseClearFStream( read);

        // check that the correct number of domains have been read in
        BCL_ExampleCheck
        (
          loop_domain_locator.GetSize(), 2
        );

        // make sure that each domain has the correct size
        BCL_ExampleCheck
        (
          loop_domain_locator.FirstElement()->GetLoopSegments().GetSize(), 3
        );
        BCL_ExampleCheck
        (
          loop_domain_locator.LastElement()->GetLoopSegments().GetSize(), 2
        );
      }

      // test CreateLocatorLoopDomainsForCoil
      {
        BCL_MessageDbg( "test CreateLocatorLoopDomainsForCoil function with ubiqutin");

        // get a protein model with missing loops
        assemble::ProteinModel model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"))
        );

        // add loops to the model
        model.AddLoops( true, false);

        // create a handler
        fold::HandlerLocatorLoopDomain handler;

        // get the loop domain locators created from "model"
        const util::ShPtrList< fold::LocatorLoopDomain> locators
        (
          handler.CreateLocatorLoopDomainsForInteriorCoil( model)
        );

        // make sure the correct number of loop domain locators were created
        BCL_ExampleCheck
        (
          locators.GetSize(), 6
        );

        // make sure the all the loop domains have one loop segment in them and it is not rigid
        for
        (
          util::ShPtrList< fold::LocatorLoopDomain>::const_iterator
            domain_itr( locators.Begin()), domain_itr_end( locators.End());
          domain_itr != domain_itr_end;
          ++domain_itr
        )
        {
          // check one segment
          BCL_ExampleCheck( ( *domain_itr)->GetLoopSegments().GetSize(), 1);

          // check segment is not rigid
          BCL_ExampleCheck( ( *domain_itr)->GetLoopSegments().Begin()->IsRigid(), false);

          // make sure the segment actually exists in the protein model by trying to locate it
          const util::SiPtr< const assemble::SSE> located_sse
          (
            ( *domain_itr)->GetLoopSegments().Begin()->GetLocatorSSE().Locate( model)
          );

          // make sure "" is defined indicating that the sse was found in model
          BCL_ExampleCheck( located_sse.IsDefined(), true);
        }
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

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

  }; //end ExampleFoldHandlerLocatorLoopDomain

  const ExampleClass::EnumType ExampleFoldHandlerLocatorLoopDomain::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldHandlerLocatorLoopDomain())
  );

} // namespace bcl
