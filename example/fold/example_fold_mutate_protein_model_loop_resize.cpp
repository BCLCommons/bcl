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
#include "fold/bcl_fold_mutate_protein_model_loop_resize.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse_terminus_residue.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_mutate_protein_model_grow_sse.h"
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_loop_resize.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Sep 6, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelLoopResize :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelLoopResize *Clone() const
    {
      return new ExampleFoldMutateProteinModelLoopResize( *this);
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

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 999;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      protein_model.AddLoops( true, false);

      util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > > phi_psi_generator
      (
        new fold::PhiPsiGeneratorRamachandran( biol::Ramachandran::GetDefaultHistogramFilename())
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::MutateProteinModelLoopResize def_constr;

      // parameter taking member variable parameters
      // tested for correctness below
      fold::MutateProteinModelLoopResize param_constr
      (
        util::ShPtr< assemble::LocatorSSETerminusResidue>
        (
          new assemble::LocatorSSETerminusResidue( 'A', 35, biol::AASequenceFlexibility::e_NTerminal)
        ),
        0.0,
        math::Range< size_t>( 4, 4),
        storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
        storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
        biol::AASequenceFlexibility::e_CTerminal,
        phi_psi_generator,
        true,
        "scheme"
      );

      // clone constructor
      util::ShPtr< fold::MutateProteinModelLoopResize> clone_constr( param_constr.Clone());
      BCL_ExampleCheck( clone_constr->GetScheme(), "scheme");

    /////////////////
    // data access //
    /////////////////

      // GetScheme
      BCL_ExampleCheck( param_constr.GetScheme(), "scheme");

    ///////////////
    // operators //
    ///////////////

      // test c terminal operations
      BCL_MessageDbg( "test c terminal operations");
      {
        util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > locator
        (
          new assemble::LocatorSSETerminusResidue( 'A', 35, biol::AASequenceFlexibility::e_NTerminal)
        );

        fold::MutateProteinModelLoopResize shrink
        (
          locator,
          0.0,
          math::Range< size_t>( 4, 4),
          storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
          storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
          biol::AASequenceFlexibility::e_CTerminal,
          phi_psi_generator
        );

        math::MutateResult< assemble::ProteinModel> mutated_model_result( shrink( protein_model));

        BCL_ExampleAssert( mutated_model_result.GetArgument().IsDefined(), true);

        const assemble::ProteinModel &mutated_model( *mutated_model_result.GetArgument());

        const util::SiPtr< const assemble::SSE> shrunk_sse( locator->Locate( mutated_model));

        BCL_ExampleAssert( shrunk_sse.IsDefined(), true);

        BCL_MessageDbg( "shrunk sse is " + util::Format()( shrunk_sse->GetIdentification()));
        BCL_ExampleCheck( shrunk_sse->GetIdentification(), "COIL A   35 GLY <==>   35 GLY");

        // give the residue coordinates
        util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > anchor_locator
        (
          new assemble::LocatorAA( 'A', 34)
        );
        fold::MutateProteinModelGrowSSE grower
        (
          locator, phi_psi_generator, anchor_locator, biol::AASequenceFlexibility::e_CTerminal
        );

        math::MutateResult< assemble::ProteinModel> grow_result( grower( mutated_model));

        BCL_ExampleAssert( grow_result.GetArgument().IsDefined(), true);

        const assemble::ProteinModel &mutated_grown_model( *grow_result.GetArgument());

        // try to grow when already at max size
        {
          BCL_MessageDbg( "try to grow when already at max size");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_CTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutated_model( resizer( protein_model));

          BCL_ExampleCheck( mutated_model.GetArgument().IsDefined(), false);
        }

        // add one residue
        {
          BCL_MessageDbg( "test operator add one residue");
          BCL_Assert( locator.IsDefined(), "locator is not defined");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_CTerminal,
            phi_psi_generator,
            true,
            "scheme"
          );

          math::MutateResult< assemble::ProteinModel> mutate_result( resizer( mutated_grown_model));

          BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

          const util::SiPtr< const assemble::SSE> new_sse( locator->Locate( *mutate_result.GetArgument()));

          BCL_ExampleAssert( new_sse.IsDefined(), true);

          BCL_MessageDbg( "resized sse is " + util::Format()( new_sse->GetIdentification()));
          BCL_ExampleCheck( new_sse->GetIdentification(), "COIL A   35 GLY <==>   36 ILE");
          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
          {
            const std::string pdb_file
            (
              AddExampleOutputPathToFilename( resizer, "MutateProteinModelLoopResize_a.pdb")
            );
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, pdb_file);

            pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
          }

          assemble::LocatorAA locator_35( 'A', 35);
          assemble::LocatorAA locator_36( 'A', 36);
          const util::SiPtr< const biol::AABase> resi_35( locator_35.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_35.IsDefined(), true);
          const util::SiPtr< const biol::AABase> resi_36( locator_36.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_36.IsDefined(), true);
          const double psi( resi_35->CalculatePsi( resi_36->GetAtom( biol::GetAtomTypes().N)));
          const double phi( resi_36->CalculatePhi( resi_35->GetAtom( biol::GetAtomTypes().C)));
          const double expected_psi( -0.654583);
          const double expected_phi( -1.17825);
          BCL_ExampleCheckWithinTolerance( expected_psi, psi, 0.00001);
          BCL_MessageDbg
          (
            "expected psi is " + util::Format()( expected_psi) + " calculated psi is "
            + util::Format()( psi)
          );
          BCL_ExampleCheckWithinTolerance( expected_phi, phi, 0.00001);
          BCL_MessageDbg
          (
            "expected phi is " + util::Format()( expected_phi) + " calculated phi is "
            + util::Format()( phi)
          );
        }

        // try to shrink when already at min size
        {
          BCL_MessageDbg( "try to shrink when already at min size");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            0.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_CTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutated_model( resizer( mutated_grown_model));

          BCL_ExampleCheck( mutated_model.GetArgument().IsDefined(), false);
        }

        // add more than one residue
        {
          BCL_MessageDbg( "test operator add more than one residue");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 4, 4),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_CTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutate_result( resizer( mutated_grown_model));

          BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

          const util::SiPtr< const assemble::SSE> new_sse( locator->Locate( *mutate_result.GetArgument()));

          BCL_ExampleAssert( new_sse.IsDefined(), true);

          BCL_MessageDbg( "resized sse is " + util::Format()( new_sse->GetIdentification()));
          BCL_ExampleCheck( new_sse->GetIdentification(), "COIL A   35 GLY <==>   39 ASP");
          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
          {
            const std::string pdb_file
            (
              AddExampleOutputPathToFilename( resizer, "MutateProteinModelLoopResize_b.pdb")
            );
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, pdb_file);

            pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
          }

          assemble::LocatorAA locator_37( 'A', 37);
          assemble::LocatorAA locator_38( 'A', 38);
          const util::SiPtr< const biol::AABase> resi_37( locator_37.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_37.IsDefined(), true);
          const util::SiPtr< const biol::AABase> resi_38( locator_38.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_38.IsDefined(), true);
          const double psi( resi_37->CalculatePsi( resi_38->GetAtom( biol::GetAtomTypes().N)));
          const double phi( resi_38->CalculatePhi( resi_37->GetAtom( biol::GetAtomTypes().C)));
          const double expected_psi( -0.130917);
          const double expected_phi( -1.17825);
          BCL_ExampleCheckWithinTolerance( expected_psi, psi, 0.00001);
          BCL_MessageDbg
          (
            "expected psi is " + util::Format()( expected_psi) + " calculated psi is "
            + util::Format()( psi)
          );
          BCL_ExampleCheckWithinTolerance( expected_phi, phi, 0.00001);
          BCL_MessageDbg
          (
            "expected phi is " + util::Format()( expected_phi) + " calculated phi is "
            + util::Format()( phi)
          );
        }
      } //< test c terminal operations

      // test n terminal operations
      BCL_MessageDbg( "test n terminal operations");
      {
        util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > locator
        (
          new assemble::LocatorSSETerminusResidue( 'A', 39, biol::AASequenceFlexibility::e_CTerminal)
        );

        fold::MutateProteinModelLoopResize shrink
        (
          locator,
          0.0,
          math::Range< size_t>( 4, 4),
          storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
          storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
          biol::AASequenceFlexibility::e_NTerminal,
          phi_psi_generator
        );

        math::MutateResult< assemble::ProteinModel> mutated_model_result( shrink( protein_model));

        BCL_ExampleAssert( mutated_model_result.GetArgument().IsDefined(), true);

        const assemble::ProteinModel &mutated_model( *mutated_model_result.GetArgument());

        const util::SiPtr< const assemble::SSE> shrunk_sse( locator->Locate( mutated_model));

        BCL_ExampleAssert( shrunk_sse.IsDefined(), true);

        BCL_MessageDbg( "shrunk sse is " + util::Format()( shrunk_sse->GetIdentification()));
        BCL_ExampleCheck( shrunk_sse->GetIdentification(), "COIL A   39 ASP <==>   39 ASP");

        // give the residue coordinates
        util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > anchor_locator
        (
          new assemble::LocatorAA( 'A', 40)
        );
        fold::MutateProteinModelGrowSSE grower
        (
          locator, phi_psi_generator, anchor_locator, biol::AASequenceFlexibility::e_NTerminal
        );

        math::MutateResult< assemble::ProteinModel> grow_result( grower( mutated_model));

        BCL_ExampleAssert( grow_result.GetArgument().IsDefined(), true);

        const assemble::ProteinModel &mutated_grown_model( *grow_result.GetArgument());

        // try to grow when already at max size
        {
          BCL_MessageDbg( "try to grow when already at max size");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_NTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutated_model( resizer( protein_model));

          BCL_ExampleCheck( mutated_model.GetArgument().IsDefined(), false);
        }

        // add one residue
        {
          BCL_MessageDbg( "test operator add one residue");
          BCL_Assert( locator.IsDefined(), "locator is not defined");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_NTerminal,
            phi_psi_generator,
            true,
            "scheme"
          );

          math::MutateResult< assemble::ProteinModel> mutate_result( resizer( mutated_grown_model));

          BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

          const util::SiPtr< const assemble::SSE> new_sse( locator->Locate( *mutate_result.GetArgument()));

          BCL_ExampleAssert( new_sse.IsDefined(), true);

          BCL_MessageDbg( "resized sse is " + util::Format()( new_sse->GetIdentification()));
          BCL_ExampleCheck( new_sse->GetIdentification(), "COIL A   38 PRO <==>   39 ASP");
          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
          {
            const std::string pdb_file
            (
              AddExampleOutputPathToFilename( resizer, "MutateProteinModelLoopResize_c.pdb")
            );
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, pdb_file);

            pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
          }

          assemble::LocatorAA locator_38( 'A', 38);
          assemble::LocatorAA locator_39( 'A', 39);
          const util::SiPtr< const biol::AABase> resi_38( locator_38.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_38.IsDefined(), true);
          const util::SiPtr< const biol::AABase> resi_39( locator_39.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_39.IsDefined(), true);
          const double psi( resi_38->CalculatePsi( resi_39->GetAtom( biol::GetAtomTypes().N)));
          const double phi( resi_39->CalculatePhi( resi_38->GetAtom( biol::GetAtomTypes().C)));
          const double expected_psi(  2.48742);
          const double expected_phi( -0.916417);
          BCL_ExampleCheckWithinTolerance( expected_psi, psi, 0.00001);
          BCL_MessageDbg
          (
            "expected psi is " + util::Format()( expected_psi) + " calculated psi is "
            + util::Format()( psi)
          );
          BCL_ExampleCheckWithinTolerance( expected_phi, phi, 0.00001);
          BCL_MessageDbg
          (
            "expected phi is " + util::Format()( expected_phi) + " calculated phi is "
            + util::Format()( phi)
          );
        }

        // try to shrink when already at min size
        {
          BCL_MessageDbg( "try to shrink when already at min size");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            0.0,
            math::Range< size_t>( 1, 1),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_NTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutated_model( resizer( mutated_grown_model));

          BCL_ExampleCheck( mutated_model.GetArgument().IsDefined(), false);
        }

        // add more than one residue
        {
          BCL_MessageDbg( "test operator add more than one residue");
          fold::MutateProteinModelLoopResize resizer
          (
            locator,
            1.0,
            math::Range< size_t>( 4, 4),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1)),
            storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 5)),
            biol::AASequenceFlexibility::e_NTerminal,
            phi_psi_generator
          );

          math::MutateResult< assemble::ProteinModel> mutate_result( resizer( mutated_grown_model));

          BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

          const util::SiPtr< const assemble::SSE> new_sse( locator->Locate( *mutate_result.GetArgument()));

          BCL_ExampleAssert( new_sse.IsDefined(), true);

          BCL_MessageDbg( "resized sse is " + util::Format()( new_sse->GetIdentification()));
          BCL_ExampleCheck( new_sse->GetIdentification(), "COIL A   35 GLY <==>   39 ASP");
          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
          {
            const std::string pdb_file
            (
              AddExampleOutputPathToFilename( resizer, "MutateProteinModelLoopResize_d.pdb")
            );
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, pdb_file);

            pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
          }

          assemble::LocatorAA locator_37( 'A', 37);
          assemble::LocatorAA locator_38( 'A', 38);
          const util::SiPtr< const biol::AABase> resi_37( locator_37.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_37.IsDefined(), true);
          const util::SiPtr< const biol::AABase> resi_38( locator_38.Locate( *mutate_result.GetArgument()));
          BCL_ExampleAssert( resi_38.IsDefined(), true);
          const double psi( resi_37->CalculatePsi( resi_38->GetAtom( biol::GetAtomTypes().N)));
          const double phi( resi_38->CalculatePhi( resi_37->GetAtom( biol::GetAtomTypes().C)));
          const double expected_psi( -0.130917);
          const double expected_phi( -1.17825);
          BCL_ExampleCheckWithinTolerance( expected_psi, psi, 0.00001);
          BCL_MessageDbg
          (
            "expected psi is " + util::Format()( expected_psi) + " calculated psi is "
            + util::Format()( psi)
          );
          BCL_ExampleCheckWithinTolerance( expected_phi, phi, 0.00001);
          BCL_MessageDbg
          (
            "expected phi is " + util::Format()( expected_phi) + " calculated phi is "
            + util::Format()( phi)
          );
        }
      } //< test n terminal operations

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

  }; //end ExampleFoldMutateProteinModelLoopResize

  const ExampleClass::EnumType ExampleFoldMutateProteinModelLoopResize::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelLoopResize())
  );

} // namespace bcl
