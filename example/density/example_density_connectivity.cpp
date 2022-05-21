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
#include "density/bcl_density_connectivity.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_connectivity.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityConnectivity :
    public ExampleInterface
  {
  public:

    ExampleDensityConnectivity *Clone() const
    { return new ExampleDensityConnectivity( *this);}

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

    /////////////////
    // preparation //
    /////////////////

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));

      BCL_MessageStd( "reading pdb: " + pdb_filename);

      //set parameters
      const double resolution( 5.0), voxelsize( 1.5);
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
        )
      );

      //read atoms from pdb
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_filename);
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);

      //extract all atoms for conversion to densitymap
      const util::ShPtrList< pdb::Line> atom_lines( pdb.GetLines( pdb::GetLineTypes().ATOM));
      util::ShPtrVector< biol::Atom> atoms;
      for
      (
        util::ShPtrList< pdb::Line>::const_iterator
          line_itr( atom_lines.Begin()), line_itr_end( atom_lines.End());
        line_itr != line_itr_end; ++line_itr
      )
      {
        atoms.PushBack( factory.AtomFromLine( **line_itr));
      }

      // density map calculated from SimpleAtoms, with target resolution, voxelsize and smoothingkernel of choice
      const density::Map density_map( simulator->operator ()( atoms));

      const std::string output_filename_mrc
      (
        AddExampleOutputPathToFilename( density::Connectivity(), "1IE9.mrc")
      );

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, output_filename_mrc);
      density_map.WriteMRC( write);
      io::File::CloseClearFStream( write);

      // generate a protein model from the pdb file
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      const assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      // initialize SiPtrList of all the bodies of helices in the protein
      util::SiPtrVector< const assemble::SSE> vector_of_sses( model.GetSSEs());
      util::SiPtrVector< const assemble::SSEGeometryInterface> vector_of_bodies( model.GetSSEs());

      // convert that SiPtrVector into a ShPtrVector
      util::ShPtrVector< assemble::SSEGeometryInterface> shptr_vector_of_bodies;
      for
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator vector_itr( vector_of_bodies.Begin()),
          vector_itr_end( vector_of_bodies.End());
        vector_itr != vector_itr_end; ++vector_itr
      )
      {
        shptr_vector_of_bodies.PushBack
        (
          util::ShPtr< assemble::SSEGeometryInterface>
          (
            ( *vector_itr)->Clone()
          )
        );
      }

      // initialize a ShPtrList from the ShPtrVector
      const util::ShPtrList< assemble::SSEGeometryInterface> list_of_bodies
      (
        shptr_vector_of_bodies.Begin(),
        shptr_vector_of_bodies.End()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test DensityConnectivity default constructor
      density::Connectivity connectivity_object_default;

      // write out connectivity_object_default
      BCL_MessageStd
      (
        "writing out default constructed connectivity object" +
        util::Format()( connectivity_object_default)
      );

      // test DensityConnectivity constructor from various arguments
      density::Connectivity connectivity_object
      (
        util::ShPtr< assemble::SSEGeometryInterface>( vector_of_sses( 0)->Clone()),
        true,
        util::ShPtr< assemble::SSEGeometryInterface>( vector_of_sses( 1)->Clone()),
        false,
        5.5,
        10.28
      );

    /////////////////
    // data access //
    /////////////////

      // check GetConnectivity function
      BCL_ExampleCheck( connectivity_object.GetConnectivity(), 5.5);

      // check GetDistance function
      BCL_ExampleCheck( connectivity_object.GetDistance(), 10.28);

      // check GetOrientationA function
      BCL_ExampleAssert( connectivity_object.GetOrientationA(), true);

      // check GetOrientationB function
      BCL_ExampleAssert( connectivity_object.GetOrientationB(), false);

      // check whether Body A was read in correctly
      BCL_ExampleCheck( *connectivity_object.GetBodyA(), *vector_of_sses( 0));

      // check whether Body B was read in correctly
      BCL_ExampleCheck( *connectivity_object.GetBodyB(), *vector_of_sses( 1));

    //////////////////////
    // input and output //
    //////////////////////

      // test the write and read
      WriteBCLObject( connectivity_object);
      density::Connectivity connectivity_object_to_be_read;
      ReadBCLObject( connectivity_object_to_be_read);

      // check GetConnectivity function
      BCL_ExampleCheck( connectivity_object_to_be_read.GetConnectivity(), connectivity_object.GetConnectivity());

      // check GetDistance function
      BCL_ExampleCheck( connectivity_object_to_be_read.GetDistance(), connectivity_object.GetDistance());

      // check GetOrientationA function
      BCL_ExampleCheck( connectivity_object_to_be_read.GetOrientationA(), connectivity_object.GetOrientationA());

      // check whether Body A was read in correctly
      BCL_ExampleIndirectCheck
      (
        coord::EqualWithinTolerance( *connectivity_object_to_be_read.GetBodyA(), *connectivity_object.GetBodyA()),
        true,
        "Body A of read object " + util::Format()( *connectivity_object_to_be_read.GetBodyA()) +
        " should be identical with first body in connectivity_object " +
        util::Format()( *connectivity_object.GetBodyA())
      );

      // check whether Body B was read in correctly
      BCL_ExampleIndirectCheck
      (
        coord::EqualWithinTolerance( *connectivity_object_to_be_read.GetBodyB(), *connectivity_object.GetBodyB()),
        true,
        "Body B of read object " + util::Format()( *connectivity_object_to_be_read.GetBodyB()) +
        " should be identical with second body in connectivity_object " +
        util::Format()( *connectivity_object.GetBodyB())
      );

    ////////////////
    // operations //
    ////////////////

      // create list of connectivities from density map, list of sses
      storage::List< density::Connectivity> connectivities_1IE9
      (
        density::Connectivity::DetermineConnectivities( density_map, list_of_bodies)
      );

      BCL_MessageStd( "connectivities 1IE9 based on entire helix bodies from list of bodies");
      for
      (
        storage::List< density::Connectivity>::const_iterator
          list_itr( connectivities_1IE9.Begin()), list_itr_end( connectivities_1IE9.End());
        list_itr != list_itr_end; ++list_itr
      )
      {
        BCL_MessageStd
        (
          "connectivity: " + util::Format()( list_itr->GetConnectivity()) +
            " , distance: " + util::Format()( list_itr->GetDistance())
        );
      }

      // calculate all connectivities by iterating over all pairs of helices and passing entire helix bodies
      BCL_MessageStd( "connectivities 1IE9 based on entire helix bodies ");
      size_t helix_counter_a( 0);
      size_t helix_counter_b( 0);
      storage::List< density::Connectivity> connectivities_1IE9_individual;

      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator iter_a( shptr_vector_of_bodies.Begin()),
          iter_a_end( shptr_vector_of_bodies.End());
        iter_a != iter_a_end; ++iter_a
      )
      {
        helix_counter_b = helix_counter_a + 1;
        for
        (
          util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator iter_b( iter_a + 1);
          iter_b != iter_a_end; ++iter_b
        )
        {
          connectivities_1IE9_individual =
            density::Connectivity::DetermineConnectivities( density_map, *iter_a, *iter_b);

          // iterate through the 4 resulting connectivities
          for
          (
            storage::List< density::Connectivity>::const_iterator
              connectivity_iter( connectivities_1IE9_individual.Begin()),
              connectivity_iter_end( connectivities_1IE9_individual.End());
            connectivity_iter != connectivity_iter_end;
            ++connectivity_iter
          )
          {
            if( connectivity_iter->GetDistance() < 12.0)
            {
              BCL_MessageStd
              (
                util::Format()( helix_counter_a) + " " + util::Format()( helix_counter_b) + " " +
                  util::Format()( connectivity_iter->GetOrientationA()) + " " +
                  util::Format()( connectivity_iter->GetOrientationB()) + " " +
                  util::Format()( connectivity_iter->GetConnectivity()) + " " +
                  util::Format()( connectivity_iter->GetDistance())
              );
            }
          }
          ++helix_counter_b;
        }

        ++helix_counter_a;
      }

      // connectivities_1IE9_individual always contains 4 connectivities, gets overridden for every pair of bodies
      // but the last 4 connectivities that it contains should also be in connectivities_1IE9
      // (which contains 4 connectivities for all pairs of bodies)
      // check whether one (the first) of these connectivities in connectivities_1IE9_individual is also in
      // connectivities_1IE9
      bool connectivity_found( false);
      for
      (
        storage::List< density::Connectivity>::const_iterator list_itr( connectivities_1IE9.Begin()),
          list_itr_end( connectivities_1IE9.End()); list_itr != list_itr_end; ++list_itr
      )
      {
        if( DensityConnectivityEqualWithinTolerance( *list_itr, connectivities_1IE9_individual.FirstElement()))
        {
          connectivity_found = true;
        }
      }

      BCL_ExampleIndirectCheck( connectivity_found, true, "element is not in list");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

    static bool DensityConnectivityEqualWithinTolerance
    (
      const density::Connectivity &CONNECTIVITY_LHS, const density::Connectivity &CONNECTIVITY_RHS
    )
    {
      return
      (
        coord::EqualWithinTolerance( *CONNECTIVITY_LHS.GetBodyA(), *CONNECTIVITY_RHS.GetBodyA()) &&
        coord::EqualWithinTolerance( *CONNECTIVITY_LHS.GetBodyB(), *CONNECTIVITY_RHS.GetBodyB()) &&
        CONNECTIVITY_LHS.GetOrientationA() == CONNECTIVITY_RHS.GetOrientationA() &&
        CONNECTIVITY_LHS.GetOrientationB() == CONNECTIVITY_RHS.GetOrientationB() &&
        math::EqualWithinTolerance( CONNECTIVITY_LHS.GetConnectivity(), CONNECTIVITY_RHS.GetConnectivity()) &&
        math::EqualWithinTolerance( CONNECTIVITY_LHS.GetDistance(), CONNECTIVITY_RHS.GetDistance())
      );
    }

  }; // end ExampleDensityConnectivity

  const ExampleClass::EnumType ExampleDensityConnectivity::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityConnectivity())
  );

} // namespace bcl
