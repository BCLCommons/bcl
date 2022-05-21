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
#include "score/bcl_score_body_connectivity_density.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "restraint/bcl_restraint_handler_body.h"

// external includes - sorted alphabetically

namespace bcl
{
  // explicit template instantiation
  template class util::ShPtrVector< assemble::SSEGeometryInterface>;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_body_connectivity_density.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreBodyConnectivityDensity :
    public ExampleInterface
  {
  public:

    ExampleScoreBodyConnectivityDensity *Clone() const
    { return new ExampleScoreBodyConnectivityDensity( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
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
      const util::ShPtr< density::SimulateInterface> simulator( density::GetSimulators().CreateSimulator( density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution));

      //read atoms from pdb
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_filename);
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      // create "factory" to create protein model with amino acids of type AABackBone
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
        AddExampleOutputPathToFilename( score::BodyConnectivityDensity(), "1IE9.mrc")
      );

      BCL_MessageStd( "writing mrc: " + output_filename_mrc);

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, output_filename_mrc);
      density_map.WriteMRC( write);
      io::File::CloseClearFStream( write);

      // generate a protein model from the pdb file
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 12;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      const assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      // initialize SiPtrList of all the bodies of helices in the protein
      util::SiPtrVector< const assemble::SSE> vector_of_bodies( model.GetSSEs());

      // initialize empty ShPtrVector of SSEGeometryInterfaces
      util::ShPtrVector< assemble::SSEGeometryInterface> shptr_vector_of_bodies;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator vector_itr( vector_of_bodies.Begin()),
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

      const util::ShPtrList< assemble::SSEGeometryInterface> list_of_bodies
      (
        shptr_vector_of_bodies.Begin(),
        shptr_vector_of_bodies.End()
      );

      // create a ShPtrVector of bodies which will be created from "model_sses" and will be used to create the
      // restraint::Body
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > bodies( new util::ShPtrVector< assemble::SSEGeometryInterface>());

      // create a ShPtrVector of just two bodies which will be created from the first and second sse of "model_sses"
      // and will be used to create the small restraint::Body
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > two_bodies( new util::ShPtrVector< assemble::SSEGeometryInterface>());

      // add the coord::Bodies of "model_sses" to "bodies"
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          sse_itr( shptr_vector_of_bodies.Begin()), sse_itr_end( shptr_vector_of_bodies.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        bodies->PushBack( *sse_itr);
      }

      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > coord_bodies
      (
        new util::ShPtrVector< assemble::SSEGeometryInterface>( *bodies)
      );

      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > coord_two_bodies
      (
        new util::ShPtrVector< assemble::SSEGeometryInterface>( *two_bodies)
      );

      // push back first and second sse from vector_of_bodies into two_bodies
      util::ShPtr< assemble::SSEGeometryInterface> temp_first_sse( ( vector_of_bodies( 0))->Clone());
      util::ShPtr< assemble::SSEGeometryInterface> temp_second_sse( ( vector_of_bodies( 1))->Clone());
      two_bodies->PushBack( temp_first_sse);
      two_bodies->PushBack( temp_second_sse);

      // print the number of bodies
      BCL_MessageStd( "there are " + util::Format()( bodies->GetSize()) + " bodies in bodies");
      BCL_MessageStd( "there are " + util::Format()( bodies->GetSize()) + " bodies in two_bodies");

      BCL_ExampleMustOpenInputFile( read, pdb_filename);

      // create restraint::Body out of "bodies" and way to check the occupancies
      util::ShPtr< restraint::Body> body_restraint
      (
        restraint::HandlerBody
        (
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin()
          )
        ).CreateRestraintsBody( read).FirstElement()
      );
      io::File::CloseClearFStream( read);

      // create restraint::Body out of "bodies" and way to check the occupancies
      util::ShPtr< restraint::Body> body_restraint_small
      (
        new restraint::Body
        (
          // initialize with "bodies"
          coord_two_bodies,
          // initialize with ShPtr to a restraint::ContainsBodyOrigin to determine occupancies
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin()
          )
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::BodyConnectivityDensity def_constr;

      // constructor taking a ShPtr< restraint::Body>, a util::ShPtr< coord::DensityMap> and a length reduction
      score::BodyConnectivityDensity connectivity_score_object
      (
        body_restraint,
        util::ShPtr< density::Map>( density_map.Clone())
      );
      // constructor taking a ShPtr< restraint::Body>, a util::ShPtr< coord::DensityMap> and a length reduction
      // for small object
      score::BodyConnectivityDensity connectivity_score_object_small
      (
        body_restraint_small,
        util::ShPtr< density::Map>( density_map.Clone())
      );

      // check whether restraint body was read in correctly
      BCL_ExampleIndirectAssert( connectivity_score_object.GetBodyRestraint(), body_restraint, "I/O");

      // check that map of connectivity objects and scores was read in correctly (has the right size at least)
      BCL_ExampleCheck( connectivity_score_object.GetScores().GetSize(), 144);

    ////////////////
    // operations //
    ////////////////

      // check GetScores function for default constructed object
      BCL_ExampleIndirectCheck
      (
        def_constr.GetScores().IsEmpty(),
        true,
        "Default constructor should not create connectivities"
      );

      // get the bodies directly from the restraints (necessary as the fragment sizes have been changed in the restraint)
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > bodies_from_restraint( body_restraint->GetBody());

      // get the corresponding SSEGeometryInterfaces
      const util::ShPtrVector< assemble::SSEGeometryInterface> vector_of_bodies_from_restraint
      (
        bodies_from_restraint->Begin(),
        bodies_from_restraint->End()
      );

      // get the corresponding SSEGeometryInterfaces
      const util::ShPtrList< assemble::SSEGeometryInterface> list_of_bodies_from_restraint
      (
        bodies_from_restraint->Begin(),
        bodies_from_restraint->End()
      );

      size_t helix_counter_a( 0);
      size_t helix_counter_b( 0);
      storage::List< density::Connectivity> connectivities_individual;

      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator iter_a( vector_of_bodies_from_restraint.Begin()),
          iter_a_end( vector_of_bodies_from_restraint.End());
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
          connectivities_individual =
            density::Connectivity::DetermineConnectivities
            (
              density_map, *iter_a, *iter_b
            );

          // iterate through the 4 resulting connectivities
          for
          (
            storage::List< density::Connectivity>::const_iterator
              connectivity_iter( connectivities_individual.Begin()),
              connectivity_iter_end( connectivities_individual.End());
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

      // create list of connectivities from density map and list of sses directly
      storage::List< density::Connectivity> connectivities_1IE9
      (
        density::Connectivity::DetermineConnectivities( density_map, list_of_bodies_from_restraint)
      );

      // check that map of connectivities was constructed correctly
      BCL_ExampleCheck( connectivities_1IE9.GetSize(), 144);

      // check GetScores function for connectivity_score_object by finding the connectivities of connectivities_1IE9 in
      // the connectivity_score_object map
      const storage::Map< density::Connectivity, double, density::Connectivity::LessThan>
        score_map( connectivity_score_object.GetScores());

      // check that map of connectivity objects and scores was retrieved correctly
      BCL_ExampleCheck( score_map.GetSize(), 144);

      bool all_connectivities_were_in_map( true);

      // iterate over all the connectivities in connectivities_1IE9
      for
      (
        storage::List< density::Connectivity>::const_iterator
          list_itr( connectivities_1IE9.Begin()), list_itr_end( connectivities_1IE9.End());
        list_itr != list_itr_end; ++list_itr
      )
      {
        const storage::Map< density::Connectivity, double, density::Connectivity::LessThan>::const_iterator
        connectivity_itr( score_map.Find( *list_itr));

        BCL_MessageDbg
        (
          "score found in map: " + util::Format()( connectivity_itr->second) +
            " , distance: " + util::Format()( connectivity_itr->first.GetDistance()) +
            " , connectivity: " + util::Format()( connectivity_itr->first.GetConnectivity())
        );

        all_connectivities_were_in_map = all_connectivities_were_in_map && connectivity_itr != score_map.End();
      }

      BCL_ExampleIndirectCheck
      (
        all_connectivities_were_in_map, true,
        "connectivity from list should be in score_map"
      );
      if( !all_connectivities_were_in_map)
      {
        BCL_MessageStd( "The actual map: " + util::Format()( score_map));
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "Test write and read");

      // test the write and read
      WriteBCLObject( connectivity_score_object_small);
      score::BodyConnectivityDensity connectivity_score_object_to_be_read;
      ReadBCLObject( connectivity_score_object_to_be_read);

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBodyAssignment

  const ExampleClass::EnumType ExampleScoreBodyConnectivityDensity::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreBodyConnectivityDensity())
  );

} // namespace bcl
