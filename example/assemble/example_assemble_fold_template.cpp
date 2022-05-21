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
#include "assemble/bcl_assemble_fold_template.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_fold_template.cpp
  //! @details This example checks all functions in the fold template class.
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by linders on Nov 6, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleFoldTemplate :
    public ExampleInterface
  {
  public:

    ExampleAssembleFoldTemplate *Clone() const
    {
      return new ExampleAssembleFoldTemplate( *this);
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

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      const std::string default_pdb_id( "DFLT");
      assemble::FoldTemplate def_construct;
      BCL_Example_Check
      (
        def_construct.GetGeometries().GetSize() == 0 &&
        def_construct.GetPDBID() == default_pdb_id,
        "Default constructor failed"
      );

      // test constructor from geometry vectors
      BCL_MessageStd( "testing vector construct");
      util::ShPtrVector< assemble::SSEGeometryPhiPsi> geometries;
      const util::SiPtrVector< const assemble::SSE> sse_vector( protein_model.GetSSEs());
      const util::ShPtr< assemble::CollectorTopologyInterface> sp_collector( new assemble::CollectorTopologyCombined());
      // iterate through sse_vector in PROTEIN_MODEL
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // add geometry to data vector
        geometries.PushBack( util::ShPtr< assemble::SSEGeometryPhiPsi>( new assemble::SSEGeometryPhiPsi( **sse_itr)));
      }
      assemble::FoldTemplate vector_construct( geometries, sp_collector);

      BCL_Example_Check
      (
        vector_construct.GetGeometries().GetSize() == geometries.GetSize() &&
        vector_construct.GetPDBID() == default_pdb_id,
        "Geometry vector constructor failed"
      );

      // test constructor with protein model (get functions checked below)
      BCL_MessageStd( "testing protein model construct and PDB ID");
      assemble::FoldTemplate model_construct( protein_model, sp_collector, "1UBI");
      BCL_MessageStd
      (
        "fold template contains the following geometries:\n" + model_construct.GetGeometryIndentifications()
      );

      // test copy constructor
      BCL_MessageStd( "testing copy constructor");
      assemble::FoldTemplate copy_construct( model_construct);
      BCL_Example_Check
      (
        EqualBodies( model_construct.GetGeometries(), copy_construct.GetGeometries()) &&
        model_construct.GetPDBID() == copy_construct.GetPDBID(),
        "copy should be " + util::Format()( model_construct) + " but instead is " + util::Format()( copy_construct)
      );

      // test clone constructor
      BCL_MessageStd( "testing clone constructor");
      util::ShPtr< assemble::FoldTemplate> clone_construct( model_construct.Clone());
      BCL_Example_Check
      (
        EqualBodies( model_construct.GetGeometries(), clone_construct->GetGeometries()) &&
        model_construct.GetPDBID() == clone_construct->GetPDBID(),
        "clone should be " + util::Format()( model_construct) + " but instead is " + util::Format()( *clone_construct)
      );

    /////////////////
    // data access //
    /////////////////

      const std::string correct_static_class_name( "bcl::assemble::FoldTemplate");

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::FoldTemplate>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetGeometryIndentifications (check size and first element)
      const std::string identifications( model_construct.GetGeometryIndentifications());
      // extract first string
      const std::string first_identification( util::SplitString( identifications, "\n").FirstElement());
      const std::string correct_identification( "STRAND A   41 GLN <==>   44 ILE");
      BCL_Example_Check
      (
        first_identification == correct_identification,
        "GetGeometryIndentifications is not returning the proper first identification "
      );

      // check GetHelicalGeometries
      const size_t number_of_helices( 1);
      BCL_Example_Check
      (
        model_construct.GetHelicalGeometries().GetSize() == number_of_helices,
        "GetHelicalGeometries size is " + util::Format()( model_construct.GetHelicalGeometries().GetSize()) +
        " instead of " + util::Format()( number_of_helices)
      );

      // check GetStrandGeometries
      const size_t number_of_strands( 4);
      BCL_Example_Check
      (
        model_construct.GetStrandGeometries().GetSize() == number_of_strands,
        "GetStrandGeometries size is " + util::Format()( model_construct.GetStrandGeometries().GetSize()) +
        " instead of " + util::Format()( number_of_strands)
      );

      // check GetGeometries
      BCL_Example_Check
      (
        model_construct.GetGeometries().GetSize() == geometries.GetSize(),
        "GetGeometries doesn't return the same number of geometries as the original geometries"
      );

      // check GetPDBID
      const std::string known_pdb_id( "1UBI");
      BCL_Example_Check
      (
        model_construct.GetPDBID() == known_pdb_id,
        "GetPDBID returns " + model_construct.GetPDBID() + " instead of " + known_pdb_id
      );

      // check IsTopologyInitialized
      BCL_Example_Check
      (
        model_construct.IsTopologyInitialized() == true,
        "IsTopologyInitialized returns " + util::Format()( model_construct.IsTopologyInitialized()) + " instead of true"
      );

      // check GetTopology
      BCL_Example_Check
      (
        model_construct.GetTopology().GetGraph().GetNumberEdges() == 14,
        "GetTopology returns  graph with " +
          util::Format()( model_construct.GetTopology().GetGraph().GetNumberEdges()) + " edges instead of 14"
      );

      // check SetTopology
      vector_construct.SetTopology( sp_collector->CalculateTopology( sse_vector));
      BCL_Example_Check
      (
        vector_construct.GetTopology().GetGraph().GetNumberVertices() > 0,
        "SetTopology doesn't work properly!!"
      );

      // check GetCenter
      const linal::Vector3D known_center_of_mass( 31.2519, 31.4694, 13.4158);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( model_construct.GetCenter(), known_center_of_mass),
        "GetCenter gives " + util::Format()( model_construct.GetCenter()) +
        " instead of " + util::Format()( known_center_of_mass)
      );

      // check GetRadiusOfGyration
      const double known_rg( 8.35665);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( model_construct.GetRadiusOfGyration(), known_rg),
        "Radius of gyration gives " + util::Format()( model_construct.GetRadiusOfGyration()) +
        " instead of " + util::Format()( known_rg)
      );

      // check GetSubdomain
      assemble::FoldTemplate subdomain( model_construct.GetSubDomain( 1, 2));
      BCL_Example_Check
      (
        subdomain.GetHelicalGeometries().GetSize() == 1 && subdomain.GetStrandGeometries().GetSize() == 2,
        "GetSubdomain did not produce the right number of helices (1) and strands (2):\n" +
        util::Format()( subdomain.GetGeometryIndentifications())
      );

      // check second GetSubdomain
      assemble::FoldTemplate subdomain_b( model_construct.GetSubDomain( sse_vector));
      BCL_Example_Check
      (
        subdomain_b.GetHelicalGeometries().GetSize() == 1 && subdomain_b.GetStrandGeometries().GetSize() == 4,
        "GetSubdomain did not produce the right number of helices (1) and strands (4):\n" +
        util::Format()( subdomain_b.GetGeometryIndentifications())
      );

    ////////////////
    // operations //
    ////////////////

      // check Translate
      model_construct.Translate( -( model_construct.GetCenter()));
      // center should be at 0, 0, 0 now
      BCL_Example_Check
      (
        math::EqualWithinTolerance( model_construct.GetCenter(), linal::Vector3D( 0.0, 0.0, 0.0), 0.001, 0.000001),
        "Translate doesn't shift it into the origin, but to " + util::Format()( model_construct.GetCenter())
      );

      // check Transform
      math::TransformationMatrix3D helix_body( model_construct.GetHelicalGeometries().LastElement()->GetOrientation());
      const linal::Vector3D translate( 1.0, 1.0, 1.0);
      model_construct.Transform( helix_body.Invert());
      model_construct.Transform( math::TransformationMatrix3D( translate));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( model_construct.GetHelicalGeometries().LastElement()->GetCenter(), translate),
        "Transform gives " +
        util::Format()( model_construct.GetHelicalGeometries().LastElement()->GetCenter()) +
        " instead of " + util::Format()( translate)
      );

      // check Rotate
      // get rotation from helix geometry
      math::RotationMatrix3D rotation_matrix
      (
        model_construct.GetHelicalGeometries().LastElement()->GetOrientation().GetRotation()
      );
      // transform the entire model to the origin
      model_construct.Transform( math::Inverse( model_construct.GetOrientation()));
      // rotate model by a certain rotation matrix
      model_construct.Rotate( rotation_matrix);
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          model_construct.GetOrientation().GetRotation().GetAxis( coord::GetAxes().e_X),
          rotation_matrix.GetAxis( coord::GetAxes().e_X)
        ),
        "Rotate doesn't work properly!! "
      );

      // check FitSSEs
      assemble::Domain domain( model_construct.FitSSEs( sse_vector));
      BCL_Example_Check
      (
        domain.GetSSEs().GetSize() == 5,
        "FitSSEs failed to place SSEs correctly!! "
      );

    ///////////////
    // operators //
    ///////////////

      // test = operator
      def_construct = model_construct;
      BCL_Example_Check
      (
        EqualBodies( model_construct.GetGeometries(), def_construct.GetGeometries()) &&
        model_construct.GetPDBID() == def_construct.GetPDBID(),
        "assigned template should be " + util::Format()( model_construct) +
        " but instead is " + util::Format()( def_construct)
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( vector_construct);

      // read the object back in
      assemble::FoldTemplate fold_template_read;
      ReadBCLObject( fold_template_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_Example_Check
      (
        EqualBodies( vector_construct.GetGeometries(), fold_template_read.GetGeometries()) &&
        vector_construct.GetPDBID() == fold_template_read.GetPDBID(),
        "the written and read FoldTemplate classes differ from each other \n" +
        util::Format()( vector_construct) + "\nvs\n" + util::Format()( fold_template_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    //! @brief checks whether two body vectors are equal within the default tolerance
    //! @param BODIES_A ShPtrVector of SSEGeometries to be compared
    //! @param BODIES_B ShPtrVector of SSEGeometries to be compared
    //! @return bool whether two body vectors are equal within the default tolerance
    static bool EqualBodies
    (
      const util::ShPtrVector< assemble::SSEGeometryPhiPsi> &BODIES_A,
      const util::ShPtrVector< assemble::SSEGeometryPhiPsi> &BODIES_B
    )
    {
      // if the vectors are not the same length
      if( BODIES_A.GetSize() != BODIES_B.GetSize())
      {
        // return false
        return false;
      }

      // iterate through the vectors
      for
      (
        util::ShPtrVector< assemble::SSEGeometryPhiPsi>::const_iterator
          body_a_itr( BODIES_A.Begin()), body_a_itr_end( BODIES_A.End()),
          body_b_itr( BODIES_B.Begin()), body_b_itr_end( BODIES_B.End());
        body_a_itr != body_a_itr_end && body_b_itr != body_b_itr_end;
        ++body_a_itr, ++body_b_itr
      )
      {
        // check if the origins are different
        if
        (
          !math::SimilarWithinTolerance
          (
            ( *body_a_itr)->GetOrientation(), ( *body_b_itr)->GetOrientation(), 0.01, 0.01
          )
        )
        {
          // return false
          return false;
        }
      }

      // if this point is reached, the vectors are equal
      return true;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleFoldTemplate

  const ExampleClass::EnumType ExampleAssembleFoldTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleFoldTemplate())
  );

} // namespace bcl
