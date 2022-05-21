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
#include "assemble/bcl_assemble_domain.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "biol/bcl_biol_atom.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_domain.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleDomain :
    public ExampleInterface
  {
  public:

    ExampleAssembleDomain *Clone() const
    {
      return new ExampleAssembleDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
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
      assemble::Domain def_construct;
      BCL_Example_Check
      (
        def_construct.GetData().IsEmpty(), "Default constructor failed"
      );

      // test constructor from a Set of ShPtr to SSEs
      BCL_MessageStd( "testing constructor from a Set of ShPtr to SSEs");
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> sse_set;
      util::ShPtrVector< assemble::SSE> sse_vector;
      const util::SiPtrVector< const assemble::SSE> sses_from_model( protein_model.GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses_from_model.Begin()),
          sse_itr_end( sses_from_model.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        sse_set.Insert( util::ShPtr< assemble::SSE>( ( *sse_itr)->Clone()));
        sse_vector.PushBack( util::ShPtr< assemble::SSE>( ( *sse_itr)->Clone()));
      }
      assemble::Domain set_construct( sse_set);
      BCL_Example_Check
      (
        set_construct.GetData().GetSize() == sses_from_model.GetSize(),
        "Constructor from a Set of ShPtr to SSEs gave: " + util::Format()( set_construct.GetData()) + " instead of: " +
        util::Format()( sses_from_model)
      );

      // test constructor from ShPtrVector of SSEs
      BCL_MessageStd( "testing constructor from ShPtrVector of SSEs");
      assemble::Domain vector_construct( sse_vector);
      BCL_Example_Check
      (
        vector_construct.GetData().GetSize() == sses_from_model.GetSize(),
        "Constructor from ShPtrVector of SSEs gave: " + util::Format()( vector_construct.GetData()) + " instead of: " +
        util::Format()( sses_from_model)
      );

      // test copy constructor
      BCL_MessageStd( "testing copy constructor");
      assemble::Domain copy_construct( set_construct);
      BCL_Example_Check
      (
        copy_construct.GetData().InternalData() == set_construct.GetData().InternalData(),
        "Copy construct gave: " + util::Format()( copy_construct.GetData()) + " instead of: " +
        util::Format()( set_construct.GetData())
      );

      // test Clone function
      BCL_MessageStd( "testing Clone function");
      util::ShPtr< assemble::Domain> clone_construct( set_construct.Clone());
      BCL_Example_Check
      (
        clone_construct->GetData().InternalData() == set_construct.GetData().InternalData(),
        "Clone function gave: " + util::Format()( clone_construct->GetData()) + " instead of: " +
        util::Format()( set_construct.GetData())
      );

      // test HardCopy function
      BCL_MessageStd( "testing HardCopy function");
      util::ShPtr< assemble::Domain> hard_copy_construct( set_construct.HardCopy());
      BCL_Example_Check
      (
        hard_copy_construct->GetData().GetSize() == set_construct.GetData().GetSize(),
        "HardCopy function gave: " + util::Format()( hard_copy_construct->GetData()) + " instead of: " +
        util::Format()( set_construct.GetData())
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::Domain");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::Domain>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::Domain>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::Domain>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetNumberSSEs
      const size_t known_number_of_sses( 5);
      BCL_Example_Check
      (
        set_construct.GetNumberSSEs() == known_number_of_sses,
        "GetNumberSSEs gives " + util::Format()( set_construct.GetNumberSSEs()) + " but should give " +
        util::Format()( known_number_of_sses)
      );

      // check GetNumberSSE
      const size_t known_number_of_strands( 4);
      BCL_Example_Check
      (
        set_construct.GetNumberSSE( biol::GetSSTypes().STRAND) == known_number_of_strands,
        "GetNumberSSE gives " + util::Format()( set_construct.GetNumberSSE( biol::GetSSTypes().STRAND)) +
        " but should give " + util::Format()( known_number_of_strands)
      );

      // check GetSSEs
      BCL_Example_Check
      (
        set_construct.GetSSEs( biol::GetSSTypes().STRAND).GetSize() == known_number_of_strands,
        "GetSSEs has size " + util::Format()( set_construct.GetSSEs( biol::GetSSTypes().STRAND).GetSize()) +
        " but should have " + util::Format()( known_number_of_strands)
      );

      // GetData checked w/ constructors

      // check SetData
      def_construct.SetData( sse_set);
      BCL_Example_Check
      (
        def_construct.GetData().GetSize() == sses_from_model.GetSize(),
        "SetData has size: " + util::Format()( def_construct.GetData()) + " instead of: " +
        util::Format()( sses_from_model)
      );

      // check GetAtoms
      const size_t known_number_of_atoms( 210);
      BCL_Example_Check
      (
        set_construct.GetAtoms().GetSize() == known_number_of_atoms,
        "GetAtoms has size " + util::Format()( set_construct.GetAtoms().GetSize()) +
        " but should have " + util::Format()( known_number_of_atoms)
      );

      // check GetAtomCoordinates
      BCL_Example_Check
      (
        set_construct.GetAtomCoordinates().GetSize() == known_number_of_atoms,
        "GetAtomCoordinates has size " + util::Format()( set_construct.GetAtomCoordinates().GetSize()) +
        " but should have " + util::Format()( known_number_of_atoms)
      );

      // check GetAtomCoordinates of specific atom type
      const size_t known_number_of_residues( 42);
      BCL_Example_Check
      (
        set_construct.GetAtomCoordinates( biol::GetAtomTypes().CA).GetSize() == known_number_of_residues,
        "GetAtomCoordinates has size " +
        util::Format()( set_construct.GetAtomCoordinates( biol::GetAtomTypes().CA).GetSize()) + " but should have " +
        util::Format()( known_number_of_residues)
      );

      // check GetCenter
      const linal::Vector3D known_center( 31.2505, 31.3569, 13.7259);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( set_construct.GetCenter(), known_center),
        "GetCenter gives " + util::Format()( set_construct.GetCenter()) +
        " but should give " + util::Format()( known_center)
      );

      // check CreateSequenceFromSSEs
      const size_t total_known_residues( 72);
      BCL_Example_Check
      (
        set_construct.CreateSequenceFromSSEs()->GetSize() == total_known_residues,
        "CreateSequenceFromSSEs has size " + util::Format()( set_construct.CreateSequenceFromSSEs()->GetSize()) +
        " but should have " + util::Format()( total_known_residues)
      );

      // check GetAminoAcids 
      BCL_Example_Check
      (
        set_construct.GetAminoAcids().GetSize() == known_number_of_residues,
        "GetAminoAcids has size " + util::Format()( set_construct.GetAminoAcids().GetSize()) +
        " but should have " + util::Format()( known_number_of_residues)
      );

    ////////////////
    // operations //
    ////////////////

      // check Translate
      const linal::Vector3D translate( set_construct.GetCenter());
      const linal::Vector3D origin( 0.0, 0.0, 0.0);
      set_construct.Translate( -translate);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( set_construct.GetCenter(), origin, 0.01, 0.01),
        "Translate moves to " + util::Format()( set_construct.GetCenter()) +
        " but should be " + util::Format()( origin)
      );

      // check Transform
      set_construct.Translate( translate);
      math::TransformationMatrix3D transform( protein_model.GetOrientation());
      set_construct.Transform( transform.Invert());
      BCL_Example_Check
      (
        math::EqualWithinTolerance( set_construct.GetCenter(), origin, 0.1, 0.1),
        "Transform moves to " + util::Format()( set_construct.GetCenter()) +
        " but should be " + util::Format()( origin)
      );

      // check Rotate
      set_construct.Rotate( transform.GetRotation());
      BCL_Example_Check
      (
        math::EqualWithinTolerance( set_construct.GetCenter(), origin, 0.1, 0.1),
        "Transform moves to " + util::Format()( set_construct.GetCenter()) +
        " but should be " + util::Format()( origin)
      );

      // check Insert
      util::ShPtr< assemble::SSE> first_sse( sse_set.Begin()->HardCopy());
      sse_set.RemoveElement( sse_set.Begin());
      assemble::Domain truncated_domain( sse_set);
      BCL_Example_Check
      (
        truncated_domain.Insert( first_sse) && truncated_domain.GetNumberSSEs() == known_number_of_sses,
        "After Insert there are " + util::Format()( truncated_domain.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses)
      );

      // check Replace
      BCL_Example_Check
      (
        set_construct.Replace( first_sse) && set_construct.GetNumberSSEs() == known_number_of_sses,
        "After Replace there are " + util::Format()( set_construct.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses)
      );

      // check ReplaceWithOverlapping
      BCL_Example_Check
      (
        set_construct.ReplaceWithOverlapping( first_sse) && set_construct.GetNumberSSEs() == known_number_of_sses,
        "After ReplaceWithOverlapping there are " + util::Format()( set_construct.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses)
      );

      // check Remove
      BCL_Example_Check
      (
        copy_construct.Remove( *first_sse) && copy_construct.GetNumberSSEs() == known_number_of_sses - 1,
        "After Remove there are " + util::Format()( copy_construct.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses - 1)
      );

      // check SetToIdealConformation
      vector_construct.SetToIdealConformation();
      BCL_Example_Check
      (
        math::EqualWithinTolerance( vector_construct.GetCenter(), known_center, 0.1, 0.1),
        "SetToIdealConformation gives a center of " + util::Format()( vector_construct.GetCenter()) +
        " but should give " + util::Format()( known_center)
      );

      // check ChopSSElements
      assemble::Domain chopped_domain( set_construct);
      const size_t known_number_of_chopped_sses( 9);
      chopped_domain.ChopSSEs( storage::VectorND< 3, size_t>( 5, 3, 3));
      BCL_Example_Check
      (
        chopped_domain.GetNumberSSEs() == known_number_of_chopped_sses,
        "After ChopSSElements there are " + util::Format()( chopped_domain.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_chopped_sses)
      );

      // check DoesContain
      BCL_Example_Check
      (
        set_construct.DoesContain( *first_sse),
        "DoesContain failed to find an SSE in the domain"
      );

      // check DoesContainOverlapping
      BCL_Example_Check
      (
        set_construct.DoesContainOverlapping( *first_sse),
        "DoesContainOverlapping failed to find an SSE in the domain"
      );

      // check Join
      // break up the second sse
      const assemble::SSE &second_sse( *sse_vector( 1));
      util::ShPtr< assemble::SSE> sp_fragment_a
      (
        new assemble::SSE( second_sse.SubSequence( 0, 2), second_sse.GetType())
      );
      util::ShPtr< assemble::SSE> sp_fragment_b
      (
        new assemble::SSE( second_sse.SubSequence( 2, 2), second_sse.GetType())
      );

      assemble::Domain join_construct( vector_construct);
      // replace the second sse with two adjoining fragments
      join_construct.ReplaceWithOverlapping( sp_fragment_a);
      join_construct.ReplaceWithOverlapping( sp_fragment_b);
      join_construct.Join( biol::GetSSTypes().STRAND);
      BCL_Example_Check
      (
        join_construct.GetNumberSSEs() == known_number_of_sses,
        "After Join there are " + util::Format()( join_construct.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses)
      );

      // check GetSSEsWithShortLoops
      const size_t known_short_loop_sses( 1);
      BCL_Example_Check
      (
        set_construct.GetSSEsWithShortLoops( *first_sse, 7).GetSize() == known_short_loop_sses,
        "GetSSEsWithShortLoops gives " + util::Format()( set_construct.GetSSEsWithShortLoops( *first_sse, 7).GetSize())
        + " SSEs, but there should be " + util::Format()( known_short_loop_sses)
      );

      // check GetSSEsWithShortLoops
      const size_t total_short_loop_sses( 6);
      BCL_Example_Check
      (
        set_construct.GetSSEsWithShortLoops( 7).GetSize() == total_short_loop_sses,
        "GetSSEsWithShortLoops gives " + util::Format()( set_construct.GetSSEsWithShortLoops( 7).GetSize()) +
        " SSEs, but there should be " + util::Format()( total_short_loop_sses)
      );

      // check GetNeighborSSEs
      BCL_Example_Check
      (
        set_construct.GetNeighborSSEs( *first_sse).Second()->IsDefined(),
        "GetNeighborSSEs gives an undefined entry"
      );

      // check GetOverlappingSSEs
      BCL_Example_Check
      (
        set_construct.GetOverlappingSSEs( *first_sse).GetSize() == 1,
        "GetOverlappingSSEs did not return the overlapping SSE"
      );

      // check FilterByMinSSESizes
      assemble::Domain filtered_domain( set_construct);
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 7;
      filtered_domain.FilterByMinSSESizes( ssetype_min_size);
      BCL_Example_Check
      (
        filtered_domain.GetNumberSSEs() == known_number_of_sses - 1,
        "FilterByMinSSESizes gives " + util::Format()( filtered_domain.GetNumberSSEs()) +
        " SSEs, but there should be " + util::Format()( known_number_of_sses - 1)
      );

    ///////////////
    // operators //
    ///////////////

      // check = operator
      def_construct = set_construct;
      BCL_Example_Check
      (
        def_construct.GetData().InternalData() == set_construct.GetData().InternalData(),
        "= operator gives " + util::Format()( def_construct.GetData()) +
        " instead of " + util::Format()( set_construct.GetData())
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( set_construct);

      // read the object back in
      assemble::Domain domain_read;
      ReadBCLObject( domain_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_Example_Check
      (
        set_construct.GetNumberSSEs() == domain_read.GetNumberSSEs(),
        "the written and read Domain classes differ from each other \n" +
        util::Format()( set_construct) + "\nvs\n" + util::Format()( domain_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleDomain

  const ExampleClass::EnumType ExampleAssembleDomain::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleDomain())
  );

} // namespace bcl

