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
#include "find/bcl_find_pick_body_random.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_pick_body_random.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindPickBodyRandom :
    public ExampleInterface
  {
  public:

    ExampleFindPickBodyRandom *Clone() const
    { return new ExampleFindPickBodyRandom( *this);}

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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // get the strand secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> all_sses
      (
        protein_model.GetChains()( 0)->GetSSEs()
      );

      // create a ShPtrVector of bodies which will be created from "all_sses"
      util::ShPtrVector< assemble::SSEGeometryInterface> bodies;

      // add the coord::Bodies of "all_sses" to "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> temp_sse( ( *sse_itr)->Clone());
        bodies.PushBack( temp_sse);
      }

      // create assemble::PickBodyRandom
      find::PickBodyRandom body_pick;

      // use "body_pick" to get a random body out of "bodies"
      util::ShPtr< coord::GeometryInterface> random_body( body_pick.Pick( bodies));

      BCL_ExampleAssert( body_pick.Pick( bodies).IsDefined(), true);

      BCL_MessageStd
      (
        "the body that was randomly selected is : \n"
        + util::Format()( *random_body)
      );

      // initialize a storage::Set to hold the bodies that are picked
      storage::Set< util::ShPtr< coord::GeometryInterface> > picked_bodies;

      for( size_t counter( 0); counter < 1000 && picked_bodies.GetSize() < 7; ++counter)
      {
        picked_bodies.Insert( body_pick.Pick( bodies));
      }

      BCL_ExampleIndirectCheck
      (
        picked_bodies.GetSize(),
        7,
        "random picking should pick all seven bodies at least once within 1000 attempts"
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFindPickBodyRandom

  const ExampleClass::EnumType ExampleFindPickBodyRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleFindPickBodyRandom())
  );

} // namespace bcl
