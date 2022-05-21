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
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_handler_atom_distance_assigned.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintHandlerAtomDistanceAssigned :
    public ExampleInterface
  {
  public:

    ExampleRestraintHandlerAtomDistanceAssigned *Clone() const
    { return new ExampleRestraintHandlerAtomDistanceAssigned( *this);}

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
      // test default constructor
      BCL_MessageStd( "test default constructor");
      // create HandlerAtomDistanceAssigned "handler"
      restraint::HandlerAtomDistanceAssigned handler;

      // create string "restraint_filename" which has path for example restraint file
      const std::string restraint_filename( AddExampleInputPathToFilename( e_Biology, "atom_distance_restraint_file.txt"));

      // create stream to restraint file
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_filename);

      // test CreateRestraints function
      BCL_MessageStd( "test CreateRestraints function");
      // create ShPtrVector "restraint"
      util::ShPtrVector< restraint::AtomDistance> restraints( handler.ReadRestraints( read));
      BCL_MessageStd
      (
        "The restraints in " + restraint_filename + " : \n" + util::Format()( restraints)
      );
      io::File::CloseClearFStream( read);

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // make sure created restraints are correct
      BCL_MessageStd( "make sure created restraints are correct");
      size_t correct_number_of_restraints( 2);
      BCL_Example_Check
      (
        restraints.GetSize() == correct_number_of_restraints,
        "number of restraint should be " + util::Format()( correct_number_of_restraints)
        + " but is " + util::Format()( restraints.GetSize())
      );

      // test WriteRestraints
      BCL_MessageStd( "test WriteRestraints");
      storage::Vector
      <
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
      > restraint_list;

      biol::Atom cb( biol::GetAtomTypes().CB);
      biol::Atom ca( biol::GetAtomTypes().CA);

      restraint_list.PushBack
      (
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
        (
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >
          (
            storage::Triplet< char, int, biol::Atom>( 'A', 23, cb),
            storage::Triplet< char, int, biol::Atom>( 'B', 26, ca)
          ),
          storage::VectorND< 3, double>( 30.0, 35.0, 25.0)
        )
      );

      restraint_list.PushBack
      (
        storage::Pair
        <
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
        >
        (
          storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >
          (
            storage::Triplet< char, int, biol::Atom>( 'C', 20, ca),
            storage::Triplet< char, int, biol::Atom>( 'D', 16, cb)
          ),
          storage::VectorND< 3, double>( 35.0, 55.0, 15.5)
        )
      );

      const std::string output_restraint_filename
      (
        AddExampleOutputPathToFilename( handler, "HandlerAtomDistanceAssigned_WriteRestraints.txt")
      );
      BCL_MessageStd( "test WriteRestraints output filename is " + output_restraint_filename);
      io::OFStream write_restraints;

      BCL_ExampleMustOpenOutputFile( write_restraints, output_restraint_filename);

      BCL_MessageStd( "test WriteRestraints writing restraints");
      handler.WriteRestraints( write_restraints, restraint_list);
      io::File::CloseClearFStream( write_restraints);

      std::string correct_filename( output_restraint_filename + ".correct");

      BCL_ExampleCheck( io::File::FilesMatch( output_restraint_filename, correct_filename), true);

      // read in restraints
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "contact.RR"));
      util::ShPtrVector< restraint::AtomDistance> contacts( handler.ReadRestraints( read));
      io::File::CloseClearFStream( read);

      // test GetDataConstruct
      const size_t nr_restraints( 4);
      BCL_ExampleIndirectCheck( contacts.GetSize(), nr_restraints, "Inference of file format");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDistanceAssigned

  const ExampleClass::EnumType ExampleRestraintHandlerAtomDistanceAssigned::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintHandlerAtomDistanceAssigned())
  );

} // namespace bcl
