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
#include "fold/bcl_fold_mutate_aa_rotate.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_aa_rotate.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateAARotate :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateAARotate *Clone() const
    {
      return new ExampleFoldMutateAARotate( *this);
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

      // create test residue
      biol::AABackBone residue( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA)));

      // set to ideal conformation
      residue.SetToIdealConformation( biol::GetSSTypes().HELIX, math::TransformationMatrix3D());

      {
        pdb::Handler newpdb;
        size_t atom_number( 1);
        newpdb.AppendLines( pdb::Factory::WriteResiduesToLines( residue, 'A', atom_number));
        const std::string out_pdb_name( AddExampleOutputPathToFilename( residue, "mutate_aa_rotate_start.pdb"));
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        newpdb.WriteLines( write);
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      fold::MutateAARotate default_constr;

      // test constructor taking parameters

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

        {
          const biol::Atom &ca_atom( residue.GetCA());
          const biol::Atom &c_atom( residue.GetAtom( biol::GetAtomTypes().C));
          const coord::LineSegment3D rotation_axis( ca_atom.GetCoordinates(), c_atom.GetCoordinates());
          const double rotation( math::g_Pi / 4.0);
          fold::MutateAARotate mutate_rotate
          (
            rotation_axis,
            ca_atom.GetCoordinates(), //< rotation origin
            rotation
          );
          util::ShPtr< biol::AABase> rotated_aa( mutate_rotate( residue).GetArgument());
          pdb::Handler newpdb;
          size_t atom_number( 1);
          newpdb.AppendLines( pdb::Factory::WriteResiduesToLines( *rotated_aa, 'A', atom_number));
          const std::string out_pdb_name( AddExampleOutputPathToFilename( *rotated_aa, "mutate_aa_rotate_45.pdb"));
          BCL_MessageDbg( "printing to file " + out_pdb_name);
          const std::string correct_pdb_name( out_pdb_name + ".correct");
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
          newpdb.WriteLines( write);
          io::File::CloseClearFStream( write);
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( out_pdb_name, correct_pdb_name, 1.0e-6), true
          );
        }

        {
          const biol::Atom &ca_atom( residue.GetCA());
          const biol::Atom &n_atom( residue.GetAtom( biol::GetAtomTypes().N));
          const coord::LineSegment3D rotation_axis( n_atom.GetCoordinates(), ca_atom.GetCoordinates());
          const double rotation( math::g_Pi / 2.0);
          fold::MutateAARotate mutate_rotate
          (
            rotation_axis,
            residue.GetAtom( biol::GetAtomTypes().C).GetCoordinates(), //< rotation origin
            rotation
          );
          util::ShPtr< biol::AABase> rotated_aa( mutate_rotate( residue).GetArgument());
          pdb::Handler newpdb;
          size_t atom_number( 1);
          newpdb.AppendLines( pdb::Factory::WriteResiduesToLines( *rotated_aa, 'A', atom_number));
          const std::string out_pdb_name( AddExampleOutputPathToFilename( *rotated_aa, "mutate_aa_rotate_90.pdb"));
          BCL_MessageDbg( "printing to file " + out_pdb_name);
          const std::string correct_pdb_name( out_pdb_name + ".correct");
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
          newpdb.WriteLines( write);
          io::File::CloseClearFStream( write);
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( out_pdb_name, correct_pdb_name, 1.0e-6), true
          );
        }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateAARotate

  const ExampleClass::EnumType ExampleFoldMutateAARotate::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateAARotate())
  );

} // namespace bcl
