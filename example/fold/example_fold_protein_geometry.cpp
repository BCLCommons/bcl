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
#include "fold/bcl_fold_protein_geometry.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_protein_geometry.cpp
  //! @brief this example tests the implementation of fold::XlinkParameters
  //!
  //! @author fischea
  //! @date Oct 12, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldProteinGeometry :
     public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to a new ExampleFoldProteinGeometry
    ExampleFoldProteinGeometry *Clone() const
    {
      return new ExampleFoldProteinGeometry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

      // prepare test data without loops
      const std::string test_pdb( AddExampleInputPathToFilename( e_Biology, "2AP3_no_loops.pdb"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, test_pdb);
      const pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( pdb));
      const biol::AABase &anchor_1( *model.GetAminoAcids()( 50));
      const biol::AABase &anchor_2( *model.GetAminoAcids()( 59));

      // prepare test data with loops
      const std::string test_pdb_loops( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      BCL_ExampleMustOpenInputFile( read, test_pdb_loops);
      const pdb::Handler pdb_loops( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model_loops( pdb::Factory().ProteinModelFromPDB( pdb_loops));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::ProteinGeometry default_geom;

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( default_geom.GetClassIdentifier(), ( GetStaticClassName< fold::ProteinGeometry>()));

    ////////////////
    // operations //
    ////////////////

      // correct solution for the local coordinate systems
      linal::Matrix3x3< double> correct_base_1( util::GetUndefinedDouble());
      correct_base_1( 0, 0) = 0.751;
      correct_base_1( 0, 1) = -0.446;
      correct_base_1( 0, 2) = -0.487;
      correct_base_1( 1, 0) = -0.078;
      correct_base_1( 1, 1) = 0.672;
      correct_base_1( 1, 2) = -0.737;
      correct_base_1( 2, 0) = 0.655;
      correct_base_1( 2, 1) = 0.592;
      correct_base_1( 2, 2) = 0.470;
      linal::Matrix3x3< double> correct_base_2( util::GetUndefinedDouble());
      correct_base_2( 0, 0) = 0.193;
      correct_base_2( 0, 1) = -0.855;
      correct_base_2( 0, 2) = -0.481;
      correct_base_2( 1, 0) = -0.906;
      correct_base_2( 1, 1) = 0.033;
      correct_base_2( 1, 2) = -0.422;
      correct_base_2( 2, 0) = 0.377;
      correct_base_2( 2, 1) = 0.518;
      correct_base_2( 2, 2) = -0.768;

      // compute the local coordinate systems
      linal::Matrix3x3< double> base_1( fold::ProteinGeometry::GetLocalCoordinateSystem( anchor_1));
      linal::Matrix3x3< double> base_2( fold::ProteinGeometry::GetLocalCoordinateSystem( anchor_2));

      // check the correctness of the computation for base_1
      bool correct_1( true);
      for( size_t x( 0); x < 3; ++x)
      {
        for( size_t y( 0); y < 3; ++y)
        {
          if( !math::EqualWithinTolerance( correct_base_1( x, y), base_1( x, y), 0.001, 0.001))
          {
            correct_1 = false;
            break;
          }
        }
      }
      BCL_Example_Check( correct_1, "base_1 computed incorrectly.")

      // check the correctness of the computation for base_2
      bool correct_2( true);
      for( size_t x( 0); x < 3; ++x)
      {
        for( size_t y( 0); y < 3; ++y)
        {
          if( !math::EqualWithinTolerance( correct_base_2( x, y), base_2( x, y), 0.001, 0.001))
          {
            correct_2 = false;
            break;
          }
        }
      }
      BCL_Example_Check( correct_2, "base_2 computed incorrectly.");

      // test if templates can be combined
      const biol::AASequence &loop_1( *model_loops.GetSSEs( biol::GetSSTypes().COIL)( 5));
      const util::SiPtrVector< const biol::AABase> aas_1( loop_1.GetMembers());

      // combine sequences
      const util::ShPtr< biol::AASequence> sp_combined_1
      (
        fold::ProteinGeometry::CombineSequences( loop_1, loop_1, 2, 3)
      );
      const util::ShPtr< biol::AASequence> sp_combined_2
      (
        fold::ProteinGeometry::CombineSequences( loop_1, loop_1, 3, 4)
      );

      // check if combined sequence 1 is the same as the original sequence
      BCL_ExampleCheckWithinTolerance
      (
        loop_1.GetFirstAA()->GetCenter(), sp_combined_1->GetFirstAA()->GetCenter(), 0.01
      );
      BCL_ExampleCheck( loop_1.GetSize(), sp_combined_1->GetSize());
      for( size_t i( 122); i <= 125; ++i)
      {
        const storage::VectorND< 2, double> dihedral_orig( loop_1.CalculatePhiPsi( i));
        const storage::VectorND< 2, double> dihedral_comb( sp_combined_1->CalculatePhiPsi( i));
        BCL_ExampleCheckWithinTolerance( dihedral_orig( 0), dihedral_comb( 0), 0.01);
        BCL_ExampleCheckWithinTolerance( dihedral_orig( 1), dihedral_comb( 1), 0.01);
      }

      // check if combined sequence 2 is the same as the original sequence
      BCL_ExampleCheckWithinTolerance
      (
        loop_1.GetFirstAA()->GetCenter(), sp_combined_2->GetFirstAA()->GetCenter(), 0.01
      );
      BCL_ExampleCheck( loop_1.GetSize(), sp_combined_2->GetSize());
      for( size_t i( 122); i <= 125; ++i)
      {
        const storage::VectorND< 2, double> dihedral_orig( loop_1.CalculatePhiPsi( i));
        const storage::VectorND< 2, double> dihedral_comb( sp_combined_2->CalculatePhiPsi( i));
        BCL_ExampleCheckWithinTolerance( dihedral_orig( 0), dihedral_comb( 0), 0.01);
        BCL_ExampleCheckWithinTolerance( dihedral_orig( 1), dihedral_comb( 1), 0.01);
      }

      return 0;
    }

  }; // class ExampleFoldProteinGeometry

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldProteinGeometry::s_Instance
  (
     GetExamples().AddEnum( ExampleFoldProteinGeometry())
  );

} // namespace bcl
