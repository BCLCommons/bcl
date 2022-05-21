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
#include "fold/bcl_fold_mutate_protein_model_loop_domain_ccd.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_loop_domain_closure.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_loop_domain_ccd.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelLoopDomainCCD :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelLoopDomainCCD *Clone() const
    {
      return new ExampleFoldMutateProteinModelLoopDomainCCD( *this);
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstantPhiPsi
    //! TODO: add a brief comment
    //! TODO: add an general comment to this class
    //!
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class RandomNumberGenerator :
      public random::DistributionInterface
    {

    private:

    //////////
    // data //
    //////////

      bool m_RotatePhi;
      double m_RandomDouble;
      uint64_t m_RandomSizeT;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      RandomNumberGenerator( const bool ROTATE_PHI, const double RANDOM_DOUBLE, const size_t RANDOM_SIZE_T) :
        m_RotatePhi( ROTATE_PHI),
        m_RandomDouble( RANDOM_DOUBLE),
        m_RandomSizeT( RANDOM_SIZE_T)
      {
      }

      RandomNumberGenerator *Clone() const
      {
        return new RandomNumberGenerator( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name if this function is overwritten
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the seed
      //! @param SEED seed to be used
      //! @return the seed, that is was set too
      uint64_t SetSeed( const uint64_t SEED)
      {
        return 0;
      }

      //! @brief get the seed
      //! @return the seed, that was used to start the rng
      uint64_t GetSeed() const
      {
        return 0;
      }

      //! @brief default range for random double
      //! @return range, in which Double() generates a random number
      const math::Range< double> &GetDoubleRange() const
      {
        static const math::Range< double> s_default_range( math::RangeBorders::e_LeftClosed, 0.0, 1.0, math::RangeBorders::e_RightOpen);
        return s_default_range;
      }

      uint64_t Unsigned64BitInt() const
      {
        return m_RandomSizeT;
      }

      //! @brief random double
      //! @return random number in double range
      double Double() const
      {
        return m_RandomDouble;
      }

      //! @brief generate a random boolean
      //! @return random true or false
      bool Boolean() const
      {
        return m_RotatePhi;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class RandomNumberGenerator

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test rotation around psi of N-terminal anchor residue (the last residue in the N-terminal anchor sse)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_a.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( true, 0.1, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test rotation around phi of first residue in loop domain
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_b.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( true, 0.28, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test rotation around psi of first residue in loop domain
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_c.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( false, 0.28, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test rotation around phi of residue in interior (not first or last residue) of loop domain
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_d.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( true, 0.4, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test rotation around psi of last residue of loop domain
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_e.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( false, 0.8, out_pdb_name, -1.08035881, 0.0, "2LZM.pdb");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test with LoopDomainCToN rotation around the anchor sse anchor aa
      {
        RandomNumberGenerator rng( true, 0.8, 1);

        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // 21-23
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 21, 23), false));

        util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> > loop_domain_locator
        (
          new fold::LocatorLoopDomain( loop_segments, false)
        );

        const std::string input_pdb( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_c_to_n_21_23.pdb"));
        util::ShPtr< assemble::ProteinModel> model( Proteins::GetModel( input_pdb).Clone());

        util::ShPtr< assemble::ProteinModelData> pmd( model->GetProteinModelData()->Clone());
        util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > pmd_ldl_locators
        (
          new util::ShPtrList< fold::LocatorLoopDomain>
          (
            1, util::ShPtr< fold::LocatorLoopDomain>( new fold::LocatorLoopDomain( loop_segments, false))
          )
        );

        pmd->Insert( assemble::ProteinModelData::e_LoopDomainLocators, pmd_ldl_locators);
        model->SetProteinModelData( pmd);

        score::ProteinModelLoopDomainClosure score;

        const double start_score( score( *model));

        BCL_MessageStd( "before ccd score is " + util::Format()( start_score));

        util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > target_residue
        (
          new assemble::LocatorAA( 'A', 20)
        );

        fold::MutateProteinModelLoopDomainCCD loop_domain_ccd
        (
          loop_domain_locator,
          rng,
          math::Range< double>( 0.0, 1.0)
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_ccd( *model).GetArgument());
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_f.pdb")
        );
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );

        const double end_score( score( *new_model));
        const double expected_score( -0.72517);
        BCL_MessageStd( "after ccd score is " + util::Format()( end_score));
        BCL_ExampleCheckWithinTolerance( end_score, expected_score, 0.00001);
      }

      // test with LoopDomainCToN
      {
        RandomNumberGenerator rng( true, 0.7, 1);

        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // 21-23
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 21, 23), false));

        util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> > loop_domain_locator
        (
          new fold::LocatorLoopDomain( loop_segments, false)
        );

        const std::string input_pdb( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_c_to_n_21_23.pdb"));
        util::ShPtr< assemble::ProteinModel> model( Proteins::GetModel( input_pdb).Clone());

        util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > target_residue
        (
          new assemble::LocatorAA( 'A', 20)
        );

        storage::Set< biol::AtomType> superimpose_atom_types
        (
          storage::Set< biol::AtomType>::Create( biol::GetAtomTypes().N, biol::GetAtomTypes().CA, biol::GetAtomTypes().C)
        );

        util::ShPtr< assemble::ProteinModelData> pmd( model->GetProteinModelData()->Clone());
        util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > pmd_ldl_locators
        (
          new util::ShPtrList< fold::LocatorLoopDomain>
          (
            1, util::ShPtr< fold::LocatorLoopDomain>( new fold::LocatorLoopDomain( loop_segments, false))
          )
        );

        pmd->Insert( assemble::ProteinModelData::e_LoopDomainLocators, pmd_ldl_locators);
        model->SetProteinModelData( pmd);

        score::ProteinModelLoopDomainClosure score;

        const double start_score( score( *model));

        BCL_MessageStd( "before ccd score is " + util::Format()( start_score));

        fold::MutateProteinModelLoopDomainCCD loop_domain_ccd
        (
          loop_domain_locator,
          rng,
          math::Range< double>( 0.0, 1.0)
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_ccd( *model).GetArgument());
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_protein_model_loop_domain_ccd_g.pdb")
        );
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        const double end_score( score( *new_model));
        const double expected_score( -0.45761);

        BCL_MessageStd( "after ccd score is " + util::Format()( end_score));
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
        BCL_ExampleCheckWithinTolerance( end_score, expected_score, 0.00001);
      }

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

    void DoMutationSingleLoopSegment
    (
      const bool MUTATE_PHI, const double RANDOM_DOUBLE, const std::string &OUTPUT_PDB_FILENAME,
      const double PHI, const double PSI, const std::string &INPUT_PDB_FILENAME
    ) const
    {
      RandomNumberGenerator rng( MUTATE_PHI, RANDOM_DOUBLE, 1);

      // create list of loop segments
      storage::List< fold::LocatorLoopSegment> loop_segments;

      // 35-38
      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));

      // nterminal sse locator
      assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, INPUT_PDB_FILENAME)));

      util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> > sp_loop_domain_locator
      (
        new fold::LocatorLoopDomain( loop_segments)
      );

      util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > target_residue
      (
        new assemble::LocatorAA( 'A', 39)
      );

      fold::MutateProteinModelLoopDomainCCD loop_domain_ccd
      (
        sp_loop_domain_locator,
        rng,
        math::Range< double>( 0.0, 1.0)
      );

      util::ShPtr< assemble::ProteinModel> new_model( loop_domain_ccd( model).GetArgument());
      BCL_MessageDbg( "printing to file " + OUTPUT_PDB_FILENAME);
      Proteins::WriteModelToPDB( *new_model, OUTPUT_PDB_FILENAME);
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelLoopDomainCCD

  const ExampleClass::EnumType ExampleFoldMutateProteinModelLoopDomainCCD::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelLoopDomainCCD())
  );

} // namespace bcl
