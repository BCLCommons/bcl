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
#include "fold/bcl_fold_mutate_loop_domain_dihedral.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_collector_loop_domain_random.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutate_protein_model_loop_domain_ccd.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_loop_domain_dihedral.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateLoopDomainDihedral :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateLoopDomainDihedral *Clone() const
    {
      return new ExampleFoldMutateLoopDomainDihedral( *this);
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

      double m_RandomDouble;
      uint64_t m_Random64BitInt;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      RandomNumberGenerator( const double RANDOM_DOUBLE, const size_t RANDOM_SIZE_T) :
        m_RandomDouble( RANDOM_DOUBLE),
        m_Random64BitInt( RANDOM_SIZE_T)
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
        return m_Random64BitInt;
      }

      //! @brief random double
      //! @return random number in double range
      double Double() const
      {
        return m_RandomDouble;
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstantPhiPsi
    //!
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class ConstantPhiPsi :
      public math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> >
    {

    private:

    //////////
    // data //
    //////////

      double m_Phi;
      double m_Psi;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ConstantPhiPsi() :
        m_Phi(),
        m_Psi()
      {
      }

      ConstantPhiPsi( const double PHI, const double PSI) :
        m_Phi( PHI),
        m_Psi( PSI)
      {
      }

      //! @brief Clone function
      //! @return pointer to new PhiPsiGeneratorRamachandran
      ConstantPhiPsi *Clone() const
      {
        return new ConstantPhiPsi( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a t_ResultType object
      //! @param RESIDUE Argument to be used to evaluate the function
      //! @return function value of the given argument
      storage::VectorND< 2, double> operator()( const fold::MutationResidue &RESIDUE) const
      {
        return storage::VectorND< 2, double>( m_Phi, m_Psi);
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

    }; // class ConstantPhiPsi

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test mutating the last residue in the n terminal anchor sse (i.e. the psi of the loops anchor residue)
      {
        BCL_MessageDbg( "test mutating the last residue in the n terminal anchor sse (i.e. the psi of the loops anchor residue)");
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateLoopDomainDihedral(), "mutate_loop_domain_dihedral_a.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( 0.1, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb", util::GetUndefinedDouble(), math::g_Pi / 2.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutation of first residue in loop domain
      {
        BCL_MessageDbg( "test mutation of first residue in loop domain");
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_b.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( 0.28, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb", math::g_Pi, math::g_Pi / 2.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutation of residue in interior (not first or last residue) of loop domain
      {
        BCL_MessageDbg( "test mutation of residue in interior (not first or last residue) of loop domain");
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_c.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( 0.4, out_pdb_name, -1.08035881, 2.27067336, "2LZM_rotated_35_38.pdb", math::g_Pi, math::g_Pi / 2.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutation of last residue of loop domain
      {
        BCL_MessageDbg( "test mutation of last residue of loop domain");
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_d.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( 0.8, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_35_38.pdb", math::g_Pi, math::g_Pi / 2.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutation of pseudo residue phi
      {
        BCL_MessageDbg( "test mutation of pseudo residue phi");
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_e.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationSingleLoopSegment( 0.95, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_35_38.pdb", math::g_Pi / 4.0, util::GetUndefinedDouble());
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test change psi of last residue in n terminal anchor sse (residue 11)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_f.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.05, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", util::GetUndefinedDouble(), math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate first residue in loop domain (residue 12)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_g.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.12, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate residue in loop domain (residue 13)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_h.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.20, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate residue in loop domain (residue 22)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_i.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.40, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate residue in loop domain (residue 35)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_j.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.55, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate last residue in loop domain (residue 38)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_k.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.90, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate pseudo residue in loop domain (residue 39)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_l.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.97, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 4.0, util::GetUndefinedDouble());
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate residue in loop domain (residue 21)
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_m.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.30, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      // test mutate residue in loop domain (residue 23) with negative rotation
      {
        const std::string out_pdb_name
        (
          AddExampleOutputPathToFilename( fold::MutateProteinModelLoopDomainCCD(), "mutate_loop_domain_dihedral_n.pdb")
        );
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        DoMutationMultipleLoopSegment( 0.50, out_pdb_name, -1.08035881, 0.0, "2LZM_rotated_12_38.pdb", -math::g_Pi / 2.0, math::g_Pi / 4.0);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
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

    void DoMutationSingleLoopSegment
    (
      const double RANDOM_DOUBLE, const std::string &OUTPUT_PDB_FILENAME, const double PSEUDO_PHI,
      const double PSEUDO_PSI,    const std::string &INPUT_PDB_FILENAME, const double GENERATOR_PHI,
      const double GENERATOR_PSI
    ) const
    {
      RandomNumberGenerator rng( RANDOM_DOUBLE, 1);

      util::ShPtr
      <
        find::CollectorInterface< storage::List< fold::MutationResidue>, fold::LoopDomain>
      > mutation_residue_collector( new fold::CollectorLoopDomainRandom( 1, rng));

      util::ShPtr
      <
        math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> >
      > phi_psi_generator( new ConstantPhiPsi( GENERATOR_PHI, GENERATOR_PSI));

      fold::MutateLoopDomainDihedral loop_mutate( mutation_residue_collector, phi_psi_generator);

      // create list of loop segments
      storage::List< fold::LocatorLoopSegment> loop_segments;

      // 35-38
      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));

      util::ShPtr< assemble::ProteinModel> model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, INPUT_PDB_FILENAME)).Clone());

      fold::LocatorLoopDomain loop_domain_locator( loop_segments);

      util::ShPtr< fold::LoopDomain> loop_domain( loop_domain_locator.Locate( *model));

      util::ShPtr< fold::LoopDomain> mutated_loop_domain( loop_mutate( *loop_domain).GetArgument());

      const util::ShPtr< assemble::ProteinModel> new_model( mutated_loop_domain->UpdateProteinModel( *model));

      BCL_MessageDbg( "printing to file " + OUTPUT_PDB_FILENAME);
      Proteins::WriteModelToPDB( *new_model, OUTPUT_PDB_FILENAME);
    }

    void DoMutationMultipleLoopSegment
    (
      const double RANDOM_DOUBLE, const std::string &OUTPUT_PDB_FILENAME, const double PSEUDO_PHI,
      const double PSEUDO_PSI,    const std::string &INPUT_PDB_FILENAME, const double GENERATOR_PHI,
      const double GENERATOR_PSI
    ) const
    {
      RandomNumberGenerator rng( RANDOM_DOUBLE, 1);

      util::ShPtr
      <
        find::CollectorInterface< storage::List< fold::MutationResidue>, fold::LoopDomain>
      > mutation_residue_collector( new fold::CollectorLoopDomainRandom( 1, rng));

      util::ShPtr
      <
        math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> >
      > phi_psi_generator( new ConstantPhiPsi( GENERATOR_PHI, GENERATOR_PSI));

      fold::MutateLoopDomainDihedral loop_mutate( mutation_residue_collector, phi_psi_generator);

      // create list of loop segments
      storage::List< fold::LocatorLoopSegment> loop_segments;

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 12, 13), false));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 14, 20), true));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 21, 23), false));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 24, 27), true));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 28, 30), true));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 31, 34), true));

      loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));

      util::ShPtr< assemble::ProteinModel> model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, INPUT_PDB_FILENAME)).Clone());

      fold::LocatorLoopDomain loop_domain_locator( loop_segments);

      util::ShPtr< fold::LoopDomain> loop_domain( loop_domain_locator.Locate( *model));

      util::ShPtr< fold::LoopDomain> mutated_loop_domain( loop_mutate( *loop_domain).GetArgument());

      const util::ShPtr< assemble::ProteinModel> new_model
      (
        mutated_loop_domain->UpdateProteinModel( *model)
      );

      BCL_MessageDbg( "printing to file " + OUTPUT_PDB_FILENAME);
      Proteins::WriteModelToPDB( *new_model, OUTPUT_PDB_FILENAME);
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateLoopDomainDihedral

  const ExampleClass::EnumType ExampleFoldMutateLoopDomainDihedral::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateLoopDomainDihedral())
  );

} // namespace bcl
