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
#include "fold/bcl_fold_mutate_protein_model_loop_domain_grow.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_loop_domain_grow.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelLoopDomainGrow :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelLoopDomainGrow *Clone() const
    {
      return new ExampleFoldMutateProteinModelLoopDomainGrow( *this);
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

    class ConstantPhiPsi :
      public math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> >
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

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

      //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
      //! @param RESIDUE Argument to be used to evaluate the function
      //! @return function value of the given argument
      storage::VectorND< 2, double> operator()( const fold::MutationResidue &RESIDUE) const
      {
        return storage::VectorND< 2, double>( math::g_Pi, math::g_Pi / 2.0);
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
      {
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // 51-55
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 51, 55), false));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments);

        assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));

        fold::MutateProteinModelLoopDomainGrow loop_domain_grow
        (
          util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> >( loop_domain_locator.Clone()),
          util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi())
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_grow( model).GetArgument());
        const std::string out_pdb_name( AddExampleOutputPathToFilename( loop_domain_grow, "mutate_protein_model_loop_domain_grow_a.pdb"));
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      {
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // add some segments
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 39, 50), true));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments);

        assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));
        BCL_MessageDbg( "growing b");
        fold::MutateProteinModelLoopDomainGrow loop_domain_grow
        (
          util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> >( loop_domain_locator.Clone()),
          util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi())
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_grow( model).GetArgument());
        const std::string out_pdb_name( AddExampleOutputPathToFilename( loop_domain_grow, "mutate_protein_model_loop_domain_grow_b.pdb"));
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      {
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // add some segments
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 39, 50), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 51, 55), true));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments);

        assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));
        BCL_MessageDbg( "growing c");
        fold::MutateProteinModelLoopDomainGrow loop_domain_grow
        (
          util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> >( loop_domain_locator.Clone()),
          util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi())
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_grow( model).GetArgument());
        const std::string out_pdb_name( AddExampleOutputPathToFilename( loop_domain_grow, "mutate_protein_model_loop_domain_grow_c.pdb"));
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      {
        BCL_MessageDbg( "growing d");
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // add some segments
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 39, 50), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 51, 55), false));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments);

        assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));
        fold::MutateProteinModelLoopDomainGrow loop_domain_grow
        (
          util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> >( loop_domain_locator.Clone()),
          util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi())
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_grow( model).GetArgument());
        const std::string out_pdb_name( AddExampleOutputPathToFilename( loop_domain_grow, "mutate_protein_model_loop_domain_grow_d.pdb"));
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

      {
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // add some segments
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 12, 13), false));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 14, 20), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 21, 23), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 24, 27), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 28, 30), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 31, 34), true));
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments);

        assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));
        BCL_MessageDbg( "growing e");
        fold::MutateProteinModelLoopDomainGrow loop_domain_grow
        (
          util::ShPtr< find::LocatorInterface< util::ShPtr< fold::LoopDomain>, assemble::DomainInterface> >( loop_domain_locator.Clone()),
          util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi())
        );

        util::ShPtr< assemble::ProteinModel> new_model( loop_domain_grow( model).GetArgument());
        const std::string out_pdb_name( AddExampleOutputPathToFilename( loop_domain_grow, "mutate_protein_model_loop_domain_grow_e.pdb"));
        const std::string correct_pdb_name( out_pdb_name + ".correct");
        BCL_MessageDbg( "printing to file " + out_pdb_name);
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
        pdb::Factory().WriteModelToPDB( *new_model, write);
        io::File::CloseClearFStream( write);
        BCL_ExampleCheck
        (
          io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
        );
      }

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

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelLoopDomainGrow

  const ExampleClass::EnumType ExampleFoldMutateProteinModelLoopDomainGrow::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelLoopDomainGrow())
  );

} // namespace bcl
