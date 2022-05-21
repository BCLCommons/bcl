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
#include "fold/bcl_fold_mutate_aa_sequence_grow.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_aa_sequence_grow.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateAASequenceGrow :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateAASequenceGrow *Clone() const
    {
      return new ExampleFoldMutateAASequenceGrow( *this);
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
      const std::string pdb_name( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb"));

      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_name, biol::GetAAClasses().e_AABackBone));

      assemble::LocatorSSE sse_to_grow_locator( 'A', 35, 38);
      util::SiPtr< const assemble::SSE> sse_to_grow( sse_to_grow_locator.Locate( protein_model));

      assemble::LocatorSSE nterminal_sse_locator( 'A', 31, 34);
      util::SiPtr< const assemble::SSE> nterminal_sse( nterminal_sse_locator.Locate( protein_model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      fold::MutateAASequenceGrow grow_mutate
      (
        util::ShPtr< math::FunctionInterfaceSerializable< fold::MutationResidue, storage::VectorND< 2, double> > >( new ConstantPhiPsi()),
        util::ToSiPtr( *nterminal_sse)
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      biol::AASequence grown_sequence( *grow_mutate( *sse_to_grow).GetArgument());

      pdb::Handler newpdb;
      size_t atom_number( 1);
      newpdb.AppendLines( pdb::Factory::WriteAASequenceToLines( grown_sequence, atom_number));
      const std::string out_pdb_name( AddExampleOutputPathToFilename( grow_mutate, "mutate_aa_sequence_grow.pdb"));
      const std::string correct_pdb_name( out_pdb_name + ".correct");
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, out_pdb_name);
      newpdb.WriteLines( write);
      io::File::CloseClearFStream( write);
      BCL_ExampleCheck
      (
        io::File::FilesMatch( out_pdb_name, correct_pdb_name), true
      );

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

  }; //end ExampleFoldMutateAASequenceGrow

  const ExampleClass::EnumType ExampleFoldMutateAASequenceGrow::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateAASequenceGrow())
  );

} // namespace bcl
