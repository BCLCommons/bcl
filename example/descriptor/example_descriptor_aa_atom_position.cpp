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
#include "descriptor/bcl_descriptor_aa_atom_position.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_atom_position.cpp
  //!
  //! @author mendenjl
  //! @date Apr 09, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAAtomPosition :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAAtomPosition *Clone() const
    {
      return new ExampleDescriptorAAAtomPosition( *this);
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
      // construct an implementation that will read proteins from the biol directory
      const std::string pdb_fasta_dir( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object to read proteins, as either pdbs or fastas, from the directory
      util::Implementation
      <
        assemble::RetrieveProteinModelWithCache
      > retriever_all( "ProteinDirectory( " + pdb_fasta_dir + ")");

      // retrieve 1b43A
      assemble::ProteinModelWithCache protein( *retriever_all->Retrieve( "1j27A"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////
      // create a descriptor that will retrieve the SSE size
      descriptor::AAAtomPosition tm_pos_retriever;
      // call set object with protein
      tm_pos_retriever.SetObject( protein);

      // clone
      util::ShPtr< descriptor::AAAtomPosition> sp_aa_info( tm_pos_retriever.Clone());

    /////////////////
    // data access //
    /////////////////

      // test  the GetLength function
      BCL_ExampleCheck( tm_pos_retriever.GetSizeOfFeatures(), 3);

    ///////////////
    // operators //
    ///////////////

      descriptor::Iterator< biol::AABase> itr( tm_pos_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        tm_pos_retriever( itr),
        linal::MakeVector< float>( 45.330, 9.448, 8.045)
      );
      itr.GotoPosition( 63);
      // CA position from atom 64
      BCL_ExampleCheck
      (
        tm_pos_retriever( itr),
        linal::MakeVector< float>( 43.023, 13.368, 19.029)
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

  }; //end ExampleDescriptorAAAtomPosition

  const ExampleClass::EnumType ExampleDescriptorAAAtomPosition::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAAtomPosition())
  );

} // namespace bcl
