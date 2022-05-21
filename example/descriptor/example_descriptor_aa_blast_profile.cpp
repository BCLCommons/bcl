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
#include "descriptor/bcl_descriptor_aa_blast_profile.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_blast_profile.cpp
  //!
  //! @author teixeipl
  //! @date Feb 6, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAABlastProfile :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAABlastProfile *Clone() const
    {
      return new ExampleDescriptorAABlastProfile( *this);
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
      // form the blast vector, blast list and blast profile
      const double blast_array[ 20] = { 0, -2, 0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2, 6, 0, -4, -3, -3};
      const linal::Vector< double> blast_vector( 20, blast_array);
      biol::BlastProfile blast_profile( blast_vector);
      linal::Vector< float> blast_description( blast_vector);

      // form an AA and set the blast profile
      biol::AA amino_acid( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));
      amino_acid.SetBlastProfile( blast_profile);
      biol::AASequence aa_sequence;
      aa_sequence.PushBack( amino_acid);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::AABlastProfile blast;

      // clone
      util::ShPtr< descriptor::AABlastProfile> sp_blast( blast.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( blast.GetSizeOfFeatures(), 20);

    ////////////////
    // operations //
    ////////////////

      // sequence to chain to protein model to protein model with cache
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( util::CloneToShPtr( aa_sequence)));
      assemble::ProteinModelWithCache protein_model_with_cache( assemble::ProteinModel( sp_chain), false);

      descriptor::Iterator< biol::AABase> itr( blast.GetType());
      itr.SetObject( protein_model_with_cache);
      blast.SetObject( protein_model_with_cache);

      // check the descriptions produced
      BCL_ExampleCheck( blast( itr), blast_description);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // set to ignore the terminal's line width; this way this example file won't change depending on the terminal's
      // line width.  This is necessary because currently most of the WriteHelp functions used by various classes do not
      // support taking a io::FixedLineWidthWriter as an argument
      util::GetLogger().SetIgnoreTerminalLineWidth( true);

      // write current descriptors to file
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile
      (
        write,
        GetExamples().GetExamplePath() + "/../../scripts/machine_learning/descriptors/aa/AADescriptorsHelp.txt"
      );
      write << "Numeric descriptors\n";
      util::Implementation< descriptor::Base< biol::AABase, float> >::ResetHaveDisplayedHelp();
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::ResetHaveDisplayedHelp();
      io::FixedLineWidthWriter writer;
      util::Implementation< descriptor::Base< biol::AABase, float> >::WriteInstancesHelp( writer);
      writer << "\nIdentification descriptors\n";
      util::Implementation< descriptor::Base< biol::AABase, char> >::ResetHaveDisplayedHelp();
      util::Implementation< descriptor::Base< biol::AABase, char> >::WriteInstancesHelp( writer);
      write << writer.String();
      io::File::CloseClearFStream( write);

      // reset fixed-width status
      util::GetLogger().SetIgnoreTerminalLineWidth( false);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAABlastProfile

  const ExampleClass::EnumType ExampleDescriptorAABlastProfile::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAABlastProfile())
  );

} // namespace bcl
