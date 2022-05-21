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
#include "descriptor/bcl_descriptor_mapped_sequence.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa_base.h"
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "descriptor/bcl_descriptor_iterator.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_mapped_sequence.cpp
  //!
  //! @author lib14, mendenjl
  //! @date Mar 14, 2015
  //! @remarks status complete
  //! @remarks reviewed by mendenjl
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   class ExampleDescriptorMappedSequence :
     public ExampleInterface
   {
       // single instance of this example class
       static const ExampleClass::EnumType s_Instance;

  public:

      ExampleDescriptorMappedSequence *Clone() const
      {
        return new ExampleDescriptorMappedSequence( *this);
      }

  /////////////////
  // data access //
  /////////////////

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      int Run() const
      {
        // check the alias of the mapped sequence
        descriptor::MappedSequence< float> mapped_sequence;

        util::Implementation< descriptor::Base< biol::AABase, float> > implementation
        (
          "MappedSequence("
          "file extension=.ncnv,"
          "key=AASeqID,"
          "size=2,"
          "delimiter=\",\""
          ")"
        );

        // ensure that each is defined
        BCL_ExampleCheck( implementation.IsDefined(), true);

        // check the descriptor string
        BCL_ExampleCheck( implementation.GetLabel().GetValue(), "MappedSequence");
        BCL_ExampleCheck( mapped_sequence.GetAlias(), "MappedSequence");

        // get the test protein
        assemble::ProteinModelWithCache protein_model_with_cache
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Descriptor, "2KSF.pdb")),
          true
        );

        descriptor::Iterator< biol::AABase> itr( implementation->GetType());
        itr.SetObject( protein_model_with_cache);
        implementation->SetObject( protein_model_with_cache);

        // check feature size
        BCL_ExampleCheck( implementation->GetSizeOfFeatures(), 2);

        // check ncnv values in the descriptor
        BCL_ExampleCheck( implementation->operator ()( itr), linal::MakeVector< float>( 1.76389,0.87383));
        ++itr;
        BCL_ExampleCheck( implementation->operator ()( itr), linal::MakeVector< float>( 4.19484,0.553131));
        ++itr;
        BCL_ExampleCheck( implementation->operator ()( itr), linal::MakeVector< float>( 7.8518,0.442586));

        return 0;
      } // Run
   };

     const ExampleClass::EnumType ExampleDescriptorMappedSequence::s_Instance
     (
       GetExamples().AddEnum( ExampleDescriptorMappedSequence())
     );

} // namespace bcl

