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
#include "descriptor/bcl_descriptor_aa_info_from_sym_matrix_file.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_info_from_sym_matrix_file.cpp
  //!
  //! @author teixeipl, mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAInfoFromSymMatrixFile :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAInfoFromSymMatrixFile *Clone() const
    {
      return new ExampleDescriptorAAInfoFromSymMatrixFile( *this);
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
      // form an AA
      biol::AA amino_acid( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 1)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // retrieve a window of 2 characters around a string.  Use the reflecting instance so that we don't have to worry
      // about undefined characters
      util::Implementation< descriptor::Base< biol::AABase, float> > implementation( "AAInfoFromSymMatrixFile(suffix=.corr_mat_bcl)");

      // default constructor
      descriptor::AAInfoFromSymMatrixFile aa_info;

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( aa_info.GetSizeOfFeatures(), 1);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      assemble::ProteinModelWithCache protein_model_with_cache
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "mini.pdb")),
        true
      );

      descriptor::Iterator< biol::AABase> itr( aa_info.GetType());
      itr.SetObject( protein_model_with_cache);
      implementation->SetObject( protein_model_with_cache);

      // check the descriptions produced
      BCL_ExampleCheck( implementation->operator ()( itr)( 0), 1);
      ++itr;
      BCL_ExampleCheck( implementation->operator ()( itr)( 0), 2);
      ++itr;
      BCL_ExampleCheck( implementation->operator ()( itr)( 0), 4);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAAInfoFromSymMatrixFile

  const ExampleClass::EnumType ExampleDescriptorAAInfoFromSymMatrixFile::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAInfoFromSymMatrixFile())
  );

} // namespace bcl
